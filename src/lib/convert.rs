use std::collections::{HashMap, HashSet};
use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::{Context, bail};
use log::{info, warn};
use noodles::vcf::Header;
use noodles::vcf::header::record::value::Collection;
use noodles::vcf::header::record::value::map::format;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record_buf::info::field::Value as InfoField;
use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
use noodles::vcf::variant::record_buf::samples::Series as SampleSeries;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::record_buf::samples::sample::value::Array as SampleArray;
use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;

use crate::alleles::{is_skipped_alt, vep_alleles};
use crate::consequence::{VEP_RELEASE, rank, vep_major};
use crate::csq::{CsqEntry, CsqIndex, parse_entries};
use crate::fields::{Exclusions, FieldPlan, FieldSpec, Kind, Section, Shape, plan_fields};
use crate::genotype::copies_of_alt;
use crate::input::{VariantSource, open_variants};
use crate::matrix::{Cell, SparseBuilder};
use crate::output::{Column, Output, UnsValue, X_ENCODING};
use crate::pick::{
    AltAllele, Genes, Hashed, MINIMAL_OUTPUT, ManeTable, PickContext, Picked, entries_for_allele,
    missing_tie_fields, pick, read_genes, read_mane_table, worst_panel_impact,
};
use crate::var_table::{InfoValue, RecordRow, VarColumns, VarTable, widen};

const PICK_RULE: &str = "1. Keep this ALT's CSQ entries. Evidence of a split after VEP stops the \
    run. A sequence ALT keeps the entries whose Allele is VEP's spelling of the ALT. A breakend \
    keeps the entries that spell it verbatim. Otherwise a symbolic ALT or breakend keeps its \
    ALLELE_NUM entries, or the entries with VEP's spelling when ALLELE_NUM is not declared, and \
    a breakend with none of its own keeps its reference-side entries. \
    2. Keep the entries of the genes the variant lies in, preferring the --genes panel; \
    when it lies in no gene, keep every entry. \
    3. Keep the MANE Select and MANE Plus Clinical entries if any, \
    else the CANONICAL entries, else all of them. \
    4. Take the entry with the most severe Consequence term in VEP's order. \
    5. Break ties by MANE Select, then CANONICAL, then CCDS, then SWISSPROT, \
    then the longest protein, then the lowest transcript ID.";

pub struct Config {
    pub input: PathBuf,
    pub genes: Option<PathBuf>,
    pub mane_table: Option<PathBuf>,
    pub exclusions: Exclusions,
    pub threads: NonZeroUsize,
    pub command_line: String,
}

pub fn convert(config: &Config) -> anyhow::Result<Output> {
    check_required_fields(&config.exclusions)?;
    let genes = config.genes.as_deref().map(read_genes).transpose()?;
    let mane = config
        .mane_table
        .as_deref()
        .map(read_mane_table)
        .transpose()?;

    let VariantSource {
        header,
        header_lines,
        samples,
        records,
    } = open_variants(&config.input, config.threads)?;
    let csq = read_csq(&header, &config.exclusions)?;
    check_header(&header, &header_lines, csq.as_ref(), &config.exclusions)?;
    let plan = plan_fields(&header, &config.exclusions);
    check_layer_names(&plan)?;
    let uns = provenance(
        config,
        &header,
        csq.as_ref(),
        &plan,
        genes.as_ref(),
        mane.as_ref(),
    );

    let lookups = Lookups {
        csq: csq.as_ref(),
        genes: genes.as_ref().map(|genes| &genes.value),
        mane: mane.as_ref().map(|mane| &mane.value),
    };
    let mut conversion = Conversion::new(&plan, lookups, &config.exclusions, samples.len());
    if let Some(name) = conversion.var.duplicate_column() {
        bail!("two fields make the var column {name}");
    }
    warn_about_header(&header, csq.as_ref());
    if csq.is_none() {
        note_no_csq(config, &header);
    }
    if mane.is_some()
        && csq.as_ref().is_some_and(|index| {
            index.position("MANE_SELECT").is_some()
                || index.position("MANE_PLUS_CLINICAL").is_some()
        })
    {
        warn!("--mane-table is ignored: the VCF declares MANE_SELECT or MANE_PLUS_CLINICAL");
    }

    for record in records {
        let record = record?;
        conversion.push_record(&record).with_context(|| {
            let position = record.variant_start().map_or(0, usize::from);
            format!("{}:{position}", record.reference_sequence_name())
        })?;
    }

    let Conversion {
        var,
        matrix,
        n_sites_covered,
        n_sites_alt,
        skipped_alts,
        structural_without_entry,
        records_without_alt,
        ..
    } = conversion;
    let VarColumns {
        names: var_names,
        columns: var,
    } = var.finish();
    let matrices = matrix.finish();
    for (field, count) in &matrices.missing {
        warn!("{field}: missing in {count} of the cells with reads");
    }
    warn_about_counts(skipped_alts, structural_without_entry, records_without_alt);

    Ok(Output {
        obs_names: samples,
        obs: vec![
            ("n_sites_covered".to_string(), Column::Int(n_sites_covered)),
            ("n_sites_alt".to_string(), Column::Int(n_sites_alt)),
        ],
        var_names,
        var,
        x: matrices.x,
        layers: matrices.layers,
        uns,
    })
}

fn warn_about_counts(skipped_alts: u64, structural_without_entry: u64, records_without_alt: u64) {
    if skipped_alts > 0 {
        warn!("*, <*> and <NON_REF> ALT alleles skipped: {skipped_alts}");
    }
    if structural_without_entry > 0 {
        warn!("symbolic or breakend ALT alleles with no CSQ entry: {structural_without_entry}");
    }
    if records_without_alt > 0 {
        warn!("records with no ALT skipped: {records_without_alt}");
    }
}

fn check_required_fields(exclusions: &Exclusions) -> anyhow::Result<()> {
    if exclusions.contains(Section::Format, "AD") {
        bail!("FORMAT/AD is required");
    }
    if exclusions.contains(Section::Format, "GT") {
        bail!("FORMAT/GT cannot be excluded: X comes from it");
    }
    Ok(())
}

fn read_csq(header: &Header, exclusions: &Exclusions) -> anyhow::Result<Option<CsqIndex>> {
    if exclusions.contains(Section::Info, "CSQ") {
        return Ok(None);
    }
    let csq = CsqIndex::from_header(header);
    if csq.is_none() && header.infos().contains_key("CSQ") {
        bail!("INFO/CSQ has no \"Format: \" list in its Description");
    }
    Ok(csq)
}

fn check_header(
    header: &Header,
    header_lines: &[String],
    csq: Option<&CsqIndex>,
    exclusions: &Exclusions,
) -> anyhow::Result<()> {
    check_not_minimal(header, csq)?;
    if csq.is_some() && split_after_vep(header_lines) {
        bail!(
            "the header shows bcftools norm splitting multiallelics after VEP; annotate after splitting"
        );
    }
    if let Some(name) = exclusions.undeclared(header, csq).first() {
        bail!("--exclude-field {name} names nothing in this VCF's header");
    }
    let Some(ad) = header.formats().get("AD") else {
        bail!("the VCF declares no FORMAT/AD; with bcftools mpileup, add --annotate FORMAT/AD");
    };
    if ad.number() != format::Number::ReferenceAlternateBases || ad.ty() != format::Type::Integer {
        bail!("FORMAT/AD must be declared Number=R, Type=Integer");
    }
    Ok(())
}

fn names_minimal(command_line: &str) -> bool {
    command_line
        .split_whitespace()
        .any(|token| token.trim_matches(['\'', '"']) == "--minimal")
}

fn splits_multiallelics(norm_command: &str) -> bool {
    let mut tokens = norm_command.split_whitespace().peekable();
    while let Some(token) = tokens.next() {
        let value = match token {
            "-m" | "--multiallelics" => tokens.peek().copied(),
            _ => token
                .strip_prefix("--multiallelics=")
                .or_else(|| token.strip_prefix("-m")),
        };
        if value.is_some_and(|value| value.starts_with('-')) {
            return true;
        }
    }
    false
}

fn split_after_vep(header_lines: &[String]) -> bool {
    header_lines
        .iter()
        .skip_while(|line| !line.starts_with("##VEP-command-line="))
        .filter_map(|line| line.strip_prefix("##bcftools_normCommand="))
        .any(splits_multiallelics)
}

fn check_not_minimal(header: &Header, csq: Option<&CsqIndex>) -> anyhow::Result<()> {
    let Some(index) = csq else {
        return Ok(());
    };
    let on_the_command_line = vep_header(header)
        .iter()
        .filter_map(|line| line.strip_prefix("##VEP-command-line="))
        .any(names_minimal);
    if index.position("MINIMISED").is_some() || on_the_command_line {
        bail!(MINIMAL_OUTPUT);
    }
    Ok(())
}

fn check_layer_names(plan: &FieldPlan) -> anyhow::Result<()> {
    let mut sources: HashMap<String, &str> = HashMap::new();
    for spec in &plan.format {
        for name in spec.column_names() {
            if let Some(other) = sources.insert(name.clone(), &spec.id) {
                bail!(
                    "FORMAT/{other} and FORMAT/{} both make the layer {name}",
                    spec.id
                );
            }
        }
    }
    Ok(())
}

fn warn_about_header(header: &Header, csq: Option<&CsqIndex>) {
    if !header.formats().contains_key("GT") {
        warn!("the VCF declares no FORMAT/GT, so X is empty");
    }
    let missing = csq.map(missing_tie_fields).unwrap_or_default();
    if !missing.is_empty() {
        warn!("the CSQ header declares no {}", missing.join(", "));
    }
    if csq.is_some() {
        warn_about_vep_version(header);
    }
}

fn warn_about_vep_version(header: &Header) {
    let other = vep_header(header)
        .iter()
        .filter_map(|line| vep_major(line.strip_prefix("##VEP=")?))
        .find(|major| *major != VEP_RELEASE);
    if let Some(major) = other {
        warn!(
            "annotated by VEP {major}, but CSQ types and consequence ranks are pinned to \
            VEP {VEP_RELEASE}"
        );
    }
}

fn note_no_csq(config: &Config, header: &Header) {
    info!("no INFO/CSQ to read, so the variants have no annotation");
    let flags = [
        (config.genes.is_some(), "--genes".to_string()),
        (config.mane_table.is_some(), "--mane-table".to_string()),
    ];
    let exclusions = config.exclusions.without_csq_to_act_on(header);
    let unused = flags
        .into_iter()
        .filter_map(|(given, flag)| given.then_some(flag))
        .chain(
            exclusions
                .iter()
                .map(|name| format!("--exclude-field {name}")),
        );
    for flag in unused {
        warn!("{flag} does nothing: there is no INFO/CSQ to act on");
    }
}

fn provenance(
    config: &Config,
    header: &Header,
    csq: Option<&CsqIndex>,
    plan: &FieldPlan,
    genes: Option<&Hashed<Genes>>,
    mane: Option<&Hashed<ManeTable>>,
) -> Vec<(String, UnsValue)> {
    let text = |name: &str, value: String| (name.to_string(), UnsValue::Text(value));
    let list = |name: &str, values: Vec<String>| (name.to_string(), UnsValue::List(values));

    let source_vcf = config.input.file_name().map_or_else(
        || config.input.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    );
    let mut uns = vec![
        text("annatar_version", env!("CARGO_PKG_VERSION").to_string()),
        text("command_line", config.command_line.clone()),
        text("source_vcf", source_vcf),
        list("vep_header", vep_header(header)),
        list(
            "csq_fields",
            csq.map(|index| index.names.clone()).unwrap_or_default(),
        ),
        text("pick_rule", PICK_RULE.to_string()),
        text("x_encoding", X_ENCODING.to_string()),
        list("skipped_fields", plan.skipped.clone()),
    ];
    if let Some(genes) = genes {
        uns.push(list("panel_genes", genes.value.sorted()));
        uns.push(text("panel_genes_sha256", genes.sha256.clone()));
    }
    if let Some(mane) = mane {
        uns.push(text("mane_table_sha256", mane.sha256.clone()));
    }
    uns
}

fn vep_header(header: &Header) -> Vec<String> {
    header
        .other_records()
        .iter()
        .filter(|(key, _)| matches!(key.as_ref(), "VEP" | "VEP-command-line"))
        .flat_map(|(key, records)| match records {
            Collection::Unstructured(lines) => {
                lines.iter().map(|line| format!("##{key}={line}")).collect()
            }
            Collection::Structured(_) => Vec::new(),
        })
        .collect()
}

fn expected_len(shape: Shape, n_alts: usize) -> usize {
    match shape {
        Shape::Scalar => 1,
        Shape::PerAlt => n_alts,
        Shape::RefAlt => n_alts + 1,
        Shape::Fixed(n) => n,
    }
}

fn positions(shape: Shape, alt_position: usize) -> Vec<usize> {
    match shape {
        Shape::Scalar => vec![0],
        Shape::PerAlt => vec![alt_position],
        Shape::RefAlt => vec![0, alt_position + 1],
        Shape::Fixed(n) => (0..n).collect(),
    }
}

fn check_len(section: &str, id: &str, actual: usize, expected: usize) -> anyhow::Result<()> {
    if actual != expected {
        bail!("{section}/{id} has {actual} values but this record needs {expected}");
    }
    Ok(())
}

fn info_elements(value: &InfoField) -> Vec<InfoValue> {
    match value {
        InfoField::Integer(n) => vec![InfoValue::Int(Some(i64::from(*n)))],
        InfoField::Float(x) => vec![InfoValue::Float(Some(widen(*x)))],
        InfoField::Flag => vec![InfoValue::Flag(true)],
        InfoField::Character(c) => vec![InfoValue::Text(Some(c.to_string()))],
        InfoField::String(s) => vec![InfoValue::Text(Some(s.clone()))],
        InfoField::Array(InfoArray::Integer(values)) => values
            .iter()
            .map(|n| InfoValue::Int(n.map(i64::from)))
            .collect(),
        InfoField::Array(InfoArray::Float(values)) => values
            .iter()
            .map(|x| InfoValue::Float(x.map(widen)))
            .collect(),
        InfoField::Array(InfoArray::Character(values)) => values
            .iter()
            .map(|c| InfoValue::Text(c.map(|c| c.to_string())))
            .collect(),
        InfoField::Array(InfoArray::String(values)) => {
            values.iter().cloned().map(InfoValue::Text).collect()
        }
    }
}

fn read_info(
    record: &RecordBuf,
    specs: &[FieldSpec],
    n_alts: usize,
) -> anyhow::Result<Vec<Vec<InfoValue>>> {
    specs
        .iter()
        .map(|spec| {
            let Some(Some(value)) = record.info().get(spec.id.as_str()) else {
                return Ok(Vec::new());
            };
            let values = info_elements(value);
            check_len(
                "INFO",
                &spec.id,
                values.len(),
                expected_len(spec.shape, n_alts),
            )?;
            Ok(values)
        })
        .collect()
}

fn info_row(specs: &[FieldSpec], values: &[Vec<InfoValue>], alt_position: usize) -> Vec<InfoValue> {
    let missing = |kind: Kind| match kind {
        Kind::Integer => InfoValue::Int(None),
        Kind::Float => InfoValue::Float(None),
        Kind::Flag => InfoValue::Flag(false),
        Kind::Text => InfoValue::Text(None),
    };
    specs
        .iter()
        .zip(values)
        .flat_map(|(spec, values)| {
            positions(spec.shape, alt_position)
                .into_iter()
                .map(move |position| {
                    values
                        .get(position)
                        .cloned()
                        .unwrap_or_else(|| missing(spec.kind))
                })
        })
        .collect()
}

fn raw_csq(record: &RecordBuf) -> Vec<String> {
    match record.info().get("CSQ").flatten() {
        Some(InfoField::Array(InfoArray::String(entries))) => {
            entries.iter().flatten().cloned().collect()
        }
        Some(InfoField::String(entries)) => entries.split(',').map(str::to_string).collect(),
        _ => Vec::new(),
    }
}

fn joined<'a>(values: impl IntoIterator<Item = &'a String>) -> Option<String> {
    let values: Vec<&str> = values.into_iter().map(String::as_str).collect();
    (!values.is_empty()).then(|| values.join(";"))
}

struct RecordSamples<'r> {
    gt: Option<SampleSeries<'r>>,
    ad: Option<SampleSeries<'r>>,
    layers: Vec<Option<SampleSeries<'r>>>,
}

impl<'r> RecordSamples<'r> {
    fn read(
        record: &'r RecordBuf,
        specs: &[FieldSpec],
        n_obs: usize,
        n_alts: usize,
    ) -> anyhow::Result<Self> {
        let samples = record.samples();
        let layers: Vec<_> = specs.iter().map(|spec| samples.select(&spec.id)).collect();
        for (spec, series) in specs.iter().zip(&layers) {
            let expected = expected_len(spec.shape, n_alts);
            for sample in 0..n_obs {
                if let Some(value) = value_at(series.as_ref(), sample) {
                    check_len("FORMAT", &spec.id, sample_len(value), expected)?;
                }
            }
        }
        let gt = samples.select("GT");
        let ad = samples.select("AD");
        for sample in 0..n_obs {
            check_genotype(value_at(gt.as_ref(), sample), n_alts)?;
            check_depths(value_at(ad.as_ref(), sample))?;
        }
        Ok(Self { gt, ad, layers })
    }

    fn fill_alleles(&self, sample: usize, alleles: &mut Vec<Option<usize>>) {
        alleles.clear();
        if let Some(SampleValue::Genotype(genotype)) = value_at(self.gt.as_ref(), sample) {
            alleles.extend(genotype.as_ref().iter().map(Allele::position));
        }
    }

    fn allelic_depths(&self, sample: usize) -> &[Option<i32>] {
        match value_at(self.ad.as_ref(), sample) {
            Some(SampleValue::Array(SampleArray::Integer(depths))) => depths,
            _ => &[],
        }
    }

    fn cell(&self, layer: usize, sample: usize, position: usize, kind: Kind) -> Cell {
        let value = value_at(self.layers[layer].as_ref(), sample);
        match kind {
            Kind::Float => Cell::Float(float_at(value, position)),
            _ => Cell::Int(int_at(value, position)),
        }
    }
}

fn check_genotype(value: Option<&SampleValue>, n_alts: usize) -> anyhow::Result<()> {
    if let Some(SampleValue::Genotype(genotype)) = value
        && let Some(allele) = genotype
            .as_ref()
            .iter()
            .filter_map(Allele::position)
            .find(|allele| *allele > n_alts)
    {
        bail!("FORMAT/GT calls allele {allele} but the last allele here is {n_alts}");
    }
    Ok(())
}

fn check_depths(value: Option<&SampleValue>) -> anyhow::Result<()> {
    if let Some(SampleValue::Array(SampleArray::Integer(depths))) = value
        && let Some(count) = depths.iter().flatten().find(|count| **count < 0)
    {
        bail!("FORMAT/AD holds a negative count, {count}");
    }
    Ok(())
}

fn value_at<'a>(series: Option<&'a SampleSeries<'_>>, sample: usize) -> Option<&'a SampleValue> {
    series?.get(sample).flatten()
}

fn sample_len(value: &SampleValue) -> usize {
    match value {
        SampleValue::Array(SampleArray::Integer(values)) => values.len(),
        SampleValue::Array(SampleArray::Float(values)) => values.len(),
        SampleValue::Array(SampleArray::Character(values)) => values.len(),
        SampleValue::Array(SampleArray::String(values)) => values.len(),
        _ => 1,
    }
}

fn int_at(value: Option<&SampleValue>, position: usize) -> Option<i32> {
    match value? {
        SampleValue::Integer(n) => Some(*n),
        SampleValue::Array(SampleArray::Integer(values)) => values.get(position).copied().flatten(),
        _ => None,
    }
}

fn float_at(value: Option<&SampleValue>, position: usize) -> Option<f32> {
    let value = match value? {
        SampleValue::Float(x) => Some(*x),
        SampleValue::Array(SampleArray::Float(values)) => values.get(position).copied().flatten(),
        _ => None,
    };
    value.filter(|x| !x.is_nan())
}

fn reads_at(depths: &[Option<i32>], position: usize) -> u64 {
    let reads = depths.get(position).copied().flatten().unwrap_or(0);
    u64::from(reads.unsigned_abs())
}

#[derive(Default)]
struct ColumnTotals {
    n_cells_covered: u32,
    n_cells_alt: u32,
    ref_reads: u64,
    alt_reads: u64,
}

#[derive(Clone, Copy)]
struct Lookups<'a> {
    csq: Option<&'a CsqIndex>,
    genes: Option<&'a Genes>,
    mane: Option<&'a ManeTable>,
}

#[derive(Default)]
struct AltPick {
    picked: Option<Picked>,
    worst: Option<&'static str>,
    without_entry: bool,
}

impl Lookups<'_> {
    fn pick(&self, entries: &[CsqEntry], allele: &AltAllele) -> anyhow::Result<AltPick> {
        let index = match self.csq {
            Some(index) if !entries.is_empty() => index,
            _ => return Ok(AltPick::default()),
        };
        let candidates = entries_for_allele(entries, index, allele)?;
        if candidates.is_empty() {
            return Ok(AltPick {
                without_entry: true,
                ..AltPick::default()
            });
        }
        let ctx = PickContext {
            index,
            genes: self.genes,
            mane: self.mane,
        };
        Ok(AltPick {
            picked: pick(entries, &candidates, &ctx),
            worst: worst_panel_impact(entries, &candidates, &ctx),
            without_entry: false,
        })
    }
}

struct Conversion<'a> {
    plan: &'a FieldPlan,
    lookups: Lookups<'a>,
    var: VarTable,
    matrix: SparseBuilder,
    n_sites_covered: Vec<i64>,
    n_sites_alt: Vec<i64>,
    skipped_alts: u64,
    structural_without_entry: u64,
    records_without_alt: u64,
    unranked_terms: HashSet<String>,
}

fn svtype(record: &RecordBuf) -> Option<&str> {
    match record.info().get("SVTYPE").flatten() {
        Some(InfoField::String(svtype)) => Some(svtype),
        _ => None,
    }
}

impl<'a> Conversion<'a> {
    fn new(
        plan: &'a FieldPlan,
        lookups: Lookups<'a>,
        exclusions: &Exclusions,
        n_obs: usize,
    ) -> Self {
        let layers: Vec<(String, Kind)> = plan
            .format
            .iter()
            .flat_map(|spec| {
                spec.column_names()
                    .into_iter()
                    .map(|name| (name, spec.kind))
            })
            .collect();
        Self {
            plan,
            lookups,
            var: VarTable::new(
                &plan.info,
                lookups.csq,
                exclusions,
                lookups.csq.is_some() && lookups.genes.is_some(),
            ),
            matrix: SparseBuilder::new(n_obs, &layers),
            n_sites_covered: vec![0; n_obs],
            n_sites_alt: vec![0; n_obs],
            skipped_alts: 0,
            structural_without_entry: 0,
            records_without_alt: 0,
            unranked_terms: HashSet::new(),
        }
    }

    fn push_record(&mut self, record: &RecordBuf) -> anyhow::Result<()> {
        let alts = record.alternate_bases().as_ref();
        let vep = vep_alleles(record.reference_bases(), alts, svtype(record));
        let entries = match self.lookups.csq {
            Some(index) => parse_entries(&raw_csq(record), index)?,
            None => Vec::new(),
        };
        self.warn_about_unranked_terms(&entries);
        let info = read_info(record, &self.plan.info, alts.len())?;
        let n_obs = self.n_sites_covered.len();
        let samples = RecordSamples::read(record, &self.plan.format, n_obs, alts.len())?;
        if alts.is_empty() {
            self.records_without_alt += 1;
        }

        for (position, (alt, spelling)) in alts.iter().zip(&vep).enumerate() {
            if is_skipped_alt(alt) {
                self.skipped_alts += 1;
                continue;
            }
            let allele = AltAllele {
                alt,
                spelling,
                number: position + 1,
                alts,
            };
            let AltPick {
                picked,
                worst,
                without_entry,
            } = self.lookups.pick(&entries, &allele)?;
            self.structural_without_entry += u64::from(without_entry);
            let totals = self.push_column(&samples, position);
            let row = RecordRow {
                chrom: record.reference_sequence_name().to_string(),
                pos: record.variant_start().map_or(0, usize::from),
                reference: record.reference_bases().to_string(),
                alt: alt.clone(),
                id: joined(record.ids().as_ref()),
                qual: record.quality_score(),
                filter: joined(record.filters().as_ref()),
                n_cells_covered: totals.n_cells_covered,
                n_cells_alt: totals.n_cells_alt,
                ref_reads: totals.ref_reads,
                alt_reads: totals.alt_reads,
            };
            let picked = picked.map(|picked| (&entries[picked.entry], picked.tier));
            let values = info_row(&self.plan.info, &info, position);
            self.var.push(row, values, picked, worst)?;
        }
        Ok(())
    }

    fn warn_about_unranked_terms(&mut self, entries: &[CsqEntry]) {
        let Some(index) = self.lookups.csq else {
            return;
        };
        let terms = entries
            .iter()
            .filter_map(|entry| entry.get(index, "Consequence"))
            .flat_map(|consequence| consequence.split('&'));
        for term in terms {
            if rank(term).is_none() && self.unranked_terms.insert(term.to_string()) {
                warn!(
                    "VEP {VEP_RELEASE} does not rank the Consequence term {term}, so it sorts last"
                );
            }
        }
    }

    fn push_column(&mut self, samples: &RecordSamples, alt_position: usize) -> ColumnTotals {
        let specs = &self.plan.format;
        let slots: Vec<(usize, usize, Kind)> = specs
            .iter()
            .enumerate()
            .flat_map(|(layer, spec)| {
                positions(spec.shape, alt_position)
                    .into_iter()
                    .map(move |position| (layer, position, spec.kind))
            })
            .collect();
        let allele = alt_position + 1;
        let n_obs = self.n_sites_covered.len();
        let mut alleles = Vec::new();
        let mut totals = ColumnTotals::default();

        for (row, sample) in (0u32..).zip(0..n_obs) {
            samples.fill_alleles(sample, &mut alleles);
            self.matrix
                .push_copies(row, copies_of_alt(&alleles, allele));
            let depths = samples.allelic_depths(sample);
            let has_reads = depths.iter().flatten().any(|reads| *reads > 0);
            for (slot, &(layer, position, kind)) in slots.iter().enumerate() {
                let cell = samples.cell(layer, sample, position, kind);
                self.matrix.push_cell(slot, row, cell, has_reads);
            }
            let ref_reads = reads_at(depths, 0);
            let alt_reads = reads_at(depths, allele);
            if ref_reads + alt_reads == 0 {
                continue;
            }

            totals.n_cells_covered += 1;
            totals.ref_reads += ref_reads;
            totals.alt_reads += alt_reads;
            self.n_sites_covered[sample] += 1;
            if alt_reads > 0 {
                totals.n_cells_alt += 1;
                self.n_sites_alt[sample] += 1;
            }
        }
        self.matrix.end_column();
        totals
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case("norm -m -any in.vcf; Date=Wed Nov", true)]
    #[case("norm -m-any in.vcf; Date=Wed Nov", true)]
    #[case("norm -m -both in.vcf", true)]
    #[case("norm -m- in.vcf", true)]
    #[case("norm --multiallelics -any in.vcf", true)]
    #[case("norm --multiallelics=-both in.vcf", true)]
    #[case("norm in.vcf -m -any; Date=Wed Nov", true)]
    #[case("norm -m +any in.vcf", false)]
    #[case("norm -m+any in.vcf", false)]
    #[case("norm --multiallelics=+both in.vcf", false)]
    #[case("norm -d exact in.vcf", false)]
    #[case("norm -f ref.fa -o split-m-any.vcf in.vcf", false)]
    #[case("norm --multi-overlaps . in.vcf", false)]
    fn a_norm_command_splits_only_with_a_minus(#[case] command: &str, #[case] expected: bool) {
        assert_eq!(splits_multiallelics(command), expected);
    }

    fn lines(lines: &[&str]) -> Vec<String> {
        lines.iter().map(ToString::to_string).collect()
    }

    #[rstest]
    #[case(&["##VEP-command-line='vep'", "##bcftools_normCommand=norm -m -any x"], true)]
    #[case(
        &["##bcftools_normCommand=norm -f r", "##VEP-command-line='vep'", "##bcftools_normCommand=norm -m-any x"],
        true
    )]
    #[case(&["##bcftools_normCommand=norm -m -any x", "##VEP-command-line='vep'"], false)]
    #[case(&["##VEP-command-line='vep'", "##bcftools_normCommand=norm -d exact x"], false)]
    #[case(&["##VEP-command-line='vep'", "##bcftools_viewCommand=view -m -any x"], false)]
    #[case(&["##bcftools_normCommand=norm -m -any x"], false)]
    #[case(&[], false)]
    fn only_a_split_after_vep_counts(#[case] header: &[&str], #[case] expected: bool) {
        assert_eq!(split_after_vep(&lines(header)), expected);
    }

    #[rstest]
    #[case("'vep --cache --minimal --offline'", true)]
    #[case("'vep --cache --minimal'", true)]
    #[case("'--minimal --cache'", true)]
    #[case("\"vep --minimal\"", true)]
    #[case("vep --minimal", true)]
    #[case("'vep --cache --custom [PATH]/x--minimal.vcf,a,vcf'", false)]
    #[case("'vep --cache --plugin [PATH]/--minimal'", false)]
    #[case("'vep --minimal_output'", false)]
    #[case("'vep --cache'", false)]
    fn minimal_is_a_whole_token(#[case] command_line: &str, #[case] expected: bool) {
        assert_eq!(names_minimal(command_line), expected);
    }

    #[test]
    fn a_nan_float_reads_as_missing() {
        assert_eq!(float_at(Some(&SampleValue::Float(f32::NAN)), 0), None);
        assert_eq!(float_at(Some(&SampleValue::Float(0.5)), 0), Some(0.5));
        let array = SampleValue::Array(SampleArray::Float(vec![Some(f32::NAN), Some(0.5), None]));
        assert_eq!(float_at(Some(&array), 0), None);
        assert_eq!(float_at(Some(&array), 1), Some(0.5));
    }
}
