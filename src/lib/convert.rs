use std::collections::HashMap;
use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use log::warn;
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
use polars::prelude::{DataFrame, NamedFrom, Series};

use crate::alleles::{is_skipped_alt, vep_alleles};
use crate::csq::{CsqEntry, CsqIndex, parse_entries};
use crate::fields::{Exclusions, FieldPlan, FieldSpec, Kind, Section, Shape, plan_fields};
use crate::genotype::{covered, genotype_code};
use crate::h5ad::{LAYER_ENCODING, Output, UnsValue, X_ENCODING, write_h5ad};
use crate::input::{VariantSource, open_variants};
use crate::matrix::{Cell, SparseBuilder};
use crate::pick::{
    Genes, Hashed, ManeTable, PickContext, Picked, entries_for_allele, missing_tie_fields, pick,
    read_genes, read_mane_table, worst_panel_impact,
};
use crate::var_table::{InfoValue, RecordRow, VarTable, widen};

const PICK_RULE: &str = "1. Keep the CSQ entries whose Allele is VEP's allele string for this ALT. \
    2. Keep the entries of the genes the variant lies in, preferring the --genes panel; \
    when it lies in no gene, keep every entry. \
    3. Keep the MANE Select and MANE Plus Clinical entries if any, \
    else the CANONICAL entries, else all of them. \
    4. Take the entry with the most severe Consequence term in VEP's order. \
    5. Break ties by MANE Select, then CANONICAL, then CCDS, then SWISSPROT, \
    then the longest protein, then the lowest transcript ID.";

pub struct Config {
    pub input: PathBuf,
    pub output: PathBuf,
    pub genes: Option<PathBuf>,
    pub mane_table: Option<PathBuf>,
    pub exclusions: Exclusions,
    pub threads: NonZeroUsize,
    pub command_line: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Summary {
    pub cells: usize,
    pub variants: usize,
    pub skipped_alts: u64,
}

pub fn convert(config: &Config) -> anyhow::Result<Summary> {
    check_required_fields(&config.exclusions)?;
    let genes = config.genes.as_deref().map(read_genes).transpose()?;
    let mane = config
        .mane_table
        .as_deref()
        .map(read_mane_table)
        .transpose()?;

    let VariantSource {
        header,
        samples,
        records,
    } = open_variants(&config.input, config.threads)?;
    let csq = CsqIndex::from_header(&header);
    check_header(&header, csq.as_ref(), &config.exclusions)?;
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
        records_without_alt,
        ..
    } = conversion;
    let (var_names, var) = var.finish()?;
    let matrices = matrix.finish();
    for (field, count) in &matrices.missing {
        warn!("{field}: missing in {count} of the covered cells");
    }
    if skipped_alts > 0 {
        warn!("symbolic ALT alleles skipped: {skipped_alts}");
    }
    if records_without_alt > 0 {
        warn!("records with no ALT skipped: {records_without_alt}");
    }

    let summary = Summary {
        cells: samples.len(),
        variants: var_names.len(),
        skipped_alts,
    };
    let obs = DataFrame::new_infer_height(vec![
        Series::new("n_sites_covered".into(), n_sites_covered).into(),
        Series::new("n_sites_alt".into(), n_sites_alt).into(),
    ])?;
    let output = Output {
        obs_names: samples,
        obs,
        var_names,
        var,
        x: matrices.x,
        layers: matrices.layers,
        uns,
    };
    write_output(&config.output, output)?;
    Ok(summary)
}

fn check_required_fields(exclusions: &Exclusions) -> anyhow::Result<()> {
    if exclusions.contains(Section::Format, "AD") {
        bail!("FORMAT/AD is required");
    }
    if exclusions.contains(Section::Format, "GT") {
        bail!("FORMAT/GT cannot be excluded: X comes from it");
    }
    if exclusions.contains(Section::Info, "CSQ") {
        bail!("INFO/CSQ cannot be excluded: name its subfields as CSQ/subfield");
    }
    Ok(())
}

fn check_header(
    header: &Header,
    csq: Option<&CsqIndex>,
    exclusions: &Exclusions,
) -> anyhow::Result<()> {
    if csq.is_none() && header.infos().contains_key("CSQ") {
        bail!("INFO/CSQ has no \"Format: \" list in its Description");
    }
    if let Some(name) = exclusions.undeclared(header, csq).first() {
        bail!("--exclude-field {name} names nothing in this VCF's header");
    }
    let Some(ad) = header.formats().get("AD") else {
        bail!("the VCF declares no FORMAT/AD");
    };
    if ad.number() != format::Number::ReferenceAlternateBases || ad.ty() != format::Type::Integer {
        bail!("FORMAT/AD must be declared Number=R, Type=Integer");
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
}

fn partial_path(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(".partial");
    PathBuf::from(name)
}

fn write_output(path: &Path, output: Output) -> anyhow::Result<()> {
    let partial = partial_path(path);
    write_h5ad(&partial, output)
        .and_then(|()| Ok(fs::rename(&partial, path)?))
        .inspect_err(|_| {
            let _ = fs::remove_file(&partial);
        })
        .with_context(|| format!("cannot write {}", path.display()))
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
        list(
            "x_encoding",
            X_ENCODING.iter().map(ToString::to_string).collect(),
        ),
        list(
            "layer_encoding",
            LAYER_ENCODING.iter().map(ToString::to_string).collect(),
        ),
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
    match value? {
        SampleValue::Float(x) => Some(*x),
        SampleValue::Array(SampleArray::Float(values)) => values.get(position).copied().flatten(),
        _ => None,
    }
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

impl Lookups<'_> {
    fn pick(
        &self,
        entries: &[CsqEntry],
        allele: &str,
        alt: &str,
    ) -> anyhow::Result<(Option<Picked>, Option<&'static str>)> {
        let index = match self.csq {
            Some(index) if !entries.is_empty() => index,
            _ => return Ok((None, None)),
        };
        let candidates = entries_for_allele(entries, index, allele);
        if candidates.is_empty() && alt != "*" {
            bail!("no CSQ entry for ALT {alt}");
        }
        let ctx = PickContext {
            index,
            genes: self.genes,
            mane: self.mane,
        };
        let picked = pick(entries, &candidates, &ctx)?;
        Ok((picked, worst_panel_impact(entries, &candidates, &ctx)))
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
    records_without_alt: u64,
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
            var: VarTable::new(&plan.info, lookups.csq, exclusions, lookups.genes.is_some()),
            matrix: SparseBuilder::new(n_obs, &layers),
            n_sites_covered: vec![0; n_obs],
            n_sites_alt: vec![0; n_obs],
            skipped_alts: 0,
            records_without_alt: 0,
        }
    }

    fn push_record(&mut self, record: &RecordBuf) -> anyhow::Result<()> {
        let alts = record.alternate_bases().as_ref();
        let vep = vep_alleles(record.reference_bases(), alts);
        let entries = match self.lookups.csq {
            Some(index) => parse_entries(&raw_csq(record), index)?,
            None => Vec::new(),
        };
        let info = read_info(record, &self.plan.info, alts.len())?;
        let n_obs = self.n_sites_covered.len();
        let samples = RecordSamples::read(record, &self.plan.format, n_obs, alts.len())?;
        if alts.is_empty() {
            self.records_without_alt += 1;
        }

        for (position, (alt, allele)) in alts.iter().zip(&vep).enumerate() {
            if is_skipped_alt(alt) {
                self.skipped_alts += 1;
                continue;
            }
            let (picked, worst) = self.lookups.pick(&entries, allele, alt)?;
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
        let mut rows = Vec::new();
        let mut codes = Vec::new();
        let mut cells = vec![Vec::new(); slots.len()];
        let mut alleles = Vec::new();
        let mut totals = ColumnTotals::default();

        for (row, sample) in (0u32..).zip(0..n_obs) {
            samples.fill_alleles(sample, &mut alleles);
            let code = genotype_code(&alleles, allele);
            let depths = samples.allelic_depths(sample);
            if !covered(code.is_some(), depths) {
                continue;
            }
            rows.push(row);
            codes.push(code);
            for (column, &(layer, position, kind)) in cells.iter_mut().zip(&slots) {
                column.push(samples.cell(layer, sample, position, kind));
            }

            let alt_reads = reads_at(depths, allele);
            totals.n_cells_covered += 1;
            totals.ref_reads += reads_at(depths, 0);
            totals.alt_reads += alt_reads;
            self.n_sites_covered[sample] += 1;
            if alt_reads > 0 {
                totals.n_cells_alt += 1;
                self.n_sites_alt[sample] += 1;
            }
        }
        self.matrix.push_column(&rows, &codes, &cells);
        totals
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partial_files_keep_the_whole_output_name() {
        assert_eq!(
            partial_path(Path::new("dir/cells.1")),
            Path::new("dir/cells.1.partial")
        );
        assert_ne!(
            partial_path(Path::new("cells.1")),
            partial_path(Path::new("cells.2"))
        );
    }
}
