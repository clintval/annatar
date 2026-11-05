use std::collections::HashSet;

use anyhow::bail;

use crate::alleles::variant_name;
use crate::consequence::{is_coding, is_lof, most_severe};
use crate::csq::{CsqEntry, CsqIndex, is_missing};
use crate::fields::{Exclusions, FieldSpec, Kind, Section};
use crate::output::Column;
use crate::pick::Tier;
use crate::text::{CsqColumn, type_column};

const RECORD_NAMES: [&str; 11] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "id",
    "qual",
    "filter",
    "n_cells_covered",
    "n_cells_alt",
    "ref_reads",
    "alt_reads",
];

const DERIVED_NAMES: [&str; 5] = [
    "consequence_rank",
    "is_coding",
    "is_lof",
    "pick_tier",
    "worst_panel_impact",
];

#[derive(Clone, Debug, PartialEq)]
pub enum InfoValue {
    Int(Option<i64>),
    Float(Option<f64>),
    Flag(bool),
    Text(Option<String>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct RecordRow {
    pub chrom: String,
    pub pos: usize,
    pub reference: String,
    pub alt: String,
    pub id: Option<String>,
    pub qual: Option<f32>,
    pub filter: Option<String>,
    pub n_cells_covered: u32,
    pub n_cells_alt: u32,
    pub ref_reads: u64,
    pub alt_reads: u64,
}

enum InfoColumn {
    Int(Vec<Option<i64>>),
    Float(Vec<f64>),
    Flag(Vec<bool>),
    Text(Vec<Option<String>>),
}

impl InfoColumn {
    fn new(kind: Kind) -> Self {
        match kind {
            Kind::Integer => Self::Int(Vec::new()),
            Kind::Float => Self::Float(Vec::new()),
            Kind::Flag => Self::Flag(Vec::new()),
            Kind::Text => Self::Text(Vec::new()),
        }
    }

    fn push(&mut self, value: InfoValue) {
        match (self, value) {
            (Self::Int(v), InfoValue::Int(x)) => v.push(x),
            (Self::Float(v), InfoValue::Float(x)) => v.push(x.unwrap_or(f64::NAN)),
            (Self::Flag(v), InfoValue::Flag(x)) => v.push(x),
            (Self::Text(v), InfoValue::Text(x)) => v.push(x),
            _ => panic!("an INFO value does not match the kind of its column"),
        }
    }

    fn into_column(self) -> Column {
        match self {
            Self::Int(v) => Column::NullableInt(v),
            Self::Float(v) => Column::Float(v),
            Self::Flag(v) => Column::Flag(v),
            Self::Text(v) => Column::Text(v),
        }
    }
}

pub fn widen(value: f32) -> f64 {
    format!("{value}")
        .parse()
        .unwrap_or_else(|_| f64::from(value))
}

pub struct VarColumns {
    pub names: Vec<String>,
    pub columns: Vec<(String, Column)>,
}

pub struct VarTable {
    csq: Option<CsqIndex>,
    csq_columns: Vec<(String, Vec<Option<String>>)>,
    csq_positions: Vec<usize>,
    with_genes: bool,
    info: Vec<(String, InfoColumn)>,
    names: Vec<String>,
    seen: HashSet<String>,
    chrom: Vec<Option<String>>,
    pos: Vec<i64>,
    reference: Vec<Option<String>>,
    alt: Vec<Option<String>>,
    id: Vec<Option<String>>,
    qual: Vec<f64>,
    filter: Vec<Option<String>>,
    n_cells_covered: Vec<i64>,
    n_cells_alt: Vec<i64>,
    ref_reads: Vec<i64>,
    alt_reads: Vec<i64>,
    consequence_rank: Vec<Option<i64>>,
    is_coding: Vec<Option<bool>>,
    is_lof: Vec<Option<bool>>,
    pick_tier: Vec<Option<String>>,
    worst_panel_impact: Vec<Option<String>>,
}

impl VarTable {
    pub fn new(
        info: &[FieldSpec],
        csq: Option<&CsqIndex>,
        exclusions: &Exclusions,
        with_genes: bool,
    ) -> Self {
        let (csq_positions, csq_names): (Vec<usize>, Vec<String>) = csq
            .map(|index| {
                index
                    .names
                    .iter()
                    .enumerate()
                    .filter(|(_, name)| !exclusions.contains(Section::Csq, name))
                    .map(|(position, name)| (position, name.clone()))
                    .unzip()
            })
            .unwrap_or_default();

        let info = info
            .iter()
            .flat_map(|spec| {
                spec.column_names()
                    .into_iter()
                    .map(move |name| (name, spec.kind))
            })
            .map(|(name, kind)| {
                let taken = RECORD_NAMES.contains(&name.as_str())
                    || (csq.is_some() && DERIVED_NAMES.contains(&name.as_str()))
                    || csq_names.contains(&name);
                let name = if taken { format!("INFO_{name}") } else { name };
                (name, InfoColumn::new(kind))
            })
            .collect();

        Self {
            csq: csq.cloned(),
            csq_columns: csq_names.into_iter().map(|n| (n, Vec::new())).collect(),
            csq_positions,
            with_genes,
            info,
            names: Vec::new(),
            seen: HashSet::new(),
            chrom: Vec::new(),
            pos: Vec::new(),
            reference: Vec::new(),
            alt: Vec::new(),
            id: Vec::new(),
            qual: Vec::new(),
            filter: Vec::new(),
            n_cells_covered: Vec::new(),
            n_cells_alt: Vec::new(),
            ref_reads: Vec::new(),
            alt_reads: Vec::new(),
            consequence_rank: Vec::new(),
            is_coding: Vec::new(),
            is_lof: Vec::new(),
            pick_tier: Vec::new(),
            worst_panel_impact: Vec::new(),
        }
    }

    pub fn push(
        &mut self,
        row: RecordRow,
        info: Vec<InfoValue>,
        picked: Option<(&CsqEntry, Tier)>,
        worst_panel_impact: Option<&str>,
    ) -> anyhow::Result<()> {
        assert_eq!(
            info.len(),
            self.info.len(),
            "one INFO value is required per INFO column"
        );

        let (rank, coding, lof) = match picked.and_then(|(entry, _)| self.consequence(entry)) {
            Some(consequence) => {
                let (term, rank) = most_severe(consequence);
                (
                    Some(i64::from(rank)),
                    Some(is_coding(term)),
                    Some(is_lof(term)),
                )
            }
            None => (None, None, None),
        };
        let name = variant_name(&row.chrom, row.pos, &row.reference, &row.alt);
        if self.seen.contains(&name) {
            bail!(
                "variant {name} appears more than once; remove repeats with bcftools norm -d exact"
            );
        }
        let pos = i64::try_from(row.pos)?;
        let ref_reads = i64::try_from(row.ref_reads)?;
        let alt_reads = i64::try_from(row.alt_reads)?;

        for (column, value) in self.info.iter_mut().zip(info) {
            column.1.push(value);
        }
        for ((_, values), position) in self.csq_columns.iter_mut().zip(&self.csq_positions) {
            let value = picked
                .and_then(|(entry, _)| entry.values.get(*position))
                .filter(|value| !is_missing(value))
                .cloned();
            values.push(value);
        }

        self.seen.insert(name.clone());
        self.names.push(name);
        self.chrom.push(Some(row.chrom));
        self.pos.push(pos);
        self.reference.push(Some(row.reference));
        self.alt.push(Some(row.alt));
        self.id.push(row.id);
        self.qual.push(row.qual.map_or(f64::NAN, widen));
        self.filter.push(row.filter);
        self.n_cells_covered.push(i64::from(row.n_cells_covered));
        self.n_cells_alt.push(i64::from(row.n_cells_alt));
        self.ref_reads.push(ref_reads);
        self.alt_reads.push(alt_reads);
        self.consequence_rank.push(rank);
        self.is_coding.push(coding);
        self.is_lof.push(lof);
        self.pick_tier
            .push(picked.map(|(_, tier)| tier.as_str().to_string()));
        self.worst_panel_impact
            .push(worst_panel_impact.map(str::to_string));
        Ok(())
    }

    pub fn duplicate_column(&self) -> Option<&str> {
        let mut seen = HashSet::new();
        RECORD_NAMES
            .into_iter()
            .chain(self.csq_columns.iter().map(|(name, _)| name.as_str()))
            .chain(DERIVED_NAMES.into_iter().filter(|name| {
                self.csq.is_some() && (self.with_genes || *name != "worst_panel_impact")
            }))
            .chain(self.info.iter().map(|(name, _)| name.as_str()))
            .find(|name| !seen.insert(*name))
    }

    fn consequence<'a>(&self, entry: &'a CsqEntry) -> Option<&'a str> {
        entry.get(self.csq.as_ref()?, "Consequence")
    }

    pub fn finish(self) -> VarColumns {
        let record = [
            ("chrom", Column::Text(self.chrom)),
            ("pos", Column::Int(self.pos)),
            ("ref", Column::Text(self.reference)),
            ("alt", Column::Text(self.alt)),
            ("id", Column::Text(self.id)),
            ("qual", Column::Float(self.qual)),
            ("filter", Column::Text(self.filter)),
            ("n_cells_covered", Column::Int(self.n_cells_covered)),
            ("n_cells_alt", Column::Int(self.n_cells_alt)),
            ("ref_reads", Column::Int(self.ref_reads)),
            ("alt_reads", Column::Int(self.alt_reads)),
        ];
        let mut columns: Vec<(String, Column)> = record
            .into_iter()
            .map(|(name, column)| (name.to_string(), column))
            .collect();

        for (name, values) in self.csq_columns {
            let column = match type_column(&name, values) {
                CsqColumn::Integer(values) => Column::NullableInt(values),
                CsqColumn::Float(values) => Column::Float(values),
                CsqColumn::Text(values) => Column::Text(values),
            };
            columns.push((name, column));
        }

        if self.csq.is_some() {
            let derived = [
                (
                    "consequence_rank",
                    Column::NullableInt(self.consequence_rank),
                ),
                ("is_coding", Column::NullableFlag(self.is_coding)),
                ("is_lof", Column::NullableFlag(self.is_lof)),
                ("pick_tier", Column::Text(self.pick_tier)),
            ];
            columns.extend(derived.map(|(name, column)| (name.to_string(), column)));
            if self.with_genes {
                columns.push((
                    "worst_panel_impact".to_string(),
                    Column::Text(self.worst_panel_impact),
                ));
            }
        }

        for (name, column) in self.info {
            columns.push((name, column.into_column()));
        }

        VarColumns {
            names: self.names,
            columns,
        }
    }
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::*;
    use crate::consequence::rank;
    use crate::fields::Shape;
    use crate::h5ad::{data_frame, write_h5ad};
    use crate::output::Output;
    use nalgebra_sparse::CsrMatrix;
    use polars::prelude::{AnyValue, Column as FrameColumn, DataFrame, DataType};

    fn spec(id: &str, shape: Shape, kind: Kind) -> FieldSpec {
        FieldSpec {
            id: id.to_string(),
            shape,
            kind,
        }
    }

    fn index(names: &[&str]) -> CsqIndex {
        CsqIndex {
            names: names.iter().map(ToString::to_string).collect(),
        }
    }

    fn entry(values: &[&str]) -> CsqEntry {
        CsqEntry {
            values: values.iter().map(ToString::to_string).collect(),
        }
    }

    #[test]
    fn info_columns_expanding_to_one_name_are_a_duplicate() {
        let info = [
            spec("X", Shape::RefAlt, Kind::Integer),
            spec("X_ref", Shape::Scalar, Kind::Integer),
        ];
        let table = VarTable::new(&info, None, &Exclusions::default(), false);
        assert_eq!(table.duplicate_column(), Some("X_ref"));
    }

    #[test]
    fn a_csq_subfield_named_like_a_record_column_is_a_duplicate() {
        let csq = index(&["Allele", "pos"]);
        let table = VarTable::new(&[], Some(&csq), &Exclusions::default(), false);
        assert_eq!(table.duplicate_column(), Some("pos"));
    }

    #[test]
    fn worst_panel_impact_only_collides_with_genes() {
        let csq = index(&["worst_panel_impact"]);
        let without = VarTable::new(&[], Some(&csq), &Exclusions::default(), false);
        let with = VarTable::new(&[], Some(&csq), &Exclusions::default(), true);
        assert_eq!(without.duplicate_column(), None);
        assert_eq!(with.duplicate_column(), Some("worst_panel_impact"));
    }

    fn row(pos: usize, alt: &str) -> RecordRow {
        RecordRow {
            chrom: "chr1".into(),
            pos,
            reference: "A".into(),
            alt: alt.into(),
            id: None,
            qual: None,
            filter: None,
            n_cells_covered: 3,
            n_cells_alt: 1,
            ref_reads: 10,
            alt_reads: 4,
        }
    }

    fn is_categorical(column: &FrameColumn) -> bool {
        matches!(column.dtype(), DataType::Categorical(..))
    }

    fn names(df: &DataFrame) -> Vec<String> {
        df.get_column_names()
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    fn frame(columns: VarColumns) -> DataFrame {
        data_frame(columns.columns).unwrap()
    }

    fn nullable(columns: &VarColumns) -> Vec<&str> {
        columns
            .columns
            .iter()
            .filter(|(_, column)| {
                matches!(column, Column::NullableInt(_) | Column::NullableFlag(_))
            })
            .map(|(name, _)| name.as_str())
            .collect()
    }

    fn write(columns: VarColumns) -> anyhow::Result<()> {
        let rows = columns.names.len();
        let output = Output {
            obs_names: vec!["c0".into()],
            obs: Vec::new(),
            var_names: columns.names,
            var: columns.columns,
            x: CsrMatrix::try_from_csr_data(1, rows, vec![0, 0], vec![], Vec::<i8>::new()).unwrap(),
            layers: vec![],
            uns: vec![],
        };
        let dir = tempfile::tempdir()?;
        write_h5ad(&dir.path().join("out.h5ad"), output)
    }

    #[test]
    fn columns_have_spec_types() {
        let info = [
            spec("INDEL", Shape::Scalar, Kind::Flag),
            spec("AC", Shape::PerAlt, Kind::Integer),
            spec("AF", Shape::RefAlt, Kind::Float),
        ];
        let csq = index(&["Allele", "Consequence", "SYMBOL", "REVEL_score"]);
        let mut table = VarTable::new(&info, Some(&csq), &Exclusions::default(), false);
        let first = entry(&["T", "stop_gained", "TP53", "0.5"]);
        let second = entry(&["G", "missense_variant", "", ""]);
        let mut with_qual = row(100, "T");
        with_qual.qual = Some(30.5);
        let values = |flag| {
            vec![
                InfoValue::Flag(flag),
                InfoValue::Int(Some(2)),
                InfoValue::Float(Some(0.25)),
                InfoValue::Float(None),
            ]
        };
        let tier = Some((&first, Tier::Canonical));
        table.push(with_qual, values(true), tier, None).unwrap();
        let tier = Some((&second, Tier::Other));
        table
            .push(row(200, "G"), values(false), tier, None)
            .unwrap();
        let df = frame(table.finish());

        assert_eq!(df.height(), 2);
        for name in ["chrom", "ref", "alt", "id", "filter", "SYMBOL", "pick_tier"] {
            assert!(is_categorical(df.column(name).unwrap()), "{name}");
        }
        assert_eq!(df.column("pos").unwrap().dtype(), &DataType::Int64);
        assert_eq!(df.column("alt_reads").unwrap().dtype(), &DataType::Int64);
        assert_eq!(df.column("AC").unwrap().dtype(), &DataType::Int64);
        assert_eq!(df.column("INDEL").unwrap().dtype(), &DataType::Boolean);
        assert_eq!(df.column("is_lof").unwrap().dtype(), &DataType::Boolean);

        let qual = df.column("qual").unwrap().f64().unwrap();
        assert_eq!(qual.get(0), Some(30.5));
        assert!(qual.get(1).unwrap().is_nan());

        let revel = df.column("REVEL_score").unwrap().f64().unwrap();
        assert_eq!(revel.get(0), Some(0.5));
        assert!(revel.get(1).unwrap().is_nan());

        let symbol = df.column("SYMBOL").unwrap();
        assert_eq!(symbol.null_count(), 1);
        assert!(!matches!(symbol.get(0).unwrap(), AnyValue::Null));
    }

    #[test]
    fn info_collision_is_prefixed() {
        let info = [
            spec("SYMBOL", Shape::Scalar, Kind::Text),
            spec("DP", Shape::Scalar, Kind::Integer),
            spec("pick_tier", Shape::Scalar, Kind::Text),
            spec("worst_panel_impact", Shape::Scalar, Kind::Text),
            spec("alt", Shape::Fixed(2), Kind::Integer),
        ];
        let csq = index(&["SYMBOL"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), false);
        let df = frame(table.finish());
        let names = names(&df);
        for expected in [
            "INFO_SYMBOL",
            "INFO_pick_tier",
            "INFO_worst_panel_impact",
            "DP",
            "alt_1",
        ] {
            assert!(names.iter().any(|n| n == expected), "{expected}");
        }
        assert_eq!(names.iter().filter(|n| *n == "SYMBOL").count(), 1);
    }

    #[test]
    fn expanded_info_names_are_checked_for_collisions() {
        let info = [spec("AF", Shape::RefAlt, Kind::Float)];
        let csq = index(&["AF_ref"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), false);
        let df = frame(table.finish());
        let names = names(&df);
        assert!(names.contains(&"INFO_AF_ref".to_string()));
        assert!(names.contains(&"AF_alt".to_string()));
        assert_eq!(names.iter().filter(|n| *n == "AF_ref").count(), 1);
    }

    #[test]
    fn picked_values_land_in_their_rows() {
        let info = [
            spec("DP", Shape::Scalar, Kind::Integer),
            spec("INDEL", Shape::Scalar, Kind::Flag),
        ];
        let csq = index(&["Consequence", "SYMBOL", "REVEL_score"]);
        let mut table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        let missense = entry(&["missense_variant", "TP53", "0.9"]);
        let splice = entry(&["splice_donor_variant&intron_variant", "NPM1", "0.1"]);
        let intron = entry(&["intron_variant", "KIT", "0.5"]);
        let values = |dp, indel| vec![InfoValue::Int(Some(dp)), InfoValue::Flag(indel)];
        let picks = [
            (&missense, Tier::ManeSelect, "HIGH", 11, true),
            (&splice, Tier::Other, "LOW", 22, false),
            (&intron, Tier::Canonical, "MODIFIER", 33, true),
        ];
        for (i, (pick, tier, impact, dp, indel)) in picks.into_iter().enumerate() {
            table
                .push(
                    row(i + 1, "T"),
                    values(dp, indel),
                    Some((pick, tier)),
                    Some(impact),
                )
                .unwrap();
        }
        let df = frame(table.finish());

        let ints = |name: &str| -> Vec<Option<i64>> {
            df.column(name).unwrap().i64().unwrap().iter().collect()
        };
        let bools = |name: &str| -> Vec<Option<bool>> {
            df.column(name).unwrap().bool().unwrap().iter().collect()
        };
        let texts = |name: &str| -> Vec<Option<String>> {
            let column = df.column(name).unwrap().cast(&DataType::String).unwrap();
            column
                .str()
                .unwrap()
                .iter()
                .map(|v| v.map(str::to_string))
                .collect()
        };
        let ranks = ["missense_variant", "splice_donor_variant", "intron_variant"]
            .map(|term| Some(i64::from(rank(term).unwrap())));
        assert_eq!(ints("consequence_rank"), ranks);
        assert_eq!(bools("is_coding"), [Some(true), Some(false), Some(false)]);
        assert_eq!(bools("is_lof"), [Some(false), Some(true), Some(false)]);
        assert_eq!(
            texts("pick_tier"),
            [
                Some("mane_select".into()),
                Some("other".into()),
                Some("canonical".into())
            ]
        );
        assert_eq!(
            texts("worst_panel_impact"),
            [
                Some("HIGH".into()),
                Some("LOW".into()),
                Some("MODIFIER".into())
            ]
        );
        assert_eq!(
            texts("SYMBOL"),
            [Some("TP53".into()), Some("NPM1".into()), Some("KIT".into())]
        );
        let revel: Vec<Option<f64>> = df
            .column("REVEL_score")
            .unwrap()
            .f64()
            .unwrap()
            .iter()
            .collect();
        assert_eq!(revel, [Some(0.9), Some(0.1), Some(0.5)]);
        assert_eq!(ints("DP"), [Some(11), Some(22), Some(33)]);
        assert_eq!(bools("INDEL"), [Some(true), Some(false), Some(true)]);
    }

    #[test]
    fn csq_exclusion_drops_column() {
        let exclusions = Exclusions::new(vec!["INFO/CSQ/COSMIC".parse().unwrap()]);
        let info = [spec("COSMIC", Shape::Scalar, Kind::Text)];
        let csq = index(&["SYMBOL", "COSMIC"]);
        let table = VarTable::new(&info, Some(&csq), &exclusions, false);
        let df = frame(table.finish());
        let names = names(&df);
        assert!(names.contains(&"SYMBOL".to_string()));
        assert!(names.contains(&"COSMIC".to_string()));
        assert_eq!(names.iter().filter(|n| *n == "COSMIC").count(), 1);
        assert!(!names.contains(&"INFO_COSMIC".to_string()));
    }

    #[test]
    fn without_csq_derived_names_are_free_and_absent() {
        let info = [spec("pick_tier", Shape::Scalar, Kind::Text)];
        let table = VarTable::new(&info, None, &Exclusions::default(), false);
        assert_eq!(table.duplicate_column(), None);
        let columns = table.finish();
        assert_eq!(nullable(&columns).len(), 0);
        let names = names(&frame(columns));
        assert_eq!(names.len(), RECORD_NAMES.len() + 1);
        assert_eq!(names.last().unwrap(), "pick_tier");
    }

    #[test]
    fn only_columns_that_can_be_missing_are_nullable() {
        let info = [
            spec("DP", Shape::Scalar, Kind::Integer),
            spec("AF", Shape::Scalar, Kind::Float),
            spec("SOLO", Shape::Scalar, Kind::Flag),
            spec("NOTE", Shape::Scalar, Kind::Text),
        ];
        let csq = index(&["STRAND", "REVEL", "SYMBOL"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), false);
        assert_eq!(
            nullable(&table.finish()),
            ["STRAND", "consequence_rank", "is_coding", "is_lof", "DP"]
        );
    }

    #[test]
    fn a_csq_integer_that_fell_back_to_text_is_not_nullable() {
        let csq = index(&["STRAND"]);
        let mut table = VarTable::new(&[], Some(&csq), &Exclusions::default(), false);
        let plus = entry(&["plus"]);
        let tier = Some((&plus, Tier::Other));
        table.push(row(1, "T"), vec![], tier, None).unwrap();
        let columns = table.finish();
        assert!(!nullable(&columns).contains(&"STRAND"));
        assert!(is_categorical(frame(columns).column("STRAND").unwrap()));
    }

    #[test]
    fn column_order_is_fixed() {
        let info = [spec("DP", Shape::Scalar, Kind::Integer)];
        let csq = index(&["SYMBOL"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        let df = frame(table.finish());
        assert_eq!(
            names(&df),
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "id",
                "qual",
                "filter",
                "n_cells_covered",
                "n_cells_alt",
                "ref_reads",
                "alt_reads",
                "SYMBOL",
                "consequence_rank",
                "is_coding",
                "is_lof",
                "pick_tier",
                "worst_panel_impact",
                "DP",
            ]
        );
    }

    #[test]
    fn row_without_a_pick_is_empty() {
        let csq = index(&["Consequence", "SYMBOL"]);
        let mut table = VarTable::new(&[], Some(&csq), &Exclusions::default(), false);
        table.push(row(5, "T"), vec![], None, None).unwrap();
        let df = frame(table.finish());
        assert_eq!(df.height(), 1);
        for name in ["consequence_rank", "is_coding", "is_lof", "pick_tier"] {
            assert_eq!(df.column(name).unwrap().null_count(), 1, "{name}");
        }
        let symbol = df.column("SYMBOL").unwrap();
        assert!(is_categorical(symbol));
        assert_eq!(symbol.null_count(), 1);
    }

    #[test]
    fn a_duplicate_variant_is_refused_and_leaves_columns_aligned() {
        let info = [spec("DP", Shape::Scalar, Kind::Integer)];
        let csq = index(&["Consequence", "SYMBOL"]);
        let mut table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        let pick = entry(&["missense_variant", "TP53"]);
        let tier = Some((&pick, Tier::Canonical));
        table
            .push(
                row(100, "T"),
                vec![InfoValue::Int(Some(1))],
                tier,
                Some("HIGH"),
            )
            .unwrap();
        let error = table
            .push(
                row(100, "T"),
                vec![InfoValue::Int(Some(2))],
                tier,
                Some("HIGH"),
            )
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            "variant chr1:100:A:T appears more than once; \
            remove repeats with bcftools norm -d exact"
        );
        table
            .push(
                row(101, "T"),
                vec![InfoValue::Int(Some(3))],
                tier,
                Some("HIGH"),
            )
            .unwrap();
        let columns = table.finish();
        assert_eq!(columns.names.len(), 2);
        let frame = frame(columns);
        assert_eq!(frame.height(), 2);
        assert!(frame.columns().iter().all(|c| c.len() == 2));
    }

    #[test]
    fn f32_values_widen_through_their_shortest_decimal_form() {
        assert_eq!(widen(30.1), 30.1);
        assert_eq!(widen(0.1), 0.1);
        assert!(widen(f32::NAN).is_nan());
        assert_eq!(widen(f32::INFINITY), f64::INFINITY);
        assert_eq!(widen(f32::NEG_INFINITY), f64::NEG_INFINITY);
    }

    #[test]
    fn info_shapes_expand_in_order() {
        let info = [
            spec("DP4", Shape::Fixed(2), Kind::Integer),
            spec("AF", Shape::RefAlt, Kind::Float),
        ];
        let mut table = VarTable::new(&info, None, &Exclusions::default(), false);
        let values = vec![
            InfoValue::Int(Some(7)),
            InfoValue::Int(Some(8)),
            InfoValue::Float(Some(0.1)),
            InfoValue::Float(Some(0.9)),
        ];
        table.push(row(1, "T"), values, None, None).unwrap();
        let df = frame(table.finish());
        let int = |name: &str| df.column(name).unwrap().i64().unwrap().get(0);
        let float = |name: &str| df.column(name).unwrap().f64().unwrap().get(0);
        assert_eq!(int("DP4_1"), Some(7));
        assert_eq!(int("DP4_2"), Some(8));
        assert_eq!(float("AF_ref"), Some(0.1));
        assert_eq!(float("AF_alt"), Some(0.9));
    }

    #[test]
    fn names_come_from_variant_name() {
        let mut table = VarTable::new(&[], None, &Exclusions::default(), false);
        table.push(row(100, "T"), vec![], None, None).unwrap();
        table.push(row(50, "G"), vec![], None, None).unwrap();
        let names = table.finish().names;
        assert_eq!(
            names,
            [
                variant_name("chr1", 100, "A", "T"),
                variant_name("chr1", 50, "A", "G")
            ]
        );
    }

    #[test]
    fn zero_rows_finish() {
        let info = [spec("DP", Shape::Scalar, Kind::Integer)];
        let csq = index(&["SYMBOL"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        let columns = table.finish();
        assert_eq!(columns.names.len(), 0);
        let df = data_frame(columns.columns.clone()).unwrap();
        assert_eq!(df.height(), 0);
        assert_eq!(df.width(), 11 + 1 + 5 + 1);
        assert!(is_categorical(df.column("SYMBOL").unwrap()));
        write(columns).unwrap();
    }

    #[test]
    fn all_missing_text_columns_write() {
        let info = [spec("NOTE", Shape::Scalar, Kind::Text)];
        let csq = index(&["Consequence", "SYMBOL"]);
        let mut table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        table
            .push(row(1, "T"), vec![InfoValue::Text(None)], None, None)
            .unwrap();
        table
            .push(row(2, "G"), vec![InfoValue::Text(None)], None, None)
            .unwrap();
        let columns = table.finish();
        let df = data_frame(columns.columns.clone()).unwrap();
        for name in [
            "id",
            "filter",
            "pick_tier",
            "worst_panel_impact",
            "NOTE",
            "SYMBOL",
        ] {
            assert_eq!(df.column(name).unwrap().null_count(), 2, "{name}");
        }
        write(columns).unwrap();
    }

    #[test]
    fn no_csq_rows_write() {
        let mut table = VarTable::new(&[], None, &Exclusions::default(), false);
        table.push(row(1, "T"), vec![], None, None).unwrap();
        write(table.finish()).unwrap();
    }
}
