use polars::prelude::{Categories, Column, DataFrame, DataType, NamedFrom, Series};

use crate::alleles::variant_name;
use crate::consequence::{is_coding, is_lof, most_severe};
use crate::csq::{CsqEntry, CsqIndex};
use crate::fields::{Exclusions, FieldSpec, Kind, Section};
use crate::pick::Tier;
use crate::text::{TextColumn, type_column};

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

    fn into_column(self, name: &str) -> anyhow::Result<Column> {
        Ok(match self {
            Self::Int(v) => Series::new(name.into(), v).into(),
            Self::Float(v) => Series::new(name.into(), v).into(),
            Self::Flag(v) => Series::new(name.into(), v).into(),
            Self::Text(v) => categorical(name, v)?,
        })
    }
}

fn categorical(name: &str, values: Vec<Option<String>>) -> anyhow::Result<Column> {
    Ok(Series::new(name.into(), values)
        .cast(&DataType::from_categories(Categories::global()))?
        .into())
}

pub struct VarTable {
    csq: Option<CsqIndex>,
    csq_columns: Vec<(String, Vec<Option<String>>)>,
    with_genes: bool,
    info: Vec<(String, InfoColumn)>,
    names: Vec<String>,
    chrom: Vec<Option<String>>,
    pos: Vec<Option<i64>>,
    reference: Vec<Option<String>>,
    alt: Vec<Option<String>>,
    id: Vec<Option<String>>,
    qual: Vec<f64>,
    filter: Vec<Option<String>>,
    n_cells_covered: Vec<Option<i64>>,
    n_cells_alt: Vec<Option<i64>>,
    ref_reads: Vec<Option<i64>>,
    alt_reads: Vec<Option<i64>>,
    consequence_rank: Vec<Option<i64>>,
    is_coding: Vec<bool>,
    is_lof: Vec<bool>,
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
        let csq_names: Vec<String> = csq
            .map(|index| {
                index
                    .names
                    .iter()
                    .filter(|name| !exclusions.contains(Section::Csq, name))
                    .cloned()
                    .collect()
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
                    || DERIVED_NAMES.contains(&name.as_str())
                    || csq_names.contains(&name);
                let name = if taken { format!("INFO_{name}") } else { name };
                (name, InfoColumn::new(kind))
            })
            .collect();

        Self {
            csq: csq.cloned(),
            csq_columns: csq_names.into_iter().map(|n| (n, Vec::new())).collect(),
            with_genes,
            info,
            names: Vec::new(),
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
                let (term, rank) = most_severe(consequence)?;
                (Some(i64::from(rank)), is_coding(term), is_lof(term))
            }
            None => (None, false, false),
        };
        let pos = i64::try_from(row.pos)?;
        let ref_reads = i64::try_from(row.ref_reads)?;
        let alt_reads = i64::try_from(row.alt_reads)?;

        for (column, value) in self.info.iter_mut().zip(info) {
            column.1.push(value);
        }
        for (name, values) in &mut self.csq_columns {
            let value = picked
                .zip(self.csq.as_ref())
                .and_then(|((entry, _), index)| entry.get(index, name).map(str::to_string));
            values.push(value);
        }

        self.names
            .push(variant_name(&row.chrom, row.pos, &row.reference, &row.alt));
        self.chrom.push(Some(row.chrom));
        self.pos.push(Some(pos));
        self.reference.push(Some(row.reference));
        self.alt.push(Some(row.alt));
        self.id.push(row.id);
        self.qual.push(row.qual.map_or(f64::NAN, f64::from));
        self.filter.push(row.filter);
        self.n_cells_covered
            .push(Some(i64::from(row.n_cells_covered)));
        self.n_cells_alt.push(Some(i64::from(row.n_cells_alt)));
        self.ref_reads.push(Some(ref_reads));
        self.alt_reads.push(Some(alt_reads));
        self.consequence_rank.push(rank);
        self.is_coding.push(coding);
        self.is_lof.push(lof);
        self.pick_tier
            .push(picked.map(|(_, tier)| tier.as_str().to_string()));
        self.worst_panel_impact
            .push(worst_panel_impact.map(str::to_string));
        Ok(())
    }

    fn consequence<'a>(&self, entry: &'a CsqEntry) -> Option<&'a str> {
        entry.get(self.csq.as_ref()?, "Consequence")
    }

    pub fn finish(self) -> anyhow::Result<(Vec<String>, DataFrame)> {
        let mut columns: Vec<Column> = vec![
            categorical("chrom", self.chrom)?,
            Series::new("pos".into(), self.pos).into(),
            categorical("ref", self.reference)?,
            categorical("alt", self.alt)?,
            categorical("id", self.id)?,
            Series::new("qual".into(), self.qual).into(),
            categorical("filter", self.filter)?,
            Series::new("n_cells_covered".into(), self.n_cells_covered).into(),
            Series::new("n_cells_alt".into(), self.n_cells_alt).into(),
            Series::new("ref_reads".into(), self.ref_reads).into(),
            Series::new("alt_reads".into(), self.alt_reads).into(),
        ];

        for (name, values) in self.csq_columns {
            columns.push(match type_column(values) {
                TextColumn::Float(values) => Series::new(name.as_str().into(), values).into(),
                TextColumn::Text(values) => categorical(&name, values)?,
            });
        }

        columns.push(Series::new("consequence_rank".into(), self.consequence_rank).into());
        columns.push(Series::new("is_coding".into(), self.is_coding).into());
        columns.push(Series::new("is_lof".into(), self.is_lof).into());
        columns.push(categorical("pick_tier", self.pick_tier)?);
        if self.with_genes {
            columns.push(categorical("worst_panel_impact", self.worst_panel_impact)?);
        }

        for (name, column) in self.info {
            columns.push(column.into_column(&name)?);
        }

        Ok((self.names, DataFrame::new_infer_height(columns)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consequence::rank;
    use crate::fields::Shape;
    use crate::h5ad::{Output, write_h5ad};
    use nalgebra_sparse::CsrMatrix;
    use polars::prelude::AnyValue;

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

    fn is_categorical(column: &Column) -> bool {
        matches!(column.dtype(), DataType::Categorical(..))
    }

    fn names(df: &DataFrame) -> Vec<String> {
        df.get_column_names()
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    fn write(names: Vec<String>, var: DataFrame) -> anyhow::Result<()> {
        let rows = names.len();
        let output = Output {
            obs_names: vec!["c0".into()],
            obs: DataFrame::empty(),
            var_names: names,
            var,
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
        let (_, df) = table.finish().unwrap();

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
        let (_, df) = table.finish().unwrap();
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
        let (_, df) = table.finish().unwrap();
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
        let lof = entry(&["stop_gained", "TP53", "0.9"]);
        let silent = entry(&["synonymous_variant&intron_variant", "NPM1", "0.1"]);
        let values = |dp, indel| vec![InfoValue::Int(Some(dp)), InfoValue::Flag(indel)];
        table
            .push(
                row(1, "T"),
                values(11, true),
                Some((&lof, Tier::ManeSelect)),
                Some("HIGH"),
            )
            .unwrap();
        table
            .push(
                row(2, "G"),
                values(22, false),
                Some((&silent, Tier::Other)),
                Some("LOW"),
            )
            .unwrap();
        let (_, df) = table.finish().unwrap();

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
        let stop_rank = i64::from(rank("stop_gained").unwrap());
        let silent_rank = i64::from(rank("synonymous_variant").unwrap());
        assert_eq!(
            ints("consequence_rank"),
            [Some(stop_rank), Some(silent_rank)]
        );
        assert_eq!(bools("is_coding"), [Some(true), Some(true)]);
        assert_eq!(bools("is_lof"), [Some(true), Some(false)]);
        assert_eq!(
            texts("pick_tier"),
            [Some("mane_select".into()), Some("other".into())]
        );
        assert_eq!(
            texts("worst_panel_impact"),
            [Some("HIGH".into()), Some("LOW".into())]
        );
        assert_eq!(texts("SYMBOL"), [Some("TP53".into()), Some("NPM1".into())]);
        let revel: Vec<Option<f64>> = df
            .column("REVEL_score")
            .unwrap()
            .f64()
            .unwrap()
            .iter()
            .collect();
        assert_eq!(revel, [Some(0.9), Some(0.1)]);
        assert_eq!(ints("DP"), [Some(11), Some(22)]);
        assert_eq!(bools("INDEL"), [Some(true), Some(false)]);
    }

    #[test]
    fn csq_exclusion_drops_column() {
        let exclusions = Exclusions::new(vec!["CSQ/COSMIC".parse().unwrap()]);
        let info = [spec("COSMIC", Shape::Scalar, Kind::Text)];
        let csq = index(&["SYMBOL", "COSMIC"]);
        let table = VarTable::new(&info, Some(&csq), &exclusions, false);
        let (_, df) = table.finish().unwrap();
        let names = names(&df);
        assert!(names.contains(&"SYMBOL".to_string()));
        assert!(names.contains(&"COSMIC".to_string()));
        assert_eq!(names.iter().filter(|n| *n == "COSMIC").count(), 1);
        assert!(!names.contains(&"INFO_COSMIC".to_string()));
    }

    #[test]
    fn column_order_is_fixed() {
        let info = [spec("DP", Shape::Scalar, Kind::Integer)];
        let csq = index(&["SYMBOL"]);
        let table = VarTable::new(&info, Some(&csq), &Exclusions::default(), true);
        let (_, df) = table.finish().unwrap();
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
        let (_, df) = table.finish().unwrap();
        assert_eq!(df.height(), 1);
        assert_eq!(
            df.column("is_lof").unwrap().bool().unwrap().get(0),
            Some(false)
        );
        assert_eq!(
            df.column("is_coding").unwrap().bool().unwrap().get(0),
            Some(false)
        );
        for name in ["consequence_rank", "pick_tier"] {
            assert_eq!(df.column(name).unwrap().null_count(), 1, "{name}");
        }
        let symbol = df.column("SYMBOL").unwrap().f64().unwrap();
        assert!(symbol.get(0).unwrap().is_nan());
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
        let (_, df) = table.finish().unwrap();
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
        let (names, _) = table.finish().unwrap();
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
        let (names, df) = table.finish().unwrap();
        assert_eq!(names.len(), 0);
        assert_eq!(df.height(), 0);
        assert_eq!(df.width(), 11 + 1 + 5 + 1);
        write(names, df).unwrap();
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
        let (names, df) = table.finish().unwrap();
        for name in ["id", "filter", "pick_tier", "worst_panel_impact", "NOTE"] {
            assert_eq!(df.column(name).unwrap().null_count(), 2, "{name}");
        }
        write(names, df).unwrap();
    }

    #[test]
    fn no_csq_rows_write() {
        let mut table = VarTable::new(&[], None, &Exclusions::default(), false);
        table.push(row(1, "T"), vec![], None, None).unwrap();
        let (names, df) = table.finish().unwrap();
        write(names, df).unwrap();
    }
}
