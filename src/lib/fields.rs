use std::collections::HashSet;
use std::str::FromStr;

use noodles::vcf::Header;
use noodles::vcf::header::record::value::map::{format, info};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Section {
    Info,
    Format,
    Csq,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Exclusion {
    pub section: Section,
    pub id: String,
}

impl FromStr for Exclusion {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let invalid = || "expected INFO/ID, FORMAT/ID or CSQ/subfield".to_string();
        let (prefix, id) = s.split_once('/').ok_or_else(invalid)?;
        let section = match prefix {
            "INFO" => Section::Info,
            "FORMAT" => Section::Format,
            "CSQ" => Section::Csq,
            _ => return Err(invalid()),
        };
        if id.is_empty() {
            return Err(invalid());
        }
        Ok(Self {
            section,
            id: id.to_string(),
        })
    }
}

#[derive(Clone, Debug, Default)]
pub struct Exclusions(HashSet<(Section, String)>);

impl Exclusions {
    pub fn new(items: Vec<Exclusion>) -> Self {
        Self(items.into_iter().map(|e| (e.section, e.id)).collect())
    }

    pub fn contains(&self, section: Section, id: &str) -> bool {
        self.0.contains(&(section, id.to_string()))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    Scalar,
    PerAlt,
    RefAlt,
    Fixed(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Integer,
    Float,
    Flag,
    Text,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FieldSpec {
    pub id: String,
    pub shape: Shape,
    pub kind: Kind,
}

impl FieldSpec {
    pub fn column_names(&self) -> Vec<String> {
        match self.shape {
            Shape::Scalar | Shape::PerAlt => vec![self.id.clone()],
            Shape::RefAlt => vec![format!("{}_ref", self.id), format!("{}_alt", self.id)],
            Shape::Fixed(n) => (1..=n).map(|i| format!("{}_{i}", self.id)).collect(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FieldPlan {
    pub info: Vec<FieldSpec>,
    pub format: Vec<FieldSpec>,
    pub skipped: Vec<String>,
}

fn info_shape(number: info::Number) -> Option<Shape> {
    match number {
        info::Number::Count(0 | 1) => Some(Shape::Scalar),
        info::Number::Count(n) => Some(Shape::Fixed(n)),
        info::Number::AlternateBases => Some(Shape::PerAlt),
        info::Number::ReferenceAlternateBases => Some(Shape::RefAlt),
        _ => None,
    }
}

fn format_shape(number: format::Number) -> Option<Shape> {
    match number {
        format::Number::Count(0 | 1) => Some(Shape::Scalar),
        format::Number::Count(n) => Some(Shape::Fixed(n)),
        format::Number::AlternateBases => Some(Shape::PerAlt),
        format::Number::ReferenceAlternateBases => Some(Shape::RefAlt),
        _ => None,
    }
}

fn info_kind(ty: info::Type) -> Kind {
    match ty {
        info::Type::Integer => Kind::Integer,
        info::Type::Float => Kind::Float,
        info::Type::Flag => Kind::Flag,
        info::Type::Character | info::Type::String => Kind::Text,
    }
}

fn format_kind(ty: format::Type) -> Kind {
    match ty {
        format::Type::Integer => Kind::Integer,
        format::Type::Float => Kind::Float,
        format::Type::Character | format::Type::String => Kind::Text,
    }
}

pub fn plan_fields(header: &Header, exclusions: &Exclusions) -> FieldPlan {
    let mut plan = FieldPlan::default();

    for (id, map) in header.infos() {
        if id == "CSQ" || exclusions.contains(Section::Info, id) {
            continue;
        }
        match info_shape(map.number()) {
            Some(shape) => plan.info.push(FieldSpec {
                id: id.clone(),
                shape,
                kind: info_kind(map.ty()),
            }),
            None => plan.skipped.push(format!("INFO/{id}")),
        }
    }

    for (id, map) in header.formats() {
        if id == "GT" || exclusions.contains(Section::Format, id) {
            continue;
        }
        let kind = format_kind(map.ty());
        match format_shape(map.number()) {
            Some(shape) if kind != Kind::Text => plan.format.push(FieldSpec {
                id: id.clone(),
                shape,
                kind,
            }),
            _ => plan.skipped.push(format!("FORMAT/{id}")),
        }
    }

    plan
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;
    use std::path::PathBuf;

    use crate::input::open_variants;

    use super::*;

    fn bcftools_plan(exclusions: &[&str]) -> FieldPlan {
        plan_for("bcftools-1.24.vcf", exclusions)
    }

    fn plan_for(name: &str, exclusions: &[&str]) -> FieldPlan {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name);
        let source = open_variants(&path, NonZeroUsize::new(1).unwrap()).unwrap();
        let items = exclusions.iter().map(|s| s.parse().unwrap()).collect();
        plan_fields(&source.header, &Exclusions::new(items))
    }

    fn columns(specs: &[FieldSpec]) -> Vec<String> {
        specs.iter().flat_map(FieldSpec::column_names).collect()
    }

    #[test]
    fn info_columns_follow_header_order_and_expand() {
        let plan = bcftools_plan(&[]);
        assert_eq!(
            columns(&plan.info),
            [
                "INDEL", "IDV", "IMF", "DP", "VDB", "RPBZ", "MQBZ", "BQBZ", "MQSBZ", "SCBZ", "SGB",
                "MQ0F", "AC", "AN", "DP4_1", "DP4_2", "DP4_3", "DP4_4", "MQ"
            ]
        );
    }

    #[test]
    fn format_columns_split_ref_and_alt_and_skip_pl() {
        let plan = bcftools_plan(&[]);
        assert_eq!(columns(&plan.format), ["DP", "AD_ref", "AD_alt"]);
        assert_eq!(plan.skipped, ["FORMAT/PL"]);
    }

    #[test]
    fn genotype_is_neither_a_layer_nor_listed() {
        let plan = bcftools_plan(&[]);
        assert!(plan.format.iter().all(|spec| spec.id != "GT"));
        assert!(!plan.skipped.iter().any(|s| s.ends_with("/GT")));
    }

    #[test]
    fn excluded_format_field_is_absent_and_not_skipped() {
        let plan = bcftools_plan(&["FORMAT/DP", "FORMAT/PL"]);
        assert_eq!(columns(&plan.format), ["AD_ref", "AD_alt"]);
        assert_eq!(plan.skipped, Vec::<String>::new());
    }

    #[test]
    fn excluded_info_field_is_absent() {
        let plan = bcftools_plan(&["INFO/MQ"]);
        assert_eq!(plan.skipped, bcftools_plan(&[]).skipped);
        assert!(!columns(&plan.info).contains(&"MQ".to_string()));
    }

    #[test]
    fn field_shapes_cover_every_skip_and_kind() {
        let plan = plan_for("field-shapes.vcf", &[]);
        assert_eq!(columns(&plan.info), ["NOTE"]);
        assert_eq!(plan.info[0].kind, Kind::Text);
        assert_eq!(
            columns(&plan.format),
            ["AD_ref", "AD_alt", "PAIR_1", "PAIR_2"]
        );
        assert_eq!(plan.skipped, ["INFO/GL", "INFO/VAR", "FORMAT/FT"]);
    }

    #[test]
    fn exclusion_parses_section_and_id() {
        let parsed: Exclusion = "CSQ/COSMIC".parse().unwrap();
        assert_eq!(parsed.section, Section::Csq);
        assert_eq!(parsed.id, "COSMIC");
    }

    #[test]
    fn exclusion_without_section_is_rejected() {
        let error = "COSMIC".parse::<Exclusion>().unwrap_err();
        assert_eq!(error, "expected INFO/ID, FORMAT/ID or CSQ/subfield");
    }

    #[test]
    fn exclusions_match_on_section_and_id() {
        let exclusions = Exclusions::new(vec!["INFO/DP".parse().unwrap()]);
        assert!(exclusions.contains(Section::Info, "DP"));
        assert!(!exclusions.contains(Section::Format, "DP"));
    }
}
