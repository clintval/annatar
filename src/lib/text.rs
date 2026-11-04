use log::warn;

const INTEGERS: [&str; 12] = [
    "ALLELE_NUM",
    "DISTANCE",
    "GENE_PHENO",
    "HGVS_OFFSET",
    "MINIMISED",
    "MOTIF_POS",
    "OverlapBP",
    "PICK",
    "REFSEQ_OFFSET",
    "SHIFT_LENGTH",
    "STRAND",
    "TSL",
];

const FLOATS: [&str; 9] = [
    "CADD_PHRED",
    "CADD_RAW",
    "CADD_phred",
    "CADD_raw",
    "MOTIF_SCORE_CHANGE",
    "OverlapPC",
    "REVEL",
    "REVEL_score",
    "am_pathogenicity",
];

#[derive(Clone, Debug, PartialEq)]
pub enum CsqColumn {
    Integer(Vec<Option<i64>>),
    Float(Vec<f64>),
    Text(Vec<Option<String>>),
}

fn parse_every<T>(
    name: &str,
    wanted: &str,
    values: &[Option<String>],
    parse: impl Fn(&str) -> Option<T>,
) -> Option<Vec<Option<T>>> {
    values
        .iter()
        .map(|value| {
            let Some(value) = value else {
                return Some(None);
            };
            let parsed = parse(value);
            if parsed.is_none() {
                warn!(
                    "INFO/CSQ/{name}: {value} is not {wanted}, so the column is text in this file"
                );
            }
            parsed.map(Some)
        })
        .collect()
}

fn largest_part(value: &str) -> Option<f64> {
    value
        .split('&')
        .filter(|part| !part.is_empty() && *part != ".")
        .map(|part| part.parse::<f64>().ok().filter(|n| n.is_finite()))
        .try_fold(f64::NAN, |largest, n| n.map(|n| largest.max(n)))
}

pub fn type_column(name: &str, values: Vec<Option<String>>) -> CsqColumn {
    let typed = if INTEGERS.contains(&name) {
        parse_every(name, "an integer", &values, |value| value.parse().ok()).map(CsqColumn::Integer)
    } else if name == "AF" || name.ends_with("_AF") || FLOATS.contains(&name) {
        parse_every(name, "a number", &values, largest_part).map(|numbers| {
            CsqColumn::Float(
                numbers
                    .into_iter()
                    .map(|number| number.unwrap_or(f64::NAN))
                    .collect(),
            )
        })
    } else {
        None
    };
    typed.unwrap_or(CsqColumn::Text(values))
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn some(values: &[Option<&str>]) -> Vec<Option<String>> {
        values.iter().map(|v| v.map(str::to_string)).collect()
    }

    fn floats(column: CsqColumn) -> Vec<f64> {
        match column {
            CsqColumn::Float(values) => values,
            other => panic!("expected a float column, got {other:?}"),
        }
    }

    #[rstest]
    #[case("AF")]
    #[case("gnomADe_AF")]
    #[case("MAX_AF")]
    #[case("CADD_PHRED")]
    #[case("CADD_RAW")]
    #[case("CADD_phred")]
    #[case("CADD_raw")]
    #[case("MOTIF_SCORE_CHANGE")]
    #[case("OverlapPC")]
    #[case("REVEL")]
    #[case("REVEL_score")]
    #[case("am_pathogenicity")]
    fn the_dictionary_names_every_float(#[case] name: &str) {
        assert!(matches!(
            type_column(name, some(&[Some("0.5")])),
            CsqColumn::Float(_)
        ));
    }

    #[rstest]
    #[case("ALLELE_NUM")]
    #[case("DISTANCE")]
    #[case("GENE_PHENO")]
    #[case("HGVS_OFFSET")]
    #[case("MINIMISED")]
    #[case("MOTIF_POS")]
    #[case("OverlapBP")]
    #[case("PICK")]
    #[case("REFSEQ_OFFSET")]
    #[case("SHIFT_LENGTH")]
    #[case("STRAND")]
    #[case("TSL")]
    fn the_dictionary_names_every_integer(#[case] name: &str) {
        assert_eq!(
            type_column(name, some(&[Some("1"), None])),
            CsqColumn::Integer(vec![Some(1), None])
        );
    }

    #[rstest]
    #[case("PUBMED")]
    #[case("ClinVar")]
    #[case("SOMATIC")]
    #[case("MAX_AF_POPS")]
    #[case("cDNA_position")]
    #[case("Unknown_plugin")]
    fn every_other_name_is_verbatim_text(#[case] name: &str) {
        let values = some(&[Some("1&1&.&"), None]);
        assert_eq!(type_column(name, values.clone()), CsqColumn::Text(values));
    }

    #[test]
    fn a_column_with_no_values_keeps_its_type() {
        assert_eq!(
            type_column("DISTANCE", some(&[None, None])),
            CsqColumn::Integer(vec![None, None])
        );
        assert!(
            floats(type_column("REVEL", some(&[None])))
                .iter()
                .all(|v| v.is_nan())
        );
        assert_eq!(
            type_column("SYMBOL", some(&[None])),
            CsqColumn::Text(vec![None])
        );
    }

    #[test]
    fn floats_keep_the_largest_part_ignoring_empty_and_dot_parts() {
        let values = floats(type_column(
            "gnomADe_AF",
            some(&[Some("0.1&0.7"), Some("0.1&.&&0.3"), Some(".&&"), None]),
        ));
        assert_eq!(values[..2], [0.7, 0.3]);
        assert!(values[2].is_nan() && values[3].is_nan());
    }

    #[rstest]
    #[case("REVEL", "nan")]
    #[case("REVEL", "inf")]
    #[case("REVEL", "-infinity")]
    #[case("REVEL", "1e999")]
    #[case("REVEL", "0.1&x")]
    #[case("TSL", "1&2")]
    #[case("TSL", "1.5")]
    #[case("TSL", "tsl1")]
    fn a_value_failing_its_type_keeps_the_column_verbatim(#[case] name: &str, #[case] value: &str) {
        let values = some(&[Some("1"), Some(value), None]);
        assert_eq!(type_column(name, values.clone()), CsqColumn::Text(values));
    }
}
