pub fn vep_alleles(reference: &str, alternates: &[String]) -> Vec<String> {
    let first = reference.chars().next();
    let is_indel = alternates.iter().any(|alt| alt.len() != reference.len());
    let shared = is_indel
        && first.is_some_and(|base| {
            alternates
                .iter()
                .filter(|alt| !alt.contains('*'))
                .all(|alt| alt.starts_with(base))
        });
    if !shared {
        return alternates.to_vec();
    }
    alternates
        .iter()
        .map(|alt| {
            if alt.contains('*') {
                return alt.clone();
            }
            let trimmed = &alt[first.map_or(0, char::len_utf8)..];
            if trimmed.is_empty() {
                "-".to_string()
            } else {
                trimmed.to_string()
            }
        })
        .collect()
}

pub fn is_skipped_alt(alt: &str) -> bool {
    alt == "<*>" || alt == "<NON_REF>"
}

pub fn variant_name(chrom: &str, pos: usize, reference: &str, alt: &str) -> String {
    format!("{chrom}:{pos}:{reference}:{alt}")
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn owned(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| (*s).to_string()).collect()
    }

    #[rstest]
    #[case("G", &["A"], &["A"])]
    #[case("GAAA", &["G", "GA"], &["-", "A"])]
    #[case("A", &["G", "AT"], &["G", "AT"])]
    #[case("A", &["*"], &["*"])]
    #[case("A", &["<DEL>"], &["<DEL>"])]
    #[case("AC", &["AT"], &["AT"])]
    #[case("AC", &["AT", "A"], &["T", "-"])]
    #[case("AT", &["A", "*"], &["-", "*"])]
    #[case("A", &["AT", "<*>"], &["T", "<*>"])]
    #[case("AC", &["AT", "<*>"], &["T", "<*>"])]
    #[case("A", &["AT", "<NON_REF>"], &["AT", "<NON_REF>"])]
    #[case("AT", &["*"], &["*"])]
    #[case("A", &["G", "T"], &["G", "T"])]
    #[case("AT", &["A"], &["-"])]
    fn vep_alleles_trims_shared_first_base(
        #[case] reference: &str,
        #[case] alternates: &[&str],
        #[case] expected: &[&str],
    ) {
        assert_eq!(vep_alleles(reference, &owned(alternates)), owned(expected));
    }

    #[test]
    fn skipped_alts_are_the_gvcf_placeholders() {
        assert!(is_skipped_alt("<*>"));
        assert!(is_skipped_alt("<NON_REF>"));
        assert!(!is_skipped_alt("*"));
    }

    #[test]
    fn variant_name_joins_with_colons() {
        assert_eq!(
            variant_name("chr17", 7_673_704, "G", "A"),
            "chr17:7673704:G:A"
        );
    }
}
