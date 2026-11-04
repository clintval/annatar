pub fn vep_alleles(reference: &str, alternates: &[String], svtype: Option<&str>) -> Vec<String> {
    if takes_the_structural_path(reference, alternates, svtype) {
        structural_alleles(alternates, svtype)
    } else {
        sequence_alleles(reference, alternates)
    }
}

fn takes_the_structural_path(reference: &str, alternates: &[String], svtype: Option<&str>) -> bool {
    let is_acgt = |allele: &str| allele.bytes().all(|base| b"ACGT".contains(&base));
    let all_acgt = is_acgt(reference) && alternates.iter().all(|alt| is_acgt(alt));
    !all_acgt
        && (svtype.is_some()
            || alternates
                .iter()
                .any(|alt| alt != "<*>" && is_structural(alt)))
}

fn structural_alleles(alternates: &[String], svtype: Option<&str>) -> Vec<String> {
    let joined = alternates.join("/");
    let kind = svtype
        .filter(|_| !names_its_type(&joined))
        .unwrap_or(&joined);
    let term = so_term(kind);
    alternates
        .iter()
        .map(|alt| match term {
            Some(term) if !is_breakend(alt) => term.to_string(),
            _ => alt.clone(),
        })
        .collect()
}

fn names_its_type(alt: &str) -> bool {
    let alt = alt.strip_prefix('<').unwrap_or(alt);
    ["DEL", "INS", "DUP", "INV", "CNV"]
        .iter()
        .any(|kind| alt.starts_with(kind))
        || alt.strip_prefix("CN").is_some_and(starts_with_copy_number)
}

fn starts_with_copy_number(text: &str) -> bool {
    let text = text.strip_prefix('=').unwrap_or(text);
    text.starts_with(|c: char| c.is_ascii_digit())
}

fn without_repeats(kind: &str) -> &str {
    let normalized = |part: &str| {
        part.trim()
            .trim_start_matches('<')
            .trim_end_matches('>')
            .to_ascii_uppercase()
    };
    let mut parts = kind.split(['/', ',']);
    let first = parts.next().unwrap_or(kind);
    if parts.all(|part| normalized(part) == normalized(first)) {
        first
    } else {
        kind
    }
}

fn mobile_element(upper: &str) -> Option<String> {
    let (start, action) = ["INS", "DEL"]
        .into_iter()
        .filter_map(|action| Some((upper.find(&format!("{action}:ME"))?, action)))
        .min()?;
    let rest = &upper[start + "INS:ME".len()..];
    let rest = rest.strip_prefix(':').unwrap_or(rest);
    let end = rest
        .find(|c: char| !c.is_ascii_alphanumeric())
        .unwrap_or(rest.len());
    let subtype = match &rest[..end] {
        "L1" | "LINE1" => "LINE1",
        element @ ("ALU" | "HERV" | "SVA") => element,
        _ => "ME",
    };
    Some(format!("{action}_{subtype}"))
}

fn has_copy_number(upper: &str) -> bool {
    upper
        .match_indices("CN")
        .any(|(start, _)| starts_with_copy_number(&upper[start + "CN".len()..]))
}

fn abbreviation(kind: &str) -> String {
    let upper = kind.to_ascii_uppercase();
    if let Some(mobile) = mobile_element(&upper) {
        return mobile;
    }
    let abbreviation = if upper.contains("DEL") && upper.contains("DUP") {
        "CNV"
    } else if upper.contains("DUP:TANDEM") {
        "TDUP"
    } else if upper.contains("CNV:TR") {
        "TREP"
    } else if has_copy_number(&upper) {
        match kind.trim_start_matches('<').trim_end_matches('>') {
            "CN0" | "CN=0" => "DEL",
            "CN2" | "CN=2" => "DUP",
            _ => "CNV",
        }
    } else if upper.contains("CNV") {
        "CNV"
    } else if is_breakend(kind) {
        "BND"
    } else if kind.starts_with('<') || kind.ends_with('>') {
        let bare = kind.replace(['<', '>'], "");
        return bare.split(':').next().unwrap_or_default().to_string();
    } else {
        kind
    };
    abbreviation.to_string()
}

fn so_term(kind: &str) -> Option<&'static str> {
    Some(match abbreviation(without_repeats(kind)).as_str() {
        "INS" => "insertion",
        "INS_ME" => "mobile_element_insertion",
        "INS_ALU" => "Alu_insertion",
        "INS_HERV" => "HERV_insertion",
        "INS_LINE1" => "LINE1_insertion",
        "INS_SVA" => "SVA_insertion",
        "DEL" => "deletion",
        "DEL_ME" => "mobile_element_deletion",
        "DEL_ALU" => "Alu_deletion",
        "DEL_HERV" => "HERV_deletion",
        "DEL_LINE1" => "LINE1_deletion",
        "DEL_SVA" => "SVA_deletion",
        "TREP" => "tandem_repeat",
        "TDUP" => "tandem_duplication",
        "DUP" => "duplication",
        "CNV" => "copy_number_variation",
        "INV" => "inversion",
        "BND" => "chromosome_breakpoint",
        _ => return None,
    })
}

fn sequence_alleles(reference: &str, alternates: &[String]) -> Vec<String> {
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

pub fn is_breakend(alt: &str) -> bool {
    alt.contains(['[', ']']) || alt.starts_with('.') || alt.ends_with('.')
}

pub fn is_structural(alt: &str) -> bool {
    is_breakend(alt) || alt.contains('<')
}

pub fn is_skipped_alt(alt: &str) -> bool {
    matches!(alt, "*" | "<*>" | "<NON_REF>")
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
    #[case("A", &["<DEL>"], &["deletion"])]
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
        assert_eq!(
            vep_alleles(reference, &owned(alternates), None),
            owned(expected)
        );
    }

    #[rstest]
    #[case("<DEL>", "deletion")]
    #[case("<DEL:ME>", "mobile_element_deletion")]
    #[case("<DEL:ME:ALU>", "Alu_deletion")]
    #[case("<DEL:ME:HERV>", "HERV_deletion")]
    #[case("<DEL:ME:L1>", "LINE1_deletion")]
    #[case("<DEL:ME:SVA>", "SVA_deletion")]
    #[case("<INS>", "insertion")]
    #[case("<INS:ME>", "mobile_element_insertion")]
    #[case("<INS:ME:ALU>", "Alu_insertion")]
    #[case("<INS:ME:HERV>", "HERV_insertion")]
    #[case("<INS:ME:LINE1>", "LINE1_insertion")]
    #[case("<INS:ME:L1>", "LINE1_insertion")]
    #[case("<INS:ME:SVA>", "SVA_insertion")]
    #[case("<INS:ME:OTHER>", "mobile_element_insertion")]
    #[case("<DUP>", "duplication")]
    #[case("<DUP:TANDEM>", "tandem_duplication")]
    #[case("<INV>", "inversion")]
    #[case("<INV:FLIP>", "inversion")]
    #[case("<CNV>", "copy_number_variation")]
    #[case("<CNV:TR>", "tandem_repeat")]
    #[case("<CN0>", "deletion")]
    #[case("<CN=0>", "deletion")]
    #[case("<CN2>", "duplication")]
    #[case("<CN3>", "copy_number_variation")]
    #[case("<BND>", "chromosome_breakpoint")]
    fn a_symbolic_alt_is_spelled_as_its_so_term(#[case] alt: &str, #[case] term: &str) {
        assert_eq!(vep_alleles("A", &owned(&[alt]), None), owned(&[term]));
    }

    #[rstest]
    #[case(&["<DEL>", "<DUP>"], "copy_number_variation")]
    #[case(&["<CN0>", "<CN2>"], "copy_number_variation")]
    #[case(&["<INV>", "<inv>"], "inversion")]
    #[case(&["<INS:ME:L1>", "<INS:ME:ALU>"], "LINE1_insertion")]
    fn symbolic_alts_of_one_record_share_one_term(#[case] alts: &[&str], #[case] term: &str) {
        let expected = vec![term.to_string(); alts.len()];
        assert_eq!(vep_alleles("A", &owned(alts), None), expected);
    }

    #[rstest]
    #[case(&["<VNTR>"])]
    #[case(&["<DEL>", "<INV>"])]
    #[case(&["AT", "<NON_REF>"])]
    fn a_record_vep_has_no_term_for_keeps_its_alts(#[case] alts: &[&str]) {
        assert_eq!(vep_alleles("A", &owned(alts), None), owned(alts));
    }

    #[rstest]
    #[case("C[chr17:7680000[")]
    #[case("]chr17:7680000]C")]
    #[case("CT[chr17:7680000[")]
    #[case("C.")]
    #[case(".C")]
    fn a_breakend_is_never_trimmed(#[case] alt: &str) {
        assert_eq!(vep_alleles("C", &owned(&[alt]), None), owned(&[alt]));
        assert_eq!(vep_alleles("C", &owned(&[alt, "CT"]), None)[0], alt);
    }

    #[rstest]
    #[case("A", &["ACN"], "INS", &["insertion"])]
    #[case("A", &["<STR4>"], "CNV", &["copy_number_variation"])]
    #[case("A", &["<DUP>"], "DEL", &["duplication"])]
    #[case("A", &["AT"], "INS", &["T"])]
    fn an_svtype_routes_a_record_with_non_acgt_alleles(
        #[case] reference: &str,
        #[case] alternates: &[&str],
        #[case] svtype: &str,
        #[case] expected: &[&str],
    ) {
        assert_eq!(
            vep_alleles(reference, &owned(alternates), Some(svtype)),
            owned(expected)
        );
    }

    #[rstest]
    #[case("<DEL>", true)]
    #[case("<*>", true)]
    #[case("C[chr1:5[", true)]
    #[case("]chr1:5]C", true)]
    #[case("C.", true)]
    #[case(".C", true)]
    #[case("CT", false)]
    #[case("*", false)]
    fn structural_alts_are_symbolic_or_breakends(#[case] alt: &str, #[case] expected: bool) {
        assert_eq!(is_structural(alt), expected);
    }

    #[test]
    fn skipped_alts_are_placeholders_and_spanning_deletions() {
        assert!(is_skipped_alt("<*>"));
        assert!(is_skipped_alt("<NON_REF>"));
        assert!(is_skipped_alt("*"));
        assert!(!is_skipped_alt("<DEL>"));
    }

    #[test]
    fn variant_name_joins_with_colons() {
        assert_eq!(
            variant_name("chr17", 7_673_704, "G", "A"),
            "chr17:7673704:G:A"
        );
    }
}
