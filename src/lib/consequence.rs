pub const VEP_RELEASE: u32 = 116;

pub const UNRANKED: u8 = 42;

const CONSEQUENCES: [&str; 41] = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
];

const LOSS_OF_FUNCTION: [&str; 7] = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
];

const CODING: [&str; 13] = [
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
];

pub fn rank(term: &str) -> Option<u8> {
    let position = CONSEQUENCES.iter().position(|&known| known == term)?;
    u8::try_from(position + 1).ok()
}

pub fn most_severe(consequence: &str) -> (&str, u8) {
    consequence
        .split('&')
        .map(|term| (term, rank(term).unwrap_or(UNRANKED)))
        .min_by_key(|&(_, rank)| rank)
        .expect("split yields at least one term")
}

pub fn vep_major(version_line: &str) -> Option<u32> {
    let version = version_line.trim_start_matches('"').strip_prefix('v')?;
    let end = version
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(version.len());
    version[..end].parse().ok()
}

pub fn is_coding(term: &str) -> bool {
    CODING.contains(&term)
}

pub fn is_lof(term: &str) -> bool {
    LOSS_OF_FUNCTION.contains(&term)
}

pub fn impact_rank(impact: &str) -> Option<u8> {
    match impact {
        "HIGH" => Some(0),
        "MODERATE" => Some(1),
        "LOW" => Some(2),
        "MODIFIER" => Some(3),
        _ => None,
    }
}

pub fn is_up_or_downstream(term: &str) -> bool {
    matches!(term, "upstream_gene_variant" | "downstream_gene_variant")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_follows_the_vep_116_order() {
        assert_eq!(rank("transcript_ablation"), Some(1));
        assert_eq!(rank("stop_gained"), Some(4));
        assert_eq!(rank("missense_variant"), Some(13));
        assert_eq!(rank("upstream_gene_variant"), Some(32));
        assert_eq!(rank("sequence_variant"), Some(41));
        assert_eq!(rank("made_up"), None);
    }

    #[test]
    fn most_severe_picks_the_lowest_rank() {
        assert_eq!(
            most_severe("intron_variant&splice_region_variant"),
            ("splice_region_variant", 16)
        );
    }

    #[test]
    fn an_unranked_term_sorts_after_every_known_term() {
        assert_eq!(most_severe("made_up"), ("made_up", UNRANKED));
        assert_eq!(
            most_severe("made_up&sequence_variant"),
            ("sequence_variant", 41)
        );
        assert!(CONSEQUENCES.iter().all(|term| rank(term) < Some(UNRANKED)));
    }

    #[test]
    fn an_unranked_term_is_neither_coding_nor_lof() {
        assert!(!is_coding("made_up") && !is_lof("made_up"));
    }

    #[test]
    fn the_major_version_comes_from_the_vep_header_value() {
        assert_eq!(vep_major("\"v116.0\" API=\"v116\" time=\"now\""), Some(116));
        assert_eq!(vep_major("\"v116\" time=\"now\""), Some(116));
        assert_eq!(vep_major("\"v115.2\""), Some(115));
        assert_eq!(vep_major("v99"), Some(99));
        assert_eq!(vep_major("\"release\""), None);
        assert_eq!(vep_major("\"v\""), None);
    }

    #[test]
    fn lof_terms_are_identified() {
        assert!(is_lof("frameshift_variant"));
        assert!(!is_lof("missense_variant"));
    }

    #[test]
    fn coding_terms_are_identified() {
        assert!(is_coding("synonymous_variant"));
        assert!(!is_coding("intron_variant"));
    }

    #[test]
    fn impact_ranks_follow_severity() {
        assert_eq!(impact_rank("HIGH"), Some(0));
        assert_eq!(impact_rank("MODERATE"), Some(1));
        assert_eq!(impact_rank("LOW"), Some(2));
        assert_eq!(impact_rank("MODIFIER"), Some(3));
        assert_eq!(impact_rank("UNKNOWN"), None);
    }

    #[test]
    fn up_and_downstream_terms_are_identified() {
        assert!(is_up_or_downstream("upstream_gene_variant"));
        assert!(is_up_or_downstream("downstream_gene_variant"));
        assert!(!is_up_or_downstream("intergenic_variant"));
    }
}
