use anyhow::anyhow;

pub const VEP_RELEASE: u32 = 116;

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

pub fn most_severe(consequence: &str) -> anyhow::Result<(&'static str, u8)> {
    let ranks = consequence
        .split('&')
        .map(|term| {
            rank(term).ok_or_else(|| {
                anyhow!("VEP {VEP_RELEASE} does not rank the Consequence term '{term}'")
            })
        })
        .collect::<anyhow::Result<Vec<u8>>>()?;
    let best = ranks
        .into_iter()
        .min()
        .expect("split yields at least one term");
    Ok((CONSEQUENCES[usize::from(best) - 1], best))
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
        let (term, rank) = most_severe("intron_variant&splice_region_variant").unwrap();
        assert_eq!(term, "splice_region_variant");
        assert_eq!(rank, 16);
    }

    #[test]
    fn most_severe_errors_on_an_unranked_term() {
        let error = most_severe("made_up").unwrap_err();
        assert_eq!(
            error.to_string(),
            "VEP 116 does not rank the Consequence term 'made_up'"
        );
    }

    #[test]
    fn most_severe_errors_even_when_another_term_ranks() {
        assert!(most_severe("missense_variant&made_up").is_err());
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
