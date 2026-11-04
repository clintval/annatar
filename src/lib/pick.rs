use std::cmp::Reverse;
use std::collections::HashSet;
use std::fmt::Write;
use std::io::Read;
use std::path::Path;

use anyhow::{Context, anyhow, bail};
use flate2::read::MultiGzDecoder;
use sha2::{Digest, Sha256};

use crate::consequence::{impact_rank, is_up_or_downstream, most_severe};
use crate::csq::{CsqEntry, CsqIndex};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

#[derive(Clone, Debug)]
pub struct Hashed<T> {
    pub value: T,
    pub sha256: String,
}

fn sha256_hex(bytes: &[u8]) -> String {
    let mut hex = String::new();
    for byte in Sha256::digest(bytes) {
        write!(hex, "{byte:02x}").expect("writing to a String cannot fail");
    }
    hex
}

#[derive(Clone, Debug, Default)]
pub struct Genes(HashSet<String>);

impl Genes {
    pub fn contains(&self, name: &str) -> bool {
        self.0.contains(name)
    }

    pub fn sorted(&self) -> Vec<String> {
        let mut names: Vec<String> = self.0.iter().cloned().collect();
        names.sort_unstable();
        names
    }
}

pub fn read_genes(path: &Path) -> anyhow::Result<Hashed<Genes>> {
    let bytes = std::fs::read(path)
        .with_context(|| format!("cannot read the genes file {}", path.display()))?;
    let text = std::str::from_utf8(&bytes)
        .with_context(|| format!("the genes file {} is not text", path.display()))?;
    let names = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    Ok(Hashed {
        value: Genes(names),
        sha256: sha256_hex(&bytes),
    })
}

#[derive(Clone, Debug, Default)]
pub struct ManeTable {
    pub select: HashSet<String>,
    pub plus_clinical: HashSet<String>,
}

fn without_version(id: &str) -> &str {
    id.split('.').next().unwrap_or(id)
}

pub fn read_mane_table(path: &Path) -> anyhow::Result<Hashed<ManeTable>> {
    let bytes = std::fs::read(path)
        .with_context(|| format!("cannot read the MANE table {}", path.display()))?;
    let sha256 = sha256_hex(&bytes);
    let text = if bytes.starts_with(&GZIP_MAGIC) {
        let mut text = String::new();
        MultiGzDecoder::new(bytes.as_slice())
            .read_to_string(&mut text)
            .with_context(|| format!("cannot decompress the MANE table {}", path.display()))?;
        text
    } else {
        String::from_utf8(bytes)
            .with_context(|| format!("the MANE table {} is not text", path.display()))?
    };

    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow!("the MANE table {} is empty", path.display()))?;
    let columns: Vec<&str> = header
        .split('\t')
        .map(|name| name.strip_prefix('#').unwrap_or(name))
        .collect();
    let column = |name: &str| {
        columns
            .iter()
            .position(|&c| c == name)
            .ok_or_else(|| anyhow!("the MANE table {} has no {name} column", path.display()))
    };
    let transcript = column("Ensembl_nuc")?;
    let status = column("MANE_status")?;

    let mut table = ManeTable::default();
    for line in lines.filter(|line| !line.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(&id), Some(&kind)) = (fields.get(transcript), fields.get(status)) else {
            bail!("a MANE table row in {} has too few columns", path.display());
        };
        let id = without_version(id).to_string();
        match kind {
            "MANE Select" => table.select.insert(id),
            "MANE Plus Clinical" => table.plus_clinical.insert(id),
            _ => false,
        };
    }
    Ok(Hashed {
        value: table,
        sha256,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tier {
    ManeSelect,
    ManePlusClinical,
    Canonical,
    Other,
}

impl Tier {
    pub fn as_str(self) -> &'static str {
        match self {
            Tier::ManeSelect => "mane_select",
            Tier::ManePlusClinical => "mane_plus_clinical",
            Tier::Canonical => "canonical",
            Tier::Other => "other",
        }
    }
}

pub struct PickContext<'a> {
    pub index: &'a CsqIndex,
    pub genes: Option<&'a Genes>,
    pub mane: Option<&'a ManeTable>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Picked {
    pub entry: usize,
    pub tier: Tier,
}

struct Scored<'a> {
    position: usize,
    entry: &'a CsqEntry,
    rank: u8,
    inside: bool,
    mane: Option<Tier>,
    tier: Tier,
    canonical: bool,
}

fn has_gene(entry: &CsqEntry, index: &CsqIndex) -> bool {
    entry.get(index, "Gene").is_some() || entry.get(index, "SYMBOL").is_some()
}

fn in_panel(entry: &CsqEntry, index: &CsqIndex, genes: &Genes) -> bool {
    ["SYMBOL", "Gene"].iter().any(|name| {
        entry
            .get(index, name)
            .is_some_and(|gene| genes.contains(gene))
    })
}

fn mane_tier(entry: &CsqEntry, ctx: &PickContext) -> Option<Tier> {
    let declared = ctx.index.position("MANE_SELECT").is_some()
        || ctx.index.position("MANE_PLUS_CLINICAL").is_some();
    let (select, plus_clinical) = if declared {
        (
            entry.get(ctx.index, "MANE_SELECT").is_some(),
            entry.get(ctx.index, "MANE_PLUS_CLINICAL").is_some(),
        )
    } else if let (Some(table), Some(feature)) = (ctx.mane, entry.get(ctx.index, "Feature")) {
        let id = without_version(feature);
        (table.select.contains(id), table.plus_clinical.contains(id))
    } else {
        (false, false)
    };
    if select {
        Some(Tier::ManeSelect)
    } else if plus_clinical {
        Some(Tier::ManePlusClinical)
    } else {
        None
    }
}

fn score<'a>(
    position: usize,
    entry: &'a CsqEntry,
    ctx: &PickContext,
) -> anyhow::Result<Scored<'a>> {
    let consequence = entry
        .get(ctx.index, "Consequence")
        .ok_or_else(|| anyhow!("a CSQ entry has no Consequence"))?;
    let (term, rank) = most_severe(consequence)?;
    let mane = mane_tier(entry, ctx);
    let canonical = entry.get(ctx.index, "CANONICAL") == Some("YES");
    let tier = mane.unwrap_or(if canonical {
        Tier::Canonical
    } else {
        Tier::Other
    });
    Ok(Scored {
        position,
        entry,
        rank,
        inside: has_gene(entry, ctx.index) && !is_up_or_downstream(term),
        mane,
        tier,
        canonical,
    })
}

fn narrow<'a, 'b>(
    pool: Vec<&'b Scored<'a>>,
    keep: impl Fn(&Scored<'a>) -> bool,
) -> Vec<&'b Scored<'a>> {
    let kept: Vec<_> = pool.iter().copied().filter(|s| keep(s)).collect();
    if kept.is_empty() { pool } else { kept }
}

fn transcript_number(feature: Option<&str>) -> u64 {
    let id = feature.map(without_version).unwrap_or_default();
    let digits = &id[id.trim_end_matches(|c: char| c.is_ascii_digit()).len()..];
    digits.parse().unwrap_or(u64::MAX)
}

fn protein_length(entry: &CsqEntry, index: &CsqIndex) -> u64 {
    entry
        .get(index, "Protein_position")
        .and_then(|position| position.split_once('/'))
        .and_then(|(_, total)| total.parse().ok())
        .unwrap_or(0)
}

pub fn entries_for_allele(entries: &[CsqEntry], index: &CsqIndex, allele: &str) -> Vec<usize> {
    (0..entries.len())
        .filter(|&i| entries[i].get(index, "Allele") == Some(allele))
        .collect()
}

pub fn pick(
    entries: &[CsqEntry],
    candidates: &[usize],
    ctx: &PickContext,
) -> anyhow::Result<Option<Picked>> {
    let scored = candidates
        .iter()
        .map(|&i| score(i, &entries[i], ctx))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut pool: Vec<&Scored> = scored.iter().collect();
    if pool.is_empty() {
        return Ok(None);
    }

    if pool.iter().any(|s| s.inside) {
        pool.retain(|s| s.inside);
        if let Some(genes) = ctx.genes {
            pool = narrow(pool, |s| in_panel(s.entry, ctx.index, genes));
        }
    }

    if pool.iter().any(|s| s.mane.is_some()) {
        pool.retain(|s| s.mane.is_some());
    } else {
        pool = narrow(pool, |s| s.canonical);
    }

    let best_rank = pool
        .iter()
        .map(|s| s.rank)
        .min()
        .expect("pool is not empty");
    pool.retain(|s| s.rank == best_rank);

    let key = |s: &&Scored| {
        let feature = s.entry.get(ctx.index, "Feature");
        (
            s.mane != Some(Tier::ManeSelect),
            !s.canonical,
            s.entry.get(ctx.index, "CCDS").is_none(),
            s.entry.get(ctx.index, "SWISSPROT").is_none(),
            Reverse(protein_length(s.entry, ctx.index)),
            transcript_number(feature),
            feature.is_none(),
            feature.unwrap_or_default().to_string(),
        )
    };
    let winner = pool.into_iter().min_by_key(key).expect("pool is not empty");
    Ok(Some(Picked {
        entry: winner.position,
        tier: winner.tier,
    }))
}

pub fn worst_panel_impact(
    entries: &[CsqEntry],
    candidates: &[usize],
    ctx: &PickContext,
) -> Option<&'static str> {
    const IMPACTS: [&str; 4] = ["HIGH", "MODERATE", "LOW", "MODIFIER"];
    let genes = ctx.genes?;
    candidates
        .iter()
        .map(|&i| &entries[i])
        .filter(|entry| in_panel(entry, ctx.index, genes))
        .filter(|entry| {
            entry
                .get(ctx.index, "Consequence")
                .and_then(|consequence| most_severe(consequence).ok())
                .is_some_and(|(term, _)| !is_up_or_downstream(term))
        })
        .filter_map(|entry| impact_rank(entry.get(ctx.index, "IMPACT")?))
        .min()
        .map(|rank| IMPACTS[usize::from(rank)])
}

pub fn missing_tie_fields(index: &CsqIndex) -> Vec<&'static str> {
    [
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "CANONICAL",
        "CCDS",
        "SWISSPROT",
        "Protein_position",
        "SYMBOL",
        "Gene",
        "Feature",
    ]
    .into_iter()
    .filter(|name| index.position(name).is_none())
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const NAMES: [&str; 13] = [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature",
        "CANONICAL",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "CCDS",
        "SWISSPROT",
        "Protein_position",
        "BIOTYPE",
    ];

    fn index_without(dropped: &[&str]) -> CsqIndex {
        CsqIndex {
            names: NAMES
                .iter()
                .filter(|name| !dropped.contains(name))
                .map(ToString::to_string)
                .collect(),
        }
    }

    fn entry_for(index: &CsqIndex, fields: &[(&str, &str)]) -> CsqEntry {
        let mut values = vec![String::new(); index.names.len()];
        values[index.position("Allele").unwrap()] = "T".to_string();
        for (name, value) in fields {
            if let Some(position) = index.position(name) {
                values[position] = (*value).to_string();
            }
        }
        CsqEntry { values }
    }

    fn entries(index: &CsqIndex, rows: &[&[(&str, &str)]]) -> Vec<CsqEntry> {
        rows.iter().map(|row| entry_for(index, row)).collect()
    }

    fn all(entries: &[CsqEntry]) -> Vec<usize> {
        (0..entries.len()).collect()
    }

    fn genes(symbols: &[&str]) -> Genes {
        Genes(symbols.iter().map(ToString::to_string).collect())
    }

    fn table(select: &[&str], plus_clinical: &[&str]) -> ManeTable {
        ManeTable {
            select: select.iter().map(ToString::to_string).collect(),
            plus_clinical: plus_clinical.iter().map(ToString::to_string).collect(),
        }
    }

    fn pick_with(
        index: &CsqIndex,
        entries: &[CsqEntry],
        genes: Option<&Genes>,
        mane: Option<&ManeTable>,
    ) -> Option<Picked> {
        let ctx = PickContext { index, genes, mane };
        pick(entries, &all(entries), &ctx).unwrap()
    }

    #[test]
    fn prefers_panel_gene() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "intron_variant"), ("SYMBOL", "GENE_A")],
                &[("Consequence", "intron_variant"), ("SYMBOL", "TP53")],
            ],
        );
        let panel = genes(&["TP53"]);
        let picked = pick_with(&index, &rows, Some(&panel), None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn panel_matches_the_ensembl_gene_id() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "missense_variant"),
                    ("SYMBOL", "GENE_A"),
                    ("Gene", "ENSG1"),
                ],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "GENE_B"),
                    ("Gene", "ENSG2"),
                ],
            ],
        );
        let panel = genes(&["ENSG2"]);
        let picked = pick_with(&index, &rows, Some(&panel), None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn falls_back_to_every_gene_when_no_panel_gene_is_hit() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "intron_variant"), ("SYMBOL", "GENE_A")],
                &[("Consequence", "missense_variant"), ("SYMBOL", "GENE_B")],
            ],
        );
        let panel = genes(&["TP53"]);
        let picked = pick_with(&index, &rows, Some(&panel), None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn ignores_up_and_downstream_genes() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "upstream_gene_variant"),
                    ("SYMBOL", "WRAP53"),
                ],
                &[("Consequence", "missense_variant"), ("SYMBOL", "TP53")],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn a_panel_gene_only_up_or_downstream_is_not_inside() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "downstream_gene_variant"),
                    ("SYMBOL", "TP53"),
                ],
                &[("Consequence", "intron_variant"), ("SYMBOL", "GENE_A")],
            ],
        );
        let panel = genes(&["TP53"]);
        let picked = pick_with(&index, &rows, Some(&panel), None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn regulatory_entry_without_a_gene_is_not_inside_a_gene() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "regulatory_region_variant")],
                &[("Consequence", "upstream_gene_variant"), ("SYMBOL", "TP53")],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn no_gene_uses_all_entries() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "intergenic_variant")],
                &[("Consequence", "regulatory_region_variant")],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn upstream_only_entries_still_pick_one() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "downstream_gene_variant"), ("SYMBOL", "A")],
                &[("Consequence", "upstream_gene_variant"), ("SYMBOL", "B")],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn mane_before_severity() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST2"),
                ],
                &[
                    ("Consequence", "synonymous_variant"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST1"),
                    ("MANE_SELECT", "NM_000546.6"),
                ],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 1,
                tier: Tier::ManeSelect
            }
        );
    }

    #[test]
    fn mane_plus_clinical_missense_beats_mane_select_synonymous() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "synonymous_variant"),
                    ("SYMBOL", "TTN"),
                    ("MANE_SELECT", "NM_1.1"),
                ],
                &[
                    ("Consequence", "missense_variant"),
                    ("SYMBOL", "TTN"),
                    ("MANE_PLUS_CLINICAL", "NM_2.1"),
                ],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 1,
                tier: Tier::ManePlusClinical
            }
        );
    }

    #[test]
    fn canonical_when_no_mane() {
        let index = index_without(&["MANE_SELECT", "MANE_PLUS_CLINICAL"]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "stop_gained"), ("SYMBOL", "TP53")],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("CANONICAL", "YES"),
                ],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 1,
                tier: Tier::Canonical
            }
        );
    }

    #[test]
    fn tier_is_other_without_any_flag() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[&[("Consequence", "intron_variant"), ("SYMBOL", "TP53")]],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.tier, Tier::Other);
    }

    #[test]
    fn mane_table_fallback() {
        let index = index_without(&["MANE_SELECT", "MANE_PLUS_CLINICAL"]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST2.1"),
                ],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST1.4"),
                ],
            ],
        );
        let mane = table(&["ENST1"], &[]);
        let picked = pick_with(&index, &rows, None, Some(&mane)).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 1,
                tier: Tier::ManeSelect
            }
        );
    }

    #[test]
    fn mane_table_plus_clinical_tier() {
        let index = index_without(&["MANE_SELECT", "MANE_PLUS_CLINICAL"]);
        let rows = entries(
            &index,
            &[&[
                ("Consequence", "intron_variant"),
                ("SYMBOL", "TTN"),
                ("Feature", "ENST9"),
            ]],
        );
        let mane = table(&[], &["ENST9"]);
        let picked = pick_with(&index, &rows, None, Some(&mane)).unwrap();
        assert_eq!(picked.tier, Tier::ManePlusClinical);
    }

    #[test]
    fn mane_table_is_ignored_when_the_header_declares_mane_subfields() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST2"),
                ],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST1"),
                ],
            ],
        );
        let mane = table(&["ENST1"], &[]);
        let picked = pick_with(&index, &rows, None, Some(&mane)).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 0,
                tier: Tier::Other
            }
        );
    }

    fn picked_among_tied(rows: &[&[(&str, &str)]]) -> usize {
        let index = index_without(&[]);
        let full: Vec<Vec<(&str, &str)>> = rows
            .iter()
            .map(|extra| {
                let mut row = vec![("Consequence", "missense_variant"), ("SYMBOL", "TP53")];
                row.extend_from_slice(extra);
                row
            })
            .collect();
        let slices: Vec<&[(&str, &str)]> = full.iter().map(Vec::as_slice).collect();
        let built = entries(&index, &slices);
        pick_with(&index, &built, None, None).unwrap().entry
    }

    #[test]
    fn tie_mane_select_beats_a_canonical_plus_clinical_on_a_lower_id() {
        let picked = picked_among_tied(&[
            &[
                ("Feature", "ENST1"),
                ("MANE_PLUS_CLINICAL", "NM_1.1"),
                ("CANONICAL", "YES"),
            ],
            &[("Feature", "ENST2"), ("MANE_SELECT", "NM_2.1")],
        ]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn tie_canonical_beats_non_canonical_among_mane_entries() {
        let picked = picked_among_tied(&[
            &[("Feature", "ENST1"), ("MANE_PLUS_CLINICAL", "NM_1.1")],
            &[
                ("Feature", "ENST2"),
                ("MANE_PLUS_CLINICAL", "NM_2.1"),
                ("CANONICAL", "YES"),
            ],
        ]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn tie_ccds_beats_swissprot_and_a_longer_protein() {
        let picked = picked_among_tied(&[
            &[
                ("Feature", "ENST1"),
                ("SWISSPROT", "P04637"),
                ("Protein_position", "1/999"),
            ],
            &[("Feature", "ENST2"), ("CCDS", "CCDS11118")],
        ]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn tie_swissprot_beats_a_longer_protein() {
        let picked = picked_among_tied(&[
            &[("Feature", "ENST1"), ("Protein_position", "1/999")],
            &[("Feature", "ENST2"), ("SWISSPROT", "P04637")],
        ]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn tie_longest_protein_reads_the_length_after_the_slash() {
        let picked = picked_among_tied(&[
            &[("Feature", "ENST1"), ("Protein_position", "306/341")],
            &[("Feature", "ENST2"), ("Protein_position", "306-307/393")],
        ]);
        assert_eq!(picked, 1);
        let picked = picked_among_tied(&[
            &[("Feature", "ENST2"), ("Protein_position", "-/393")],
            &[("Feature", "ENST1"), ("Protein_position", "306")],
        ]);
        assert_eq!(picked, 0);
    }

    #[test]
    fn tie_transcript_number_is_numeric_and_ignores_the_version() {
        let picked = picked_among_tied(&[&[("Feature", "ENST10")], &[("Feature", "ENST9")]]);
        assert_eq!(picked, 1);
        let picked = picked_among_tied(&[
            &[("Feature", "ENST00000000002.5")],
            &[("Feature", "ENST00000000001.9")],
        ]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn tie_without_a_feature_sorts_last() {
        let picked = picked_among_tied(&[&[], &[("Feature", "ENST00000000009")]]);
        assert_eq!(picked, 1);
    }

    #[test]
    fn exact_ties_fall_to_the_lowest_transcript_id() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "missense_variant"),
                    ("SYMBOL", "A"),
                    ("Feature", "ENSTB00001"),
                ],
                &[
                    ("Consequence", "missense_variant"),
                    ("SYMBOL", "A"),
                    ("Feature", "ENSTA00001"),
                ],
            ],
        );
        let picked = pick_with(&index, &rows, None, None).unwrap();
        assert_eq!(picked.entry, 1);
    }

    #[test]
    fn missing_star_entry_is_none() {
        let index = index_without(&[]);
        let rows = entries(&index, &[&[("Consequence", "intron_variant")]]);
        let found = entries_for_allele(&rows, &index, "*");
        assert_eq!(found.len(), 0);
        let ctx = PickContext {
            index: &index,
            genes: None,
            mane: None,
        };
        assert_eq!(pick(&rows, &found, &ctx).unwrap(), None);
    }

    #[test]
    fn entries_for_allele_matches_the_allele_subfield() {
        let index = index_without(&[]);
        let mut rows = entries(&index, &[&[("Consequence", "intron_variant")][..]; 3]);
        rows[1].values[0] = "-".to_string();
        assert_eq!(entries_for_allele(&rows, &index, "T"), [0, 2]);
        assert_eq!(entries_for_allele(&rows, &index, "-"), [1]);
    }

    #[test]
    fn missing_consequence_is_an_error() {
        let index = index_without(&[]);
        let rows = entries(&index, &[&[("SYMBOL", "TP53")]]);
        let ctx = PickContext {
            index: &index,
            genes: None,
            mane: None,
        };
        let error = pick(&rows, &[0], &ctx).unwrap_err();
        assert_eq!(error.to_string(), "a CSQ entry has no Consequence");
    }

    #[test]
    fn unranked_term_is_an_error() {
        let index = index_without(&[]);
        let rows = entries(&index, &[&[("Consequence", "made_up")]]);
        let ctx = PickContext {
            index: &index,
            genes: None,
            mane: None,
        };
        assert!(pick(&rows, &[0], &ctx).is_err());
    }

    #[test]
    fn worst_panel_impact_any_transcript() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("IMPACT", "MODIFIER"),
                ],
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("IMPACT", "HIGH"),
                ],
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "OTHER"),
                    ("IMPACT", "HIGH"),
                ],
            ],
        );
        let panel = genes(&["TP53"]);
        let ctx = PickContext {
            index: &index,
            genes: Some(&panel),
            mane: None,
        };
        assert_eq!(worst_panel_impact(&rows, &[0, 1], &ctx), Some("HIGH"));
        assert_eq!(worst_panel_impact(&rows, &[0, 2], &ctx), Some("MODIFIER"));
    }

    #[test]
    fn worst_panel_impact_skips_up_and_downstream_and_needs_a_panel() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[&[
                ("Consequence", "upstream_gene_variant"),
                ("SYMBOL", "TP53"),
                ("IMPACT", "MODIFIER"),
            ]],
        );
        let panel = genes(&["TP53"]);
        let ctx = PickContext {
            index: &index,
            genes: Some(&panel),
            mane: None,
        };
        assert_eq!(worst_panel_impact(&rows, &[0], &ctx), None);
        let ctx = PickContext {
            index: &index,
            genes: None,
            mane: None,
        };
        assert_eq!(worst_panel_impact(&rows, &[0], &ctx), None);
    }

    #[test]
    fn worst_panel_impact_ignores_an_unranked_impact() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("IMPACT", "WEIRD"),
                ],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("IMPACT", "LOW"),
                ],
            ],
        );
        let panel = genes(&["TP53"]);
        let ctx = PickContext {
            index: &index,
            genes: Some(&panel),
            mane: None,
        };
        assert_eq!(worst_panel_impact(&rows, &[0, 1], &ctx), Some("LOW"));
    }

    #[test]
    fn missing_tie_fields_lists_undeclared_subfields() {
        assert_eq!(missing_tie_fields(&index_without(&[])).len(), 0);
        let index = index_without(&["CCDS", "Protein_position", "MANE_SELECT", "Gene", "Feature"]);
        assert_eq!(
            missing_tie_fields(&index),
            ["MANE_SELECT", "CCDS", "Protein_position", "Gene", "Feature"]
        );
    }

    #[test]
    fn worst_panel_impact_counts_an_unranked_term_as_not_inside() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[&[
                ("Consequence", "made_up"),
                ("SYMBOL", "TP53"),
                ("IMPACT", "HIGH"),
            ]],
        );
        let panel = genes(&["TP53"]);
        let ctx = PickContext {
            index: &index,
            genes: Some(&panel),
            mane: None,
        };
        assert_eq!(worst_panel_impact(&rows, &[0], &ctx), None);
    }

    #[test]
    fn mane_table_is_ignored_when_only_one_mane_subfield_is_declared() {
        let index = index_without(&["MANE_PLUS_CLINICAL"]);
        let rows = entries(
            &index,
            &[
                &[
                    ("Consequence", "stop_gained"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST2"),
                ],
                &[
                    ("Consequence", "intron_variant"),
                    ("SYMBOL", "TP53"),
                    ("Feature", "ENST1"),
                ],
            ],
        );
        let mane = table(&["ENST1"], &[]);
        let picked = pick_with(&index, &rows, None, Some(&mane)).unwrap();
        assert_eq!(
            picked,
            Picked {
                entry: 0,
                tier: Tier::Other
            }
        );
    }

    #[test]
    fn tier_strings_are_stable() {
        assert_eq!(Tier::ManeSelect.as_str(), "mane_select");
        assert_eq!(Tier::ManePlusClinical.as_str(), "mane_plus_clinical");
        assert_eq!(Tier::Canonical.as_str(), "canonical");
        assert_eq!(Tier::Other.as_str(), "other");
    }

    #[test]
    fn picked_entry_indexes_into_all_entries_when_candidates_are_a_subset() {
        let index = index_without(&[]);
        let rows = entries(
            &index,
            &[
                &[("Consequence", "stop_gained"), ("SYMBOL", "TP53")],
                &[("Consequence", "intron_variant"), ("SYMBOL", "TP53")],
                &[("Consequence", "missense_variant"), ("SYMBOL", "TP53")],
            ],
        );
        let ctx = PickContext {
            index: &index,
            genes: None,
            mane: None,
        };
        let picked = pick(&rows, &[1, 2], &ctx).unwrap().unwrap();
        assert_eq!(picked.entry, 2);
    }

    #[test]
    fn read_genes_skips_blanks_and_comments() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genes.txt");
        std::fs::write(&path, "# panel\nTP53\n\n  BRCA1  \nENSG00000139618\n").unwrap();
        let panel = read_genes(&path).unwrap().value;
        assert_eq!(panel.0.len(), 3);
        assert!(panel.contains("BRCA1"));
        assert!(panel.contains("ENSG00000139618"));
        assert_eq!(panel.sorted(), ["BRCA1", "ENSG00000139618", "TP53"]);
    }

    #[test]
    fn read_genes_hashes_the_bytes_it_parsed() {
        const SHA256_OF_ABC: &str =
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genes.txt");
        std::fs::write(&path, "abc").unwrap();
        let panel = read_genes(&path).unwrap();
        assert_eq!(panel.sha256, SHA256_OF_ABC);
        assert!(panel.value.contains("abc"));
    }

    #[test]
    fn read_mane_table_hashes_the_bytes_it_read() {
        let path = Path::new("tests/data/mane-summary.txt");
        let expected = sha256_hex(&std::fs::read(path).unwrap());
        assert_eq!(read_mane_table(path).unwrap().sha256, expected);
    }

    fn assert_fixture_table(mane: &ManeTable) {
        assert_eq!(mane.select.len(), 4);
        assert!(mane.select.contains("ENST00000269305"));
        assert_eq!(
            mane.plus_clinical,
            table(&[], &["ENST00000591111"]).plus_clinical
        );
    }

    #[test]
    fn read_mane_table_reads_the_plain_fixture() {
        let mane = read_mane_table(Path::new("tests/data/mane-summary.txt"))
            .unwrap()
            .value;
        assert_fixture_table(&mane);
    }

    #[test]
    fn read_mane_table_reads_gzip_by_magic_bytes() {
        use std::io::Write;
        let text = std::fs::read("tests/data/mane-summary.txt").unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mane.summary.txt");
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(&text).unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();
        assert_fixture_table(&read_mane_table(&path).unwrap().value);
    }

    #[test]
    fn read_mane_table_names_a_missing_column() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mane.txt");
        std::fs::write(
            &path,
            "#NCBI_GeneID\tsymbol\tMANE_status\n1\tA\tMANE Select\n",
        )
        .unwrap();
        let error = read_mane_table(&path).unwrap_err();
        assert!(error.to_string().contains("Ensembl_nuc"), "{error}");
    }
}
