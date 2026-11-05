# Requires scanpy and scikit-learn, which the pixi environment does not provide.
"""Simulate the paired VCF and chromatin h5ad of the README's MuData example.

Every barcode, read and tile is invented. Run from the repository root.
"""

import subprocess
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import TruncatedSVD

OUTPUT = Path(__file__).resolve().parent
SEED = 273

CELL_TYPES = {"T cell": 620, "Monocyte": 480, "Tumor": 380, "B cell": 300, "NK cell": 220}
LINEAGE = {
    "T cell": "lymphoid",
    "NK cell": "lymphoid",
    "B cell": "lymphoid",
    "Monocyte": "myeloid",
    "Tumor": "epithelial",
}
EVERY_CELL = dict.fromkeys(CELL_TYPES, 1.0)

TILES = 1600
TILES_PER_PROGRAM = 140
FRAGMENTS_PER_CELL = 170
EXTRA_READS_PER_COVERED_SITE = 0.8

CONTIGS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

CSQ_FIELDS = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp|"
    "Protein_position|ALLELE_NUM|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|CCDS|SWISSPROT|"
    "ClinVar_CLNSIG|ClinVar_CLNDN"
).split("|")

TP53 = {"SYMBOL": "TP53", "Gene": "ENSG00000141510", "SWISSPROT": "P04637"}
TP53_MANE = {
    **TP53,
    "Feature": "ENST00000269305.9",
    "CANONICAL": "YES",
    "MANE": "MANE_Select",
    "MANE_SELECT": "NM_000546.6",
    "CCDS": "CCDS11118.1",
}
TP53_OTHER = {**TP53, "Feature": "ENST00000445888.6"}
DNMT3A = {
    "SYMBOL": "DNMT3A",
    "Gene": "ENSG00000119772",
    "Feature": "ENST00000264709.7",
    "CANONICAL": "YES",
    "MANE": "MANE_Select",
    "MANE_SELECT": "NM_022552.5",
    "CCDS": "CCDS33157.1",
    "SWISSPROT": "Q9Y6K1",
}
KRAS = {"SYMBOL": "KRAS", "Gene": "ENSG00000133703", "SWISSPROT": "P01116"}
BRCA2 = {
    "SYMBOL": "BRCA2",
    "Gene": "ENSG00000139618",
    "Feature": "ENST00000380152.8",
    "CANONICAL": "YES",
    "MANE": "MANE_Select",
    "MANE_SELECT": "NM_000059.4",
    "CCDS": "CCDS9344.1",
    "SWISSPROT": "P51587",
}
ASXL1 = {
    "SYMBOL": "ASXL1",
    "Gene": "ENSG00000171456",
    "Feature": "ENST00000375687.10",
    "CANONICAL": "YES",
    "MANE": "MANE_Select",
    "MANE_SELECT": "NM_015338.6",
    "CCDS": "CCDS13201.1",
    "SWISSPROT": "Q8IXJ9",
}


def missense(transcript: dict, cdna: str, protein_id: str, protein: str, position: int, **more):
    """Return the CSQ values of a missense change on a transcript."""
    return {
        **transcript,
        "Consequence": "missense_variant",
        "IMPACT": "MODERATE",
        "HGVSc": f"{transcript['Feature']}:{cdna}",
        "HGVSp": f"{protein_id}:{protein}",
        "Protein_position": str(position),
        **more,
    }


@dataclass(frozen=True)
class Allele:
    """One ALT, the share of each cell type that carries it, and its CSQ entries."""

    alt: str
    entries: list[dict]
    carriers: dict[str, float]
    vaf: float = 0.5
    vep_allele: str = ""


@dataclass(frozen=True)
class Site:
    """One position, its expected reads in a cell of average depth, and its ALTs."""

    chrom: str
    pos: int
    ref: str
    coverage: float
    alleles: list[Allele] = field(default_factory=list)


SITES = [
    Site(
        "chr2",
        25234373,
        "C",
        0.16,
        [
            Allele(
                "T",
                [
                    missense(
                        DNMT3A,
                        "c.2645G>A",
                        "ENSP00000264709.3",
                        "p.Arg882His",
                        882,
                        ClinVar_CLNSIG="Pathogenic",
                        ClinVar_CLNDN="Acute_myeloid_leukemia",
                    )
                ],
                {"Monocyte": 0.35},
            )
        ],
    ),
    Site(
        "chr2",
        25247044,
        "G",
        0.12,
        [
            Allele(
                "A",
                [
                    {
                        **DNMT3A,
                        "Consequence": "intron_variant",
                        "IMPACT": "MODIFIER",
                        "HGVSc": "ENST00000264709.7:c.1851+38C>T",
                    }
                ],
                EVERY_CELL,
                vaf=1.0,
            )
        ],
    ),
    Site(
        "chr12",
        25245350,
        "C",
        0.14,
        [
            Allele(
                "T",
                [
                    missense(
                        {
                            **KRAS,
                            "Feature": "ENST00000311936.8",
                            "CANONICAL": "YES",
                            "MANE": "MANE_Select",
                            "MANE_SELECT": "NM_004985.5",
                            "CCDS": "CCDS8703.1",
                        },
                        "c.35G>A",
                        "ENSP00000308495.3",
                        "p.Gly12Asp",
                        12,
                        ClinVar_CLNSIG="Pathogenic",
                        ClinVar_CLNDN="Neoplasm",
                    ),
                    missense(
                        {
                            **KRAS,
                            "Feature": "ENST00000256078.10",
                            "MANE": "MANE_Plus_Clinical",
                            "MANE_PLUS_CLINICAL": "NM_033360.4",
                            "CCDS": "CCDS8702.1",
                        },
                        "c.35G>A",
                        "ENSP00000256078.5",
                        "p.Gly12Asp",
                        12,
                        ClinVar_CLNSIG="Pathogenic",
                        ClinVar_CLNDN="Neoplasm",
                    ),
                ],
                {"Tumor": 0.3},
            )
        ],
    ),
    Site(
        "chr13",
        32332592,
        "A",
        0.18,
        [
            Allele(
                "C",
                [
                    missense(
                        BRCA2,
                        "c.1114A>C",
                        "ENSP00000369497.3",
                        "p.Asn372His",
                        372,
                        ClinVar_CLNSIG="Benign",
                        ClinVar_CLNDN="not_specified",
                    )
                ],
                EVERY_CELL,
            )
        ],
    ),
    Site(
        "chr17",
        7673802,
        "C",
        0.3,
        [
            Allele(
                "A",
                [
                    missense(
                        TP53_MANE,
                        "c.818G>T",
                        "ENSP00000269305.4",
                        "p.Arg273Leu",
                        273,
                        ClinVar_CLNSIG="Likely_pathogenic",
                        ClinVar_CLNDN="Li-Fraumeni_syndrome",
                    ),
                    missense(TP53_OTHER, "c.818G>T", "ENSP00000391478.2", "p.Arg273Leu", 273),
                ],
                {"Tumor": 0.15},
                vaf=0.85,
            ),
            Allele(
                "T",
                [
                    missense(
                        TP53_MANE,
                        "c.818G>A",
                        "ENSP00000269305.4",
                        "p.Arg273His",
                        273,
                        ClinVar_CLNSIG="Pathogenic",
                        ClinVar_CLNDN="Li-Fraumeni_syndrome",
                    ),
                    missense(TP53_OTHER, "c.818G>A", "ENSP00000391478.2", "p.Arg273His", 273),
                ],
                {"Tumor": 0.85},
                vaf=0.85,
            ),
        ],
    ),
    Site(
        "chr17",
        7676154,
        "G",
        0.2,
        [
            Allele(
                "C",
                [
                    missense(
                        TP53_MANE,
                        "c.215C>G",
                        "ENSP00000269305.4",
                        "p.Pro72Arg",
                        72,
                        ClinVar_CLNSIG="Benign",
                        ClinVar_CLNDN="not_specified",
                    )
                ],
                EVERY_CELL,
            )
        ],
    ),
    Site(
        "chr20",
        32434638,
        "A",
        0.16,
        [
            Allele(
                "AG",
                [
                    {
                        **ASXL1,
                        "Consequence": "frameshift_variant",
                        "IMPACT": "HIGH",
                        "HGVSc": "ENST00000375687.10:c.1934dup",
                        "HGVSp": "ENSP00000364839.4:p.Gly646TrpfsTer12",
                        "Protein_position": "645-646",
                        "ClinVar_CLNSIG": "Pathogenic",
                        "ClinVar_CLNDN": "C-like_syndrome",
                    }
                ],
                {"Monocyte": 0.2},
                vep_allele="G",
            )
        ],
    ),
]


def header() -> list[str]:
    """Return the VCF meta-information lines, with the split recorded before VEP."""
    command = (
        "vep --allele_number --cache --canonical --ccds --format vcf --hgvs --mane --offline "
        "--symbol --uniprot --vcf --custom file=clinvar.vcf.gz,short_name=ClinVar,format=vcf,"
        "fields=CLNSIG%CLNDN"
    )
    return [
        "##fileformat=VCFv4.3",
        "##reference=GRCh38",
        '##simulated="Every cell, read and annotation in this file is invented"',
        *(f"##contig=<ID={name},length={length}>" for name, length in CONTIGS.items()),
        '##FILTER=<ID=PASS,Description="All filters passed">',
        "##bcftools_normVersion=1.24+htslib-1.24",
        "##bcftools_normCommand=norm --multiallelics -any cells.vcf; Date=Wed Nov  4 12:00:00 2026",
        '##VEP="v116.0" time="2026-11-04 12:10:00" cache="homo_sapiens/116_GRCh38" ensembl=116.0 '
        'assembly="GRCh38"',
        f"##VEP-command-line='{command}'",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl '
        f'VEP. Format: {"|".join(CSQ_FIELDS)}">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    ]


def csq(allele: Allele) -> str:
    """Return the INFO/CSQ value of one ALT."""
    common = {
        "Allele": allele.vep_allele or allele.alt,
        "Feature_type": "Transcript",
        "BIOTYPE": "protein_coding",
        "ALLELE_NUM": "1",
    }
    return ",".join(
        "|".join({**common, **entry}.get(name, "") for name in CSQ_FIELDS)
        for entry in allele.entries
    )


def call(ref: int, alt: int) -> str:
    """Return the GT:AD of a cell, with no call when no read supports this record."""
    if ref + alt == 0:
        return "./.:."
    genotype = "0/0" if alt == 0 else "1/1" if ref == 0 and alt >= 3 else "0/1"
    return f"{genotype}:{ref},{alt}"


def carried_alleles(rng: np.random.Generator, site: Site, labels: np.ndarray) -> np.ndarray:
    """Return the index of the one ALT each cell carries at a site, or -1 for none."""
    draw = rng.random(len(labels))
    carried = np.full(len(labels), -1)
    floor = np.zeros(len(labels))
    for index, allele in enumerate(site.alleles):
        ceiling = floor + np.array([allele.carriers.get(kind, 0.0) for kind in labels])
        carried[(draw >= floor) & (draw < ceiling)] = index
        floor = ceiling
    return carried


def records(rng: np.random.Generator, site: Site, labels: np.ndarray, depth: np.ndarray):
    """Return the VCF lines of a site, where split ALTs share one pileup of reads."""
    carried = carried_alleles(rng, site, labels)
    covered = rng.random(len(labels)) < 1 - np.exp(-site.coverage * depth)
    reads = covered * (1 + rng.poisson(EXTRA_READS_PER_COVERED_SITE, len(labels)))
    vaf = np.array([site.alleles[index].vaf if index >= 0 else 0.0 for index in carried])
    alt_reads = rng.binomial(reads, vaf)
    ref_reads = reads - alt_reads
    for index, allele in enumerate(site.alleles):
        alt = np.where(carried == index, alt_reads, 0)
        calls = [call(*pair) for pair in zip(ref_reads.tolist(), alt.tolist(), strict=True)]
        qual = str(30 + 2 * int(alt.sum()))
        fixed = [site.chrom, str(site.pos), ".", site.ref, allele.alt, qual, "PASS"]
        yield "\t".join([*fixed, f"CSQ={csq(allele)}", "GT:AD", *calls])


def write_vcf(rng: np.random.Generator, barcodes, labels: np.ndarray, depth: np.ndarray) -> None:
    """Write the bgzipped per-cell VCF."""
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *barcodes]
    lines = [*header(), "\t".join(columns)]
    for site in SITES:
        lines.extend(records(rng, site, labels, depth))
    with (OUTPUT / "cells.vcf.gz").open("wb") as sink:
        text = "\n".join(lines) + "\n"
        subprocess.run(["bgzip"], input=text.encode(), stdout=sink, check=True)


def program(rng: np.random.Generator) -> np.ndarray:
    """Return tile weights that open a random set of tiles."""
    weights = np.zeros(TILES)
    opened = rng.choice(TILES, TILES_PER_PROGRAM, replace=False)
    weights[opened] = rng.gamma(2.0, 1.0, TILES_PER_PROGRAM)
    return weights


def accessibility(rng: np.random.Generator, labels: np.ndarray) -> np.ndarray:
    """Return each cell's share of fragments in each tile, from nested lineage programs."""
    base = rng.gamma(0.8, 1.0, TILES)
    lineage = {name: program(rng) for name in sorted(set(LINEAGE.values()))}
    shares = np.zeros((len(labels), TILES))
    for kind in CELL_TYPES:
        cells = labels == kind
        state = rng.beta(0.8, 0.8, cells.sum())[:, None]
        drift = state * program(rng) + (1 - state) * program(rng)
        shares[cells] = base + lineage[LINEAGE[kind]] + 0.4 * program(rng) + 0.8 * drift
    return shares / shares.sum(axis=1, keepdims=True)


def tile_names(rng: np.random.Generator) -> pd.Index:
    """Return sorted, unique, genome-wide 500 bp tile names."""
    edges = np.cumsum([0, *CONTIGS.values()])
    starts = np.sort(rng.choice(edges[-1] // 500, TILES, replace=False)) * 500
    contigs = np.searchsorted(edges, starts, side="right") - 1
    names = list(CONTIGS)
    return pd.Index(
        [
            f"{names[contig]}:{start - edges[contig]}-{start - edges[contig] + 500}"
            for contig, start in zip(contigs, starts, strict=True)
        ],
        dtype=object,
    )


def lsi(counts: sparse.csr_matrix) -> np.ndarray:
    """Return TF-IDF plus SVD components, without the first, which tracks depth."""
    term = sparse.diags(1 / counts.sum(axis=1).A1) @ counts
    inverse = sparse.diags(counts.shape[0] / np.maximum(counts.getnnz(axis=0), 1))
    weighted = (term @ inverse).tocsr()
    weighted.data = np.log1p(1e4 * weighted.data)
    return TruncatedSVD(n_components=30, random_state=SEED).fit_transform(weighted)[:, 1:]


def write_atac(rng: np.random.Generator, barcodes, labels: np.ndarray, depth: np.ndarray) -> None:
    """Simulate tile counts, embed them and write the chromatin h5ad."""
    shares = accessibility(rng, labels)
    counts = sparse.csr_matrix(rng.poisson(shares * FRAGMENTS_PER_CELL * depth[:, None]))
    kinds = pd.Index(list(CELL_TYPES), dtype=object)
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(labels, categories=kinds)},
        index=pd.Index(barcodes, dtype=object),
    )
    adata = ad.AnnData(X=counts.astype(np.int32), obs=obs, var=pd.DataFrame(index=tile_names(rng)))
    adata.obsm["X_lsi"] = lsi(counts)
    sc.pp.neighbors(adata, n_neighbors=40, use_rep="X_lsi", metric="cosine", random_state=SEED)
    sc.tl.umap(adata, min_dist=0.8, random_state=SEED)
    shipped = ad.AnnData(
        X=adata.X,
        obs=adata.obs,
        var=adata.var,
        obsm={"X_umap": adata.obsm["X_umap"].astype(np.float32)},
        uns={"simulated": True},
    )
    shipped[np.argsort(barcodes)].copy().write_h5ad(
        OUTPUT / "cells.atac.h5ad", compression="gzip", compression_opts=9
    )


def main() -> None:
    """Write both example files."""
    rng = np.random.default_rng(SEED)
    labels = rng.permutation(np.repeat(list(CELL_TYPES), list(CELL_TYPES.values())))
    letters = rng.choice(list("ACGT"), (2 * len(labels), 16))
    barcodes = list(dict.fromkeys("".join(row) for row in letters))[: len(labels)]
    depth = rng.lognormal(0.0, 0.4, len(labels))
    write_vcf(rng, barcodes, labels, depth)
    write_atac(rng, barcodes, labels, depth)


if __name__ == "__main__":
    main()
