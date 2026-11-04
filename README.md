# annatar

[![Build Status](https://github.com/clintval/annatar/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/annatar/actions/workflows/ci.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Convert a VEP-annotated multi-sample VCF of single cells into an AnnData h5ad.

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    --channel conda-forge \
    --channel bioconda \
    annatar --help
```

## Introduction

The tool `annatar` reads a VEP-annotated multi-sample VCF or BCF, with cells as samples, and writes one h5ad with cells by variants.
Each ALT allele becomes a variant, numeric FORMAT fields become sparse layers, and one VEP consequence per variant becomes typed columns in `var`.
The h5ad can sit next to other single-cell modalities in a `MuData` keyed on cell barcode.

## Quick Start

Convert a VEP-annotated VCF of cells:

```bash
annatar --output 'cells.h5ad' 'cells.vep.vcf.gz'
```

Prefer the transcripts of a panel of genes, one HGNC symbol or Ensembl gene ID per line, and leave out CSQ subfields that are not needed:

```bash
annatar \
    --genes 'panel-genes.txt' \
    --exclude-field 'CSQ/Existing_variation' \
    --exclude-field 'CSQ/SOMATIC' \
    --output 'cells.h5ad' \
    'cells.vep.vcf.gz'
```

A field passed to `--exclude-field` must be declared in the VCF header, or annatar exits with an error.
The fields `FORMAT/AD`, `FORMAT/GT` and `INFO/CSQ` cannot be excluded.

Read from standard input:

```bash
bcftools view 'cells.vep.bcf' \
| annatar --output 'cells.h5ad'
```

## Output

`X` holds a genotype code per cell and variant:

| Value | Meaning |
| --- | --- |
| not stored | no genotype |
| 1 | no copies of this ALT |
| 2 | one copy of this ALT |
| 3 | two or more copies of this ALT; a haploid call counts as homozygous |

Layers are named from FORMAT fields, such as `AD_alt`. A stored `0` is a measured zero. A cell is covered when `GT` is called or `AD` sums above zero, so a caller that writes `0/0` for cells with no reads makes those cells covered.

A value that is missing at a covered cell is `NaN` in a float layer and unstored in an integer layer, so `AD_alt` can hold fewer entries than `X`. Align layers by index, not by position in `.data`.

The `var` table has one row per ALT, named `chrom:pos:ref:alt`, with the picked VEP consequence's subfields as columns. A subfield with no value in a file is an empty categorical, and types come from each file's data, so a column can differ between files.

## MuData Example

Color the cells that carry a pathogenic TP53 variant, next to a chromatin modality that shares the cell barcodes:

```python
import anndata as ad
import mudata as md

snv = ad.read_h5ad('cells.h5ad')
atac = ad.read_h5ad('cells.atac.h5ad')
mdata = md.MuData({'atac': atac, 'snv': snv})

var = mdata['snv'].var
hits = var.index[
    (var['SYMBOL'] == 'TP53')
    & var['ClinVar_CLNSIG'].isin({'Pathogenic', 'Likely_pathogenic'})
]
carriers = mdata['snv'][:, hits].layers['AD_alt'].sum(axis=1) > 0
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
