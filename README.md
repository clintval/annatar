![Banner for annatar](https://raw.githubusercontent.com/clintval/annatar/main/.github/img/annatar_banner.svg)

---

[![Build Status](https://github.com/clintval/annatar/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/annatar/actions/workflows/ci.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Convert a multi-sample VCF of single cells into an AnnData h5ad.

Install from GitHub with Rust and `cmake`, which builds the bundled HDF5:

```bash
cargo install --locked --git https://github.com/clintval/annatar
```

## Introduction

The tool `annatar` reads a multi-sample VCF or BCF, with cells as samples, and writes one h5ad with cells by variants.
Each ALT allele becomes a variant, and numeric FORMAT fields become sparse layers.
VEP annotation is optional. When the VCF has `INFO/CSQ`, one consequence per variant becomes typed columns in `var`.
The h5ad can sit next to other single-cell modalities in a `MuData` keyed on cell barcode.

## Quick Start

Convert a VCF of cells:

```bash
annatar --output 'cells.h5ad' 'cells.vcf.gz'
```

For a VEP-annotated VCF, prefer the transcripts of a panel of genes, one HGNC symbol or Ensembl gene ID per line, and leave out CSQ subfields that are not needed:

```bash
annatar \
    --genes 'panel-genes.txt' \
    --exclude-field 'CSQ/Existing_variation' \
    --exclude-field 'CSQ/SOMATIC' \
    --output 'cells.h5ad' \
    'cells.vep.vcf.gz'
```

A field passed to `--exclude-field` must be declared in the VCF header, or annatar exits with an error.
A `CSQ/` field only warns when the header has no `INFO/CSQ`, so one command serves annotated and plain VCFs.
Excluding `INFO/CSQ` drops all annotation. The fields `FORMAT/AD` and `FORMAT/GT` cannot be excluded.

Read from standard input:

```bash
bcftools view 'cells.bcf' \
| annatar --output 'cells.h5ad'
```

## Output

`X` holds a genotype code per cell and variant:

| Value | Meaning |
| --- | --- |
| not stored | no genotype, or a partial call such as `./1` |
| 1 | no copies of this ALT |
| 2 | one copy of this ALT |
| 3 | two or more copies of this ALT; a haploid call counts as homozygous |

Layers are named from FORMAT fields, such as `AD_alt`. A stored `0` is a measured zero. A cell is covered when `GT` is called or `AD` sums above zero, so a caller that writes `0/0` for cells with no reads makes those cells covered.

A value that is missing at a covered cell is `NaN` in a float layer and unstored in an integer layer, so `AD_alt` can hold fewer entries than `X`. Align layers by index, not by position in `.data`.

The `var` table has one row per ALT, named `chrom:pos:ref:alt`. The alleles `*`, `<*>` and `<NON_REF>` get no row.

VEP annotation is optional. A VCF with no `INFO/CSQ` gives no annotation columns. Otherwise the picked consequence's subfields become columns, with `consequence_rank`, `is_coding`, `is_lof` and `pick_tier`. A record with no CSQ keeps an empty row.

Column types come from subfield names and are pinned to VEP 116, so a subfield has one type in every file. Other VEP versions run with a warning. A value that fails its type warns and makes that column text in that file. A Consequence term VEP 116 does not know warns and ranks last, as 42 in `consequence_rank`.

VEP's `--allele_number` makes matching CSQ entries exact for symbolic ALTs, breakends and `--minimal` output, and multiallelic `--minimal` output needs it. Sequence alleles match on VEP's spelling regardless, so a sequence ALT never takes another allele's entry. A sequence ALT with no matching entry is an error, which annotating after splitting fixes. A symbolic ALT or breakend with no matching entry keeps an empty row and a warning.

## Combining Samples

On files with the same variants, `anndata.concat` with `merge='same'` keeps the annotation. On different variants it drops most columns, and `merge='first'` leaves the later files' variants blank, so rebuild `var` from the inputs:

```python
import anndata as ad
import pandas as pd

samples = {'a': ad.read_h5ad('a.h5ad'), 'b': ad.read_h5ad('b.h5ad')}
adata = ad.concat(samples, join='outer', label='sample', index_unique='-')
var = pd.concat([sample.var for sample in samples.values()])
adata.var = var[~var.index.duplicated()].reindex(adata.var_names)
```

Text columns can then come back as strings rather than categoricals, and read counts in `var`, such as `alt_reads`, come from the first file with that variant. An outer join on different variants also drops the stored zeros in layers, so an unstored `AD_alt` no longer means an uncovered cell.

## MuData Example

Color the cells that carry a pathogenic TP53 variant, next to a chromatin modality that shares the cell barcodes:

```python
import anndata as ad
import mudata as md
import numpy as np

snv = ad.read_h5ad('cells.h5ad')
atac = ad.read_h5ad('cells.atac.h5ad')
mdata = md.MuData({'atac': atac, 'snv': snv})

var = mdata['snv'].var
pathogenic = {'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}
hits = var.index[(var['SYMBOL'] == 'TP53') & var['ClinVar_CLNSIG'].isin(pathogenic)]
carriers = np.asarray(mdata['snv'][:, hits].layers['AD_alt'].sum(axis=1)).ravel() > 0
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
