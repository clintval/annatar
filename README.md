![Banner for annatar](.github/img/annatar_banner.svg)

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
annatar 'cells.vcf.gz' 'cells.h5ad'
```

For a VEP-annotated VCF, prefer panel genes and leave out CSQ subfields:

```bash
annatar \
    --genes 'panel-genes.txt' \
    --exclude-field 'INFO/CSQ/Existing_variation' \
    --exclude-field 'INFO/CSQ/SOMATIC' \
    'cells.vep.vcf.gz' \
    'cells.h5ad'
```

Excluding `INFO/CSQ` drops all annotation. The fields `FORMAT/AD` and `FORMAT/GT` cannot be excluded.

## Output

| Slot | Contents |
| --- | --- |
| `X` | Copies of each ALT in each cell's called genotype |
| `layers` | Numeric FORMAT fields, such as `AD_ref`, `AD_alt` and `DP` |
| `obs` | Per-cell counts, `n_sites_covered` and `n_sites_alt` |
| `var` | One row per ALT allele, named `chrom:pos:ref:alt` |
| `uns` | Provenance and the encoding of `X` |

#### Genotypes

`X` holds the number of copies of the ALT in the called genotype, so `0/1` is 1 and `1/1` is 2, and hom-ref and uncalled cells are not stored.
So `X > 0` means a called carrier and `X.sum()` counts ALT haplotypes.
A haploid ALT call is 1, so `X == 2` does not find it.

Layers are read counts, with zeros and missing values not stored.
At a variant a cell is a carrier (`X > 0`, or `AD_alt > 0` for reads), reference (`AD_ref > 0` with no `X` or `AD_alt`), or has no data.
A cell with reads but no call keeps them in the layers, and `AD_alt` counts alt reads in cells called `0/0`.

#### Variants

Every ALT gets its own row, whether it is an SNV, indel, MNV, symbolic ALT or breakend.
The placeholder alleles `*`, `<*>` and `<NON_REF>` get no row.
The VCF must declare `FORMAT/AD`.

#### VEP Annotation

VEP is optional.
When the VCF has `INFO/CSQ`, one consequence per variant is added to `var` with its VEP subfield names, along with these columns:

| Column | Meaning |
| --- | --- |
| `consequence_rank` | Severity of the consequence, where lower is more severe |
| `is_coding` | The consequence is coding |
| `is_lof` | The consequence is loss of function |
| `pick_tier` | `mane_select`, `mane_plus_clinical`, `canonical` or `other` |
| `worst_panel_impact` | Highest `IMPACT` among the `--genes` panel genes |

Column types are fixed for VEP 116, so a column has the same type in every file.
Other VEP versions run with a warning.

> [!IMPORTANT]
> Split multiallelic records before running VEP, not after.
> annatar stops on evidence of a split after VEP, including a `bcftools norm` split line in the header, even if no record changed.
> It cannot detect a split after `--pick` when neither `--allele_number` nor that header line is there.
> annatar refuses `--minimal` output, so run VEP without it.
> Pass `--allele_number` to VEP so symbolic ALTs and breakends match by position.

## Merging Libraries

Merge the h5ad files of several libraries into one AnnData.
`anndata.concat` loses `var` annotation when the files hold different variants, so rebuild `var` from the inputs:

```python
import anndata as ad
import pandas as pd

libraries = {'a': ad.read_h5ad('a.h5ad'), 'b': ad.read_h5ad('b.h5ad')}
adata = ad.concat(libraries, join='outer', label='library', index_unique='-')
var = pd.concat([library.var for library in libraries.values()])
adata.var = var[~var.index.duplicated()].reindex(adata.var_names)
```

After merging, read counts in `var`, such as `alt_reads`, describe the first file with that variant.
`X` and the layers keep every stored value, since neither stores a zero.

## MuData Example

Color the cells that carry a pathogenic TP53 variant on a chromatin UMAP.
`MuData` lines the modalities up by cell barcode, and cells with no reads at these variants stay gray:

```python
import anndata as ad
import mudata as md
import muon
import numpy as np

mu = md.MuData({'atac': ad.read_h5ad('cells.atac.h5ad'), 'variants': ad.read_h5ad('cells.h5ad')})

variants = mu['variants']
pathogenic = {'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}
hits = (variants.var['SYMBOL'] == 'TP53') & variants.var['ClinVar_CLNSIG'].isin(pathogenic)
alt = np.asarray(variants[:, hits].layers['AD_alt'].sum(axis=1)).ravel()
ref = np.asarray(variants[:, hits].layers['AD_ref'].sum(axis=1)).ravel()
variants.obs['TP53'] = np.select([alt > 0, ref > 0], ['carrier', 'reference only'], default=None)

mu.pull_obs('TP53', mods='variants')
muon.pl.embedding(mu, basis='atac:X_umap', color='variants:TP53')
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
