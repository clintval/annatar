import gzip
import os
import re
import shlex
import subprocess
import sys
import warnings
from importlib.util import find_spec
from unittest.mock import MagicMock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import BINARY
from conftest import REPO
from conftest import logged

import annatar

EXAMPLE = REPO / "docs" / "example"
VCF = EXAMPLE / "cells.vcf.gz"
ATAC = EXAMPLE / "cells.atac.h5ad"
TP53_PATHOGENIC = ["chr17:7673802:C:A", "chr17:7673802:C:T"]
TP53_BENIGN = "chr17:7676154:G:C"
TP53_STATE_BY_CELL_TYPE = {
    "T cell": {"carrier": 0, "reference only": 181, "no reads": 439},
    "Monocyte": {"carrier": 0, "reference only": 129, "no reads": 351},
    "Tumor": {"carrier": 89, "reference only": 8, "no reads": 283},
    "B cell": {"carrier": 0, "reference only": 75, "no reads": 225},
    "NK cell": {"carrier": 0, "reference only": 55, "no reads": 165},
}

pytestmark = [
    pytest.mark.skipif(not EXAMPLE.is_dir(), reason="the sdist leaves out docs/example"),
    pytest.mark.filterwarnings("ignore:The inplace keyword in DataFrame"),
]


def readme_block(language):
    readme = (REPO / "README.md").read_text()
    section = readme.split("## MuData Example\n")[1].split("\n## ")[0]
    blocks = re.findall(rf"```{language}\n(.*?)```", section, flags=re.DOTALL)
    assert len(blocks) == 1
    return blocks[0]


def run_annatar(*arguments):
    environment = {**os.environ, "NO_COLOR": "1"}
    environment.pop("RUST_LOG", None)
    command = [BINARY, *arguments]
    return subprocess.run(command, capture_output=True, text=True, env=environment, check=True)


def carriers_of(variants, names):
    alt = np.asarray(variants[:, names].layers["AD_alt"].sum(axis=1)).ravel()
    return variants.obs_names[alt > 0]


@pytest.fixture
def example(tmp_path, monkeypatch):
    (tmp_path / "docs").symlink_to(EXAMPLE.parent)
    monkeypatch.chdir(tmp_path)
    program, *arguments = shlex.split(readme_block("bash"))
    assert program == "annatar"
    run_annatar(*arguments)
    namespace = {}
    if find_spec("muon") is None:
        monkeypatch.setitem(sys.modules, "muon", MagicMock())
        exec(readme_block("python"), namespace)
        return namespace

    import matplotlib

    matplotlib.use("Agg")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(readme_block("python"), namespace)
    return namespace


def test_the_readme_example_finds_carriers_only_in_the_tumor_cluster(example):
    obs = example["mu"].obs
    states = obs["variants:TP53"].astype(object).fillna("no reads")
    table = pd.crosstab(obs["atac:cell_type"], states)
    assert table.to_dict("index") == TP53_STATE_BY_CELL_TYPE


def test_the_readme_example_selects_both_pathogenic_tp53_alleles(example):
    variants = example["variants"]
    assert list(variants.var_names[example["hits"]]) == TP53_PATHOGENIC
    assert TP53_BENIGN in variants.var_names[variants.var["SYMBOL"] == "TP53"]


def test_the_readme_plot_call_names_what_the_mudata_holds(example):
    muon = example["muon"]
    if not isinstance(muon, MagicMock):
        pytest.skip("the real muon checks this by plotting")
    call = muon.pl.embedding.call_args
    mu = example["mu"]
    assert call.args == (mu,)
    assert set(call.kwargs["color"]) <= set(mu.obs.columns)
    modality, basis = call.kwargs["basis"].split(":")
    assert basis in mu[modality].obsm


def test_carriers_agree_by_barcode_without_mudata():
    atac = ad.read_h5ad(ATAC)
    variants = annatar.read_vcf(VCF)
    assert set(atac.obs_names) == set(variants.obs_names)
    assert list(atac.obs_names) != list(variants.obs_names)
    cell_types = atac.obs.loc[carriers_of(variants, TP53_PATHOGENIC), "cell_type"]
    expected = {kind: row["carrier"] for kind, row in TP53_STATE_BY_CELL_TYPE.items()}
    assert cell_types.value_counts().to_dict() == expected


def test_a_germline_variant_is_seen_across_every_cell_type():
    atac = ad.read_h5ad(ATAC)
    variants = annatar.read_vcf(VCF)
    cell_types = atac.obs.loc[carriers_of(variants, [TP53_BENIGN]), "cell_type"]
    assert set(cell_types) == set(TP53_STATE_BY_CELL_TYPE)


def test_annatar_writes_no_warnings_on_the_vcf(tmp_path):
    run = run_annatar(VCF, tmp_path / "cells.h5ad")
    assert logged(run, "INFO")
    assert logged(run, "WARN") == []
    assert "WARN" not in run.stderr

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        annatar.read_vcf(VCF)
    assert [str(warning.message) for warning in caught] == []


def test_the_vcf_is_split_sparse_and_marked_simulated():
    with gzip.open(VCF, "rt") as lines:
        text = lines.read()
    records = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
    assert '##VEP="v116.0"' in text
    assert "##simulated=" in text
    assert not any("," in fields[4] for fields in records)
    assert [fields[4] for fields in records if fields[1] == "7673802"] == ["A", "T"]

    variants = annatar.read_vcf(VCF)
    assert (variants.var["n_cells_covered"] < variants.n_obs / 2).all()
    assert (variants.var["n_cells_alt"] > 0).all()
    assert (variants.var.loc[TP53_PATHOGENIC, "pick_tier"] == "mane_select").all()
    assert variants.var.loc["chr17:7673802:C:T", "HGVSp"].endswith("p.Arg273His")


def test_the_chromatin_file_is_marked_simulated_and_has_an_embedding():
    atac = ad.read_h5ad(ATAC)
    assert atac.uns["simulated"]
    assert atac.obsm["X_umap"].shape == (atac.n_obs, 2)
    assert set(atac.obs["cell_type"]) == set(TP53_STATE_BY_CELL_TYPE)
    assert atac.X.nnz > 0


def test_the_example_files_stay_small():
    assert 1_500 <= ad.read_h5ad(ATAC).n_obs <= 3_000
    assert VCF.stat().st_size < 1_000_000
    assert ATAC.stat().st_size < 1_000_000
