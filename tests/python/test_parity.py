import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import DATA
from conftest import read_quietly

FIXTURES = sorted(
    path for path in DATA.iterdir() if path.name.endswith((".vcf", ".vcf.gz", ".bcf"))
)
REFUSED_BY_THE_CLI = {"sites-only.vcf"}
CONVERTED = [path for path in FIXTURES if path.name not in REFUSED_BY_THE_CLI]
PANEL = DATA / "panel-genes.txt"
MANE = DATA / "mane-summary.txt"


def assert_same_matrix(left, right):
    assert type(left) is type(right)
    assert left.shape == right.shape
    assert left.dtype == right.dtype
    assert left.indptr.dtype == right.indptr.dtype
    assert left.indices.dtype == right.indices.dtype
    np.testing.assert_array_equal(left.indptr, right.indptr)
    np.testing.assert_array_equal(left.indices, right.indices)
    np.testing.assert_array_equal(left.data, right.data)


def assert_same_frame(left, right):
    pd.testing.assert_frame_equal(
        left, right, check_exact=True, check_column_type=True, check_index_type=False
    )
    assert left.index.name == right.index.name
    assert left.columns.name == right.columns.name
    for name in left.columns:
        assert left[name].dtype == right[name].dtype, name
        if isinstance(left[name].dtype, pd.CategoricalDtype):
            assert left[name].cat.categories.dtype == right[name].cat.categories.dtype, name
            assert list(left[name].cat.categories) == list(right[name].cat.categories), name
            assert left[name].cat.ordered == right[name].cat.ordered, name


def assert_same_uns(left, right):
    assert list(left) == list(right)
    for key in left.keys() - {"command_line"}:
        assert type(left[key]) is type(right[key]), key
        if isinstance(left[key], np.ndarray):
            assert left[key].dtype == right[key].dtype, key
            np.testing.assert_array_equal(left[key], right[key])
        else:
            assert left[key] == right[key], key


def assert_same_anndata(left, right):
    assert left.shape == right.shape
    assert_same_matrix(left.X, right.X)
    assert list(left.layers) == list(right.layers)
    for name in left.layers:
        assert_same_matrix(left.layers[name], right.layers[name])
    assert_same_frame(left.obs, right.obs)
    assert_same_frame(left.var, right.var)
    assert_same_uns(left.uns, right.uns)
    for slot in ("obsm", "varm", "obsp", "varp"):
        assert list(getattr(left, slot)) == list(getattr(right, slot)) == []


@pytest.mark.parametrize("path", CONVERTED, ids=lambda path: path.name)
def test_read_vcf_equals_the_cli_output(path, run_cli, cli_output):
    run = run_cli(path)
    assert run.returncode == 0, run.stderr
    assert_same_anndata(read_quietly(path), ad.read_h5ad(cli_output))


@pytest.mark.parametrize(
    ("name", "flags", "options"),
    [
        ("basic.vcf", ["--genes", PANEL], {"genes": PANEL}),
        ("no-csq.vcf", ["--genes", PANEL], {"genes": PANEL}),
        ("basic.vcf", ["--mane-table", MANE], {"mane_table": MANE}),
        ("basic.vcf", ["--exclude-field", "INFO/CSQ"], {"exclude_fields": ["INFO/CSQ"]}),
        (
            "basic.vcf",
            ["--exclude-field", "INFO/CSQ/SYMBOL"],
            {"exclude_fields": ("INFO/CSQ/SYMBOL",)},
        ),
        (
            "typed-fields.vcf",
            ["--exclude-field", "INFO/AC", "--exclude-field", "FORMAT/DP"],
            {"exclude_fields": ["INFO/AC", "FORMAT/DP"]},
        ),
        ("basic.vcf.gz", ["--threads", "2"], {"threads": 2}),
    ],
)
def test_options_equal_their_cli_flags(name, flags, options, run_cli, cli_output):
    run = run_cli(DATA / name, *flags)
    assert run.returncode == 0, run.stderr
    assert_same_anndata(read_quietly(DATA / name, **options), ad.read_h5ad(cli_output))


@pytest.mark.parametrize("name", ["basic.vcf", "no-records.vcf", "no-csq.vcf"])
def test_a_written_result_reads_back_the_same(name, tmp_path):
    adata = read_quietly(DATA / name)
    adata.write_h5ad(tmp_path / "written.h5ad")
    assert_same_anndata(adata, ad.read_h5ad(tmp_path / "written.h5ad"))
