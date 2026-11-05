import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from conftest import DATA
from conftest import read_quietly
from test_parity import assert_same_anndata

import annatar

COUNTS = ["n_cells_covered", "n_cells_alt", "ref_reads", "alt_reads"]
SHARED = "chr1:100:G:A"


def read(name, **options):
    return read_quietly(DATA / name, **options)


def written(tmp_path, name):
    path = tmp_path / f"{name}.h5ad"
    read(f"{name}.vcf").write_h5ad(path)
    return path


def test_different_variants_keep_their_annotation():
    merged = annatar.concat([read("csq-types-a.vcf"), read("csq-types-b.vcf")])
    assert sorted(merged.var_names) == [SHARED, "chr1:200:C:T", "chr1:300:CA:C"]
    assert merged.var.loc["chr1:200:C:T", "Consequence"] == "missense_variant"
    assert merged.var.loc["chr1:300:CA:C", "Consequence"] == "frameshift_variant"
    assert merged.var.loc[SHARED, "Consequence"] == "stop_gained"
    assert merged.var["SYMBOL"].notna().all()


def test_cells_are_labelled_with_their_library():
    merged = annatar.concat({"a": read("csq-types-a.vcf"), "b": read("csq-types-b.vcf")})
    assert list(merged.obs_names) == ["c1-a", "c2-a", "c1-b", "c2-b"]
    assert list(merged.obs["library"]) == ["a", "a", "b", "b"]


def test_shared_variants_get_summed_counts():
    merged = annatar.concat([read("csq-types-a.vcf"), read("csq-types-b.vcf")])
    assert merged.var.loc[SHARED, COUNTS].tolist() == [4, 2, 24, 6]
    assert merged.var.loc["chr1:200:C:T", COUNTS].tolist() == [2, 1, 12, 3]
    assert merged.var.loc["chr1:300:CA:C", COUNTS].tolist() == [2, 1, 12, 3]


def test_counts_come_from_the_merged_layers():
    a, b = read("basic.vcf"), read("multiallelic.vcf")
    merged = annatar.concat([a, b])
    alt = merged.layers["AD_alt"]
    ref = merged.layers["AD_ref"]
    assert merged.var["alt_reads"].tolist() == np.asarray(alt.sum(axis=0)).ravel().tolist()
    assert merged.var["ref_reads"].tolist() == np.asarray(ref.sum(axis=0)).ravel().tolist()
    covered = np.asarray(((ref + alt) > 0).sum(axis=0)).ravel()
    assert merged.var["n_cells_covered"].tolist() == covered.tolist()
    assert merged.var["n_cells_alt"].tolist() == np.asarray((alt > 0).sum(axis=0)).ravel().tolist()
    assert all(merged.var[name].dtype == np.int64 for name in COUNTS)


def test_paths_and_anndata_give_the_same_merge(tmp_path):
    a, b = written(tmp_path, "csq-types-a"), written(tmp_path, "csq-types-b")
    from_paths = annatar.concat({"a": a, "b": str(b)})
    from_objects = annatar.concat({"a": read("csq-types-a.vcf"), "b": read("csq-types-b.vcf")})
    pd.testing.assert_frame_equal(from_paths.var, from_objects.var, check_categorical=False)
    assert (from_paths.X != from_objects.X).nnz == 0


def test_a_list_of_paths_names_libraries_by_file(tmp_path):
    a, b = written(tmp_path, "csq-types-a"), written(tmp_path, "csq-types-b")
    merged = annatar.concat([a, b])
    assert merged.obs["library"].unique().tolist() == ["csq-types-a", "csq-types-b"]


def test_a_list_of_anndata_names_libraries_by_position():
    merged = annatar.concat([read("csq-types-a.vcf"), read("csq-types-b.vcf")])
    assert merged.obs["library"].unique().tolist() == ["0", "1"]


def test_files_with_one_name_must_be_named(tmp_path):
    (tmp_path / "x").mkdir()
    (tmp_path / "y").mkdir()
    read("no-csq.vcf").write_h5ad(tmp_path / "x" / "cells.h5ad")
    read("no-csq.vcf").write_h5ad(tmp_path / "y" / "cells.h5ad")
    with pytest.raises(ValueError, match="more than one library is named cells"):
        annatar.concat([tmp_path / "x" / "cells.h5ad", tmp_path / "y" / "cells.h5ad"])


@pytest.mark.parametrize(
    "names",
    [
        ("csq-types-a.vcf", "csq-types-fallback.vcf", "csq-types-b.vcf"),
        ("csq-types-fallback.vcf", "csq-types-b.vcf", "csq-types-a.vcf"),
    ],
)
def test_no_text_becomes_nan(names):
    libraries = [read(name) for name in names]
    merged = annatar.concat(libraries)
    held = set()
    for library in libraries:
        fresh = library.var.index.difference(held)
        held |= set(fresh)
        for name in library.var.columns:
            present = fresh[library.var.loc[fresh, name].notna()]
            assert merged.var.loc[present, name].notna().all(), name


def test_text_wins_over_numbers_in_a_column():
    merged = annatar.concat([read("csq-types-a.vcf"), read("csq-types-fallback.vcf")])
    assert isinstance(merged.var["CADD_PHRED"].dtype, pd.CategoricalDtype)
    assert merged.var.loc[SHARED, "CADD_PHRED"] == "35.5"
    assert merged.var.loc[SHARED, "STRAND"] == "-1"
    flipped = annatar.concat([read("csq-types-fallback.vcf"), read("csq-types-a.vcf")])
    assert flipped.var.loc["chr1:200:C:T", "CADD_PHRED"] == "1.5abc"
    assert flipped.var.loc["chr1:200:C:T", "STRAND"] == "plus"


def test_a_later_library_fills_what_the_first_does_not_have():
    plain = read("csq-types-a.vcf", exclude_fields=["INFO/CSQ"])
    merged = annatar.concat([plain, read("csq-types-b.vcf")])
    assert merged.var.loc[SHARED, "SYMBOL"] == "GENE1"
    assert pd.isna(merged.var.loc["chr1:200:C:T", "SYMBOL"])


def picked_differently(tmp_path, name, picked):
    source = (DATA / "csq-types-a.vcf").read_text()
    entry = "A|stop_gained|HIGH|GENE1|ENSG1|ENST1|YES|1035|916|306/393|8/11|"
    assert entry in source
    path = tmp_path / f"{name}.vcf"
    path.write_text(source.replace(entry, picked))
    return read(path)


def test_a_shared_variant_is_one_librarys_row(tmp_path):
    intronic = picked_differently(
        tmp_path, "intronic", "A|intron_variant|MODIFIER|GENE1|ENSG1|ENST1|YES|||306/393||"
    )
    coding = picked_differently(
        tmp_path, "coding", "A|stop_gained|HIGH|GENE1|ENSG1|ENST1|YES|1035|916||8/11|"
    )
    row = annatar.concat([intronic, coding]).var.loc[SHARED]
    assert (row["Consequence"], row["IMPACT"]) == ("intron_variant", "MODIFIER")
    assert row["Protein_position"] == "306/393"
    assert pd.isna(row["EXON"])
    assert row[COUNTS].tolist() == [4, 2, 24, 6]
    row = annatar.concat([coding, intronic]).var.loc[SHARED]
    assert (row["Consequence"], row["IMPACT"]) == ("stop_gained", "HIGH")
    assert row["EXON"] == "8/11"
    assert pd.isna(row["Protein_position"])
    assert row[COUNTS].tolist() == [4, 2, 24, 6]


def test_a_missing_record_value_stays_missing_on_the_first_librarys_row():
    first, second = read("csq-types-a.vcf"), read("csq-types-b.vcf")
    first.var["qual"] = np.nan
    second.var["qual"] = 2.0
    merged = annatar.concat([first, second])
    assert np.isnan(merged.var.loc[SHARED, "qual"])
    assert merged.var.loc["chr1:300:CA:C", "qual"] == 2.0


def test_the_first_library_wins_a_disagreement():
    first, second = read("csq-types-a.vcf"), read("csq-types-b.vcf")
    first.var["qual"] = 1.0
    second.var["qual"] = 2.0
    assert annatar.concat([first, second]).var.loc[SHARED, "qual"] == 1.0
    assert annatar.concat([second, first]).var.loc[SHARED, "qual"] == 2.0


def test_types_survive_the_merge():
    a, b = read("typed-fields.vcf"), read("basic.vcf")
    merged = annatar.concat([a, b])
    assert merged.X.dtype == np.int8
    assert sp.isspmatrix_csr(merged.X)
    for library in (a, b):
        for name, layer in library.layers.items():
            assert merged.layers[name].dtype == layer.dtype, name
        for name, dtype in library.var.dtypes.items():
            if not isinstance(dtype, pd.CategoricalDtype) and dtype.name != "bool":
                assert merged.var[name].dtype == dtype, name


def test_a_flag_absent_from_a_library_stays_boolean():
    merged = annatar.concat([read("typed-fields.vcf"), read("no-csq.vcf")])
    flags = [
        name for name, dtype in read("typed-fields.vcf").var.dtypes.items() if dtype.name == "bool"
    ]
    assert flags
    assert all(merged.var[name].dtype == "boolean" for name in flags)


def test_a_layer_absent_from_a_library_is_kept():
    a = read("typed-fields.vcf")
    b = read("typed-fields.vcf", exclude_fields=["FORMAT/DP"])
    assert "DP" in a.layers and "DP" not in b.layers
    merged = annatar.concat({"a": a, "b": b})
    assert "DP" in merged.layers
    assert merged[merged.obs["library"] == "a"].layers["DP"].nnz == a.layers["DP"].nnz
    assert merged[merged.obs["library"] == "b"].layers["DP"].nnz == 0


def test_inputs_are_left_alone():
    a, b = read("typed-fields.vcf"), read("no-csq.vcf")
    before = a.copy(), b.copy()
    annatar.concat([a, b])
    assert_same_anndata(a, before[0])
    assert_same_anndata(b, before[1])


def test_shared_provenance_is_kept():
    merged = annatar.concat([read("csq-types-a.vcf"), read("csq-types-b.vcf")])
    assert "x_encoding" in merged.uns
    assert "pick_rule" in merged.uns
    assert "source_vcf" not in merged.uns


def test_a_library_with_no_variants_merges():
    merged = annatar.concat([read("no-records.vcf"), read("basic.vcf")])
    assert merged.shape == (6, 5)
    assert merged.var["SYMBOL"].notna().all()


def test_numbers_of_two_types_become_floats():
    first, second = read("csq-types-a.vcf"), read("csq-types-b.vcf")
    second.var["DP"] = second.var["DP"].astype("float64") + 0.5
    merged = annatar.concat([second, first])
    assert merged.var["DP"].dtype == "float64"
    assert merged.var.loc[SHARED, "DP"] == 8.5
    assert merged.var.loc["chr1:200:C:T", "DP"] == 8.0


def test_a_column_that_cannot_share_a_type_becomes_text():
    first, second = read("csq-types-a.vcf"), read("csq-types-b.vcf")
    second.var["DP"] = [True, False]
    merged = annatar.concat([first, second])
    assert isinstance(merged.var["DP"].dtype, pd.CategoricalDtype)
    assert merged.var.loc[SHARED, "DP"] == "8"
    assert merged.var.loc["chr1:300:CA:C", "DP"] == "False"


def test_a_merge_writes_and_reads_back(tmp_path):
    libraries = [read("csq-types-a.vcf"), read("csq-types-fallback.vcf"), read("typed-fields.vcf")]
    merged = annatar.concat(libraries)
    merged.write_h5ad(tmp_path / "merged.h5ad")
    back = annatar.concat([tmp_path / "merged.h5ad"])
    assert back.shape == merged.shape
    pd.testing.assert_series_equal(
        back.var["STRAND"], merged.var["STRAND"], check_categorical=False
    )


def test_nothing_to_merge_raises():
    with pytest.raises(ValueError):
        annatar.concat([])
