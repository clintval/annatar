import hashlib
import subprocess
import warnings
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
BINARY = REPO / "target" / "release" / "annatar"
DATA = REPO / "tests" / "data"


def read(path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", message=".*future.infer_string.*")
        return ad.read_h5ad(path)


def convert(tmp_path, name, *flags):
    output = tmp_path / f"{Path(name).name}.h5ad"
    subprocess.run(
        [BINARY, "--output", output, *flags, DATA / name],
        check=True,
        capture_output=True,
    )
    return read(output)


def dense(matrix):
    return np.asarray(matrix.todense())


def assert_same_tables(left, right):
    assert (left.X != right.X).nnz == 0
    assert left.layers.keys() == right.layers.keys()
    for key in left.layers:
        assert np.array_equal(left.layers[key].indptr, right.layers[key].indptr)
        assert np.array_equal(left.layers[key].indices, right.layers[key].indices)
        assert np.array_equal(left.layers[key].data, right.layers[key].data, equal_nan=True)
    pd.testing.assert_frame_equal(left.var, right.var)
    pd.testing.assert_frame_equal(left.obs, right.obs)


def test_fixture_encodings(fixture_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", message=".*future.infer_string.*")
        a = ad.read_h5ad(fixture_path)
        assert a.shape == (3, 2) and a.X.dtype == np.int8 and sp.issparse(a.X)
        assert a.layers["AD_alt"].nnz == 2 and a.layers["AD_alt"][0, 0] == 0
        assert np.isnan(a.layers["VAF"].toarray()[0, 0])
        assert isinstance(a.var["SYMBOL"].dtype, pd.CategoricalDtype)
        assert pd.isna(a.var["SYMBOL"].iloc[1])
        assert np.isnan(a.var["REVEL_score"].iloc[1])
        assert str(a.var["AC"].dtype) == "Int64" and pd.isna(a.var["AC"].iloc[1])
        assert a.var["INDEL"].dtype == bool
        assert list(a.uns["x_encoding"])[1].startswith("1:")
        assert a.uns["source_vcf"] == "in.vcf"
    md.MuData({"variants": a})


def test_multiallelic_gives_one_column_per_alt(tmp_path):
    a = convert(tmp_path, "multiallelic.vcf")
    assert list(a.var_names) == [
        "chr1:100:C:T",
        "chr1:100:C:G",
        "chr1:200:GAAA:G",
        "chr1:200:GAAA:GA",
    ]
    t, g = a.var_names.get_loc("chr1:100:C:T"), a.var_names.get_loc("chr1:100:C:G")
    c1 = a.obs_names.get_loc("c1")
    assert dense(a.X)[c1, t] == 1 and dense(a.X)[c1, g] == 2
    assert list(a.var["Consequence"].iloc[:2]) == ["missense_variant", "stop_gained"]
    assert list(a.var["AC"].iloc[:2]) == [2, 3]
    assert list(a.var["alt_reads"].iloc[:2]) == [11, 9]


def test_multiallelic_indels_match_trimmed_vep_alleles(tmp_path):
    a = convert(tmp_path, "multiallelic.vcf")
    assert list(a.var["Allele"].iloc[2:]) == ["-", "A"]


def test_allele_number_matches_by_position_in_the_alt_list(tmp_path):
    a = convert(tmp_path, "allele-number.vcf")
    assert a.var.loc["chr1:100:CG:CT", "Allele"] == "T"
    assert a.var.loc["chr1:200:AT:A", "Consequence"] == "stop_gained"
    assert a.var.loc["chr1:200:AT:ATT", "Consequence"] == "frameshift_variant"
    assert a.var.loc["chr1:300:C:C[chr1:900[", "Allele"] == "C[chr1:900["
    assert a.var.loc["chr1:300:C:C[chr1:900[", "Consequence"] == "coding_sequence_variant"
    assert a.var.loc["chr1:400:C:<DEL>", "Allele"] == "deletion"
    assert a.var.loc["chr1:500:A:G", "Consequence"] == "missense_variant"


@pytest.mark.parametrize(
    ("variant", "allele", "consequence"),
    [
        ("chr1:100:C:<DEL>", "deletion", "feature_truncation"),
        ("chr1:200:C:<DUP>", "duplication", "feature_elongation"),
        ("chr1:300:C:<INS>", "insertion", "intron_variant"),
        ("chr1:400:C:C[chr1:900[", "C[chr1:900[", "coding_sequence_variant"),
        ("chr1:500:g:a", "A", "stop_gained"),
        ("chr1:700:A:<DEL>", "copy_number_variation", "feature_truncation"),
        ("chr1:700:A:<DUP>", "copy_number_variation", "feature_truncation"),
        ("chr1:800:a:ag", "G", "stop_gained"),
        ("chr1:800:a:at", "T", "missense_variant"),
        ("chr1:900:C:C[chr1:901[", "C[chr1:901[", "coding_sequence_variant"),
        ("chr1:900:C:C[chr1:902[", "C[chr1:902[", "intron_variant"),
        ("chr1:1000:C:<INS:ME:L1>", "LINE1_insertion", "intron_variant"),
    ],
)
def test_alleles_match_in_the_spelling_vep_writes(tmp_path, variant, allele, consequence):
    a = convert(tmp_path, "allele-spelling.vcf")
    assert a.var.loc[variant, "Allele"] == allele
    assert a.var.loc[variant, "Consequence"] == consequence


def test_a_lowercase_variant_keeps_the_case_of_the_vcf_in_its_name(tmp_path):
    a = convert(tmp_path, "allele-spelling.vcf")
    assert a.var.loc["chr1:500:g:a", "ref"] == "g"
    assert a.var.loc["chr1:500:g:a", "alt"] == "a"


@pytest.mark.parametrize(
    ("name", "variant", "allele"),
    [("symbolic-csq.vcf", "chr1:100:A:<DEL>", "deletion"), ("breakend-csq.vcf", "chr1:100:A:A[chr1:500[", "A[chr1:500[")],
)
def test_vep_true_symbolic_and_breakend_records_convert(tmp_path, name, variant, allele):
    a = convert(tmp_path, name)
    assert a.var.loc[variant, "Allele"] == allele
    assert a.var.loc[variant, "pick_tier"] == "canonical"


def test_split_snvs_pick_their_own_allele_whatever_the_stale_numbers_say(tmp_path):
    a = convert(tmp_path, "split-after-vep-snv.vcf")
    assert a.var.loc["chr1:100:G:A", "Allele"] == "A"
    assert a.var.loc["chr1:100:G:A", "Consequence"] == "stop_gained"
    assert a.var.loc["chr1:100:G:T", "Allele"] == "T"
    assert a.var.loc["chr1:100:G:T", "Consequence"] == "intron_variant"


def test_split_breakends_ignore_numbers_past_the_alt_count(tmp_path):
    a = convert(tmp_path, "split-after-vep-breakend.vcf")
    assert a.var.loc["chr1:100:C:C[chr1:901[", "Consequence"] == "coding_sequence_variant"
    assert a.var.loc["chr1:100:C:C[chr1:902[", "Consequence"] == "intron_variant"


def test_a_breakend_with_no_entry_for_its_number_keeps_the_reference_side(tmp_path):
    a = convert(tmp_path, "breakend-reference-side.vcf")
    assert a.var["Allele"].tolist() == ["C."]
    assert a.var["Consequence"].tolist() == ["feature_truncation"]


def test_minimal_output_with_no_numbers_needs_one_allele_on_the_record(tmp_path):
    a = convert(tmp_path, "minimal-no-number.vcf")
    assert a.var.loc["chr1:100:CG:CT", "Allele"] == "T"
    assert a.var.loc["chr1:300:C:<DEL>", "Allele"] == "deletion"
    breakend = a.var.loc["chr1:200:C:C[chr1:900["]
    assert pd.isna(breakend["Allele"]) and pd.isna(breakend["pick_tier"])
    assert breakend["alt_reads"] == 3


def test_a_structural_alt_with_no_entry_keeps_an_empty_row(tmp_path):
    a = convert(tmp_path, "structural-no-entry.vcf")
    assert a.var.loc["chr1:100:A:G", "Consequence"] == "missense_variant"
    assert pd.isna(a.var.loc["chr1:100:A:<DEL>", "Consequence"])
    assert a.var.loc["chr1:100:A:<DEL>", "alt_reads"] == 3
    assert a.var.loc["chr1:200:C:C[chr1:900[", "Consequence"] == "intron_variant"
    assert pd.isna(a.var.loc["chr1:200:C:C[chr1:901[", "Consequence"])


def test_haploid_call_counts_as_homozygous(tmp_path):
    a = convert(tmp_path, "haploid.vcf")
    assert dense(a.X)[:, 0].tolist() == [3, 1, 0]


ANNOTATION = {
    "Allele",
    "Consequence",
    "SYMBOL",
    "consequence_rank",
    "is_coding",
    "is_lof",
    "pick_tier",
    "worst_panel_impact",
}


def test_no_csq_has_no_csq_or_derived_columns_and_works(tmp_path):
    a = convert(tmp_path, "no-csq.vcf")
    assert a.shape == (3, 2)
    assert not ANNOTATION & set(a.var.columns)
    assert list(a.uns["csq_fields"]) == []
    assert a.var["DP"].tolist() == [18, 4]


def test_a_panel_on_a_plain_vcf_adds_no_columns(tmp_path):
    a = convert(tmp_path, "no-csq.vcf", "--genes", DATA / "panel-genes.txt")
    assert not ANNOTATION & set(a.var.columns)


def test_excluding_info_csq_matches_a_plain_vcf(tmp_path):
    a = convert(tmp_path, "basic.vcf", "--exclude-field", "INFO/CSQ")
    assert a.shape == (3, 5)
    assert not ANNOTATION & set(a.var.columns)
    assert list(a.uns["csq_fields"]) == []
    assert a.var["alt_reads"].sum() > 0


def test_no_records_has_every_sample_and_no_variants(tmp_path):
    a = convert(tmp_path, "no-records.vcf")
    assert a.shape == (3, 0)
    assert list(a.obs_names) == ["c1", "c2", "c3"]
    assert a.obs["n_sites_covered"].tolist() == [0, 0, 0]
    assert a.layers["AD_ref"].shape == (3, 0)


def test_no_gt_has_an_empty_x_and_still_has_layers(tmp_path):
    a = convert(tmp_path, "no-gt.vcf")
    assert a.X.nnz == 0
    assert dense(a.layers["AD_alt"]).tolist() == [[3], [0], [0]]
    assert a.layers["AD_alt"].nnz == 2
    assert a.obs["n_sites_covered"].tolist() == [1, 1, 0]


@pytest.mark.parametrize("name", ["basic.vcf", "multiallelic.vcf", "question.vcf"])
def test_layers_share_one_sparsity_pattern(tmp_path, name):
    a = convert(tmp_path, name)
    ad_ref, ad_alt, dp = (a.layers[k] for k in ["AD_ref", "AD_alt", "DP"])
    assert np.array_equal(ad_ref.indptr, ad_alt.indptr)
    assert np.array_equal(ad_ref.indptr, dp.indptr)
    assert np.array_equal(ad_ref.indices, ad_alt.indices)
    assert np.array_equal(ad_ref.indices, dp.indices)


@pytest.mark.parametrize("name", ["basic.vcf", "multiallelic.vcf", "question.vcf"])
def test_obs_counts_match_the_matrices(tmp_path, name):
    a = convert(tmp_path, name)
    ad_ref, ad_alt = a.layers["AD_ref"], a.layers["AD_alt"]
    assert a.obs["n_sites_covered"].tolist() == np.diff(ad_ref.indptr).tolist()
    assert a.obs["n_sites_alt"].tolist() == np.asarray((ad_alt > 0).sum(axis=1)).ravel().tolist()


def test_a_missing_ref_count_still_counts_the_cell(tmp_path):
    a = convert(tmp_path, "missing-ref-count.vcf")
    assert a.layers["AD_ref"][0].nnz == 1
    assert a.obs["n_sites_covered"].tolist() == [2, 1, 0]
    assert a.obs["n_sites_alt"].tolist() == [2, 0, 0]
    assert a.var["n_cells_covered"].tolist() == [2, 1]
    assert a.var["n_cells_alt"].tolist() == [1, 1]
    assert a.var["ref_reads"].tolist() == [7, 4]
    assert a.var["alt_reads"].tolist() == [5, 2]


def test_skipped_alts_keep_the_original_positions(tmp_path):
    a = convert(tmp_path, "symbolic.vcf")
    column = a.var_names.get_loc("chr1:200:A:G")
    assert dense(a.layers["AD_alt"])[0, column] == 3
    assert dense(a.X)[0, column] == 2
    assert a.var["AC"].iloc[column] == 1


def test_star_alts_get_no_row_and_keep_their_slots(tmp_path):
    a = convert(tmp_path, "star.vcf")
    assert list(a.var_names) == ["chr1:100:AT:A", "chr1:200:A:G"]
    assert a.var["AC"].tolist() == [1, 2]
    assert dense(a.layers["AD_alt"])[0].tolist() == [3, 3]
    assert dense(a.layers["AF"])[0].tolist() == pytest.approx([0.3, 0.3])
    assert dense(a.X)[0].tolist() == [2, 2]
    assert a.obs["n_sites_alt"].tolist() == [2, 0]
    assert a.var["Consequence"].tolist() == ["frameshift_variant", "missense_variant"]


def test_placeholders_count_toward_the_shared_base_trim(tmp_path):
    a = convert(tmp_path, "symbolic.vcf")
    assert a.var.loc["chr1:100:AC:AT", "Allele"] == "T"
    assert a.var.loc["chr1:100:AC:AT", "Consequence"] == "frameshift_variant"


def test_a_symbolic_alt_on_a_record_without_csq_keeps_its_row(tmp_path):
    a = convert(tmp_path, "partly-annotated.vcf")
    assert pd.isna(a.var.loc["chr1:300:A:<DEL>", "SYMBOL"])
    assert a.var.loc["chr1:300:A:<DEL>", "n_cells_covered"] == 2


def test_a_record_without_csq_keeps_its_row_with_empty_annotation(tmp_path):
    a = convert(tmp_path, "partly-annotated.vcf")
    assert a.var["SYMBOL"].iloc[0] == "GENE1"
    assert pd.isna(a.var["SYMBOL"].iloc[1])
    assert pd.isna(a.var["pick_tier"].iloc[1])
    assert pd.isna(a.var["consequence_rank"].iloc[1])
    assert a.var["alt_reads"].tolist() == [3, 2, 3]


def test_csq_text_stays_verbatim(tmp_path):
    a = convert(tmp_path, "csq-types-a.vcf")
    expected = {
        "PUBMED": "12345&67890",
        "ClinVar": "12360",
        "SOMATIC": "0&1&0",
        "PHENO": "1&1&0",
        "Existing_variation": "rs1&COSV2&CM3",
        "Protein_position": "306/393",
        "EXON": "8/11",
        "HGNC_ID": "HGNC:11998",
        "MAX_AF_POPS": "gnomADe_AFR",
        "Unknown_plugin": "3.5",
    }
    for name, value in expected.items():
        assert a.var[name].iloc[0] == value, name
        assert isinstance(a.var[name].dtype, pd.CategoricalDtype), name
    assert pd.isna(a.var["PUBMED"].iloc[1])


@pytest.mark.parametrize("name", ["csq-types-a.vcf", "csq-types-b.vcf"])
def test_positions_are_text_in_every_file(tmp_path, name):
    a = convert(tmp_path, name)
    for column in ["cDNA_position", "CDS_position"]:
        assert isinstance(a.var[column].dtype, pd.CategoricalDtype), column
        assert a.var[column].str.split("-").notna().all(), column


def test_csq_integers_are_nullable_integers_by_name(tmp_path):
    a = convert(tmp_path, "csq-types-a.vcf")
    for name in ["STRAND", "TSL", "DISTANCE", "GENE_PHENO", "ALLELE_NUM"]:
        assert str(a.var[name].dtype) == "Int64", name
    assert a.var["STRAND"].tolist() == [-1, 1]
    assert a.var["TSL"].iloc[0] == 1 and pd.isna(a.var["TSL"].iloc[1])
    assert a.var["DISTANCE"].isna().all()


def test_csq_floats_keep_the_largest_value(tmp_path):
    a = convert(tmp_path, "csq-types-a.vcf")
    for name in ["AF", "gnomADe_AF", "MAX_AF", "CADD_PHRED", "REVEL_score"]:
        assert a.var[name].dtype == np.float64, name
    assert a.var["gnomADe_AF"].iloc[0] == 0.000002
    assert a.var["AF"].iloc[0] == 0.25
    assert a.var["MAX_AF"].iloc[0] == 0.1
    assert np.isnan(a.var["gnomADe_AF"].iloc[1])


def test_a_float_column_with_no_values_still_compares(tmp_path):
    a = convert(tmp_path, "csq-types-a.vcf")
    assert a.var["REVEL_score"].isna().all()
    assert not (a.var["REVEL_score"] > 0.5).any()


@pytest.mark.parametrize(
    ("column", "values"),
    [
        ("STRAND", ["1", "plus"]),
        ("TSL", ["1&2", "2"]),
        ("REVEL_score", ["0.5", "nan"]),
        ("gnomADe_AF", ["0.1&0.2", "inf"]),
        ("CADD_PHRED", ["1.5", "1.5abc"]),
    ],
)
def test_a_value_failing_its_type_turns_that_column_into_verbatim_text(tmp_path, column, values):
    a = convert(tmp_path, "csq-types-fallback.vcf")
    assert a.var[column].tolist() == values
    assert isinstance(a.var[column].dtype, pd.CategoricalDtype)


def dtypes(adata):
    return adata.var.dtypes.astype(str).to_dict()


def test_every_var_column_has_one_dtype_with_and_without_missing_values(tmp_path):
    complete = convert(tmp_path, "dtypes-complete.vcf")
    missing = convert(tmp_path, "dtypes-missing.vcf")
    assert complete.var.notna().all().all()
    assert missing.var.iloc[1].isna().any()
    assert dtypes(complete) == dtypes(missing)


@pytest.mark.parametrize("name", ["dtypes-complete.vcf", "dtypes-missing.vcf"])
def test_integer_and_boolean_columns_follow_whether_they_can_be_missing(tmp_path, name):
    types = dtypes(convert(tmp_path, name))
    for column in ["consequence_rank", "DP", "AC", "STRAND", "DISTANCE"]:
        assert types[column] == "Int64", column
    for column in ["pos", "n_cells_covered", "n_cells_alt", "ref_reads", "alt_reads"]:
        assert types[column] == "int64", column
    for column in ["is_coding", "is_lof"]:
        assert types[column] == "boolean", column
    assert types["SOLO"] == "bool"


def test_is_coding_and_is_lof_are_missing_where_nothing_was_picked(tmp_path):
    a = convert(tmp_path, "dtypes-missing.vcf")
    assert a.var["is_coding"].iloc[0] and not a.var["is_lof"].iloc[0]
    assert pd.isna(a.var["is_coding"].iloc[1]) and pd.isna(a.var["is_lof"].iloc[1])


def test_nullable_columns_keep_the_var_column_order(tmp_path):
    a = convert(tmp_path, "dtypes-complete.vcf")
    names = list(a.var.columns)
    assert names[:4] == ["chrom", "pos", "ref", "alt"]
    assert names.index("consequence_rank") < names.index("is_coding") < names.index("is_lof") < names.index("pick_tier")
    assert names.index("STRAND") < names.index("DISTANCE") < names.index("REVEL_score") < names.index("consequence_rank")
    assert names[-5:] == ["DP", "AC", "SCORE", "NOTE", "SOLO"]


def test_an_empty_file_has_the_same_var_dtypes(tmp_path):
    source = (DATA / "dtypes-complete.vcf").read_text().splitlines(keepends=True)
    empty = tmp_path / "empty.vcf"
    empty.write_text("".join(line for line in source if line.startswith("#")))
    output = tmp_path / "empty.h5ad"
    subprocess.run([BINARY, "--output", output, empty], check=True, capture_output=True)
    assert dtypes(read(output)) == dtypes(convert(tmp_path, "dtypes-complete.vcf"))


def concat_inputs(tmp_path, other):
    first = convert(tmp_path, "csq-types-a.vcf")
    recounted = tmp_path / "recounted.vcf"
    recounted.write_text((DATA / other).read_text().replace("0/1:5,3", "0/1:4,6"))
    return {"a": first, "b": convert(tmp_path, recounted)}


def test_concat_of_identical_variants_keeps_the_annotation(tmp_path):
    samples = concat_inputs(tmp_path, "csq-types-a.vcf")
    same = ad.concat(samples, join="outer", merge="same", index_unique="-")
    assert set(samples["a"].var.columns) - set(same.var.columns) == {"ref_reads", "alt_reads"}
    assert dtypes(same) == {k: v for k, v in dtypes(samples["a"]).items() if k in same.var}
    first = ad.concat(samples, join="outer", merge="first", index_unique="-")
    assert first.var["alt_reads"].tolist() == samples["a"].var["alt_reads"].tolist()
    assert dtypes(first) == dtypes(samples["a"])


def test_concat_of_different_variants_cannot_keep_the_annotation_alone(tmp_path):
    samples = concat_inputs(tmp_path, "csq-types-b.vcf")
    same = ad.concat(samples, join="outer", merge="same", index_unique="-")
    assert "SYMBOL" not in same.var and "Consequence" not in same.var
    first = ad.concat(samples, join="outer", merge="first", index_unique="-")
    assert set(first.var_names) == {"chr1:100:G:A", "chr1:200:C:T", "chr1:300:CA:C"}
    assert first.var["SYMBOL"].isna().tolist() == [False, False, True]
    assert first.var["pos"].dtype == np.float64


@pytest.mark.filterwarnings("ignore:The 'future.infer_string' option")
def test_the_readme_concat_recipe_keeps_every_value(tmp_path):
    samples = concat_inputs(tmp_path, "csq-types-b.vcf")
    adata = ad.concat(samples, join="outer", label="sample", index_unique="-")
    var = pd.concat([sample.var for sample in samples.values()])
    adata.var = var[~var.index.duplicated()].reindex(adata.var_names)

    assert adata.var["cDNA_position"].tolist() == ["1035", "862", "519-520"]
    assert not isinstance(adata.var["cDNA_position"].dtype, pd.CategoricalDtype)
    assert adata.var["Consequence"].notna().all()
    expected = dtypes(samples["a"])
    for column in ["STRAND", "consequence_rank", "is_lof", "REVEL_score", "pos", "alt_reads"]:
        assert dtypes(adata)[column] == expected[column], column

    output = tmp_path / "concat.h5ad"
    adata.write_h5ad(output)
    assert read(output).var["cDNA_position"].tolist() == ["1035", "862", "519-520"]


@pytest.mark.filterwarnings("ignore:The 'future.infer_string' option")
def test_an_outer_join_on_different_variants_loses_stored_zeros(tmp_path):
    samples = concat_inputs(tmp_path, "csq-types-b.vcf")
    assert (samples["a"].layers["AD_alt"].data == 0).sum() == 2
    adata = ad.concat(samples, join="outer", index_unique="-")
    assert (adata.layers["AD_alt"].data == 0).sum() == 0
    assert adata.X.nnz == samples["a"].X.nnz + samples["b"].X.nnz


def test_an_unranked_term_sorts_after_every_known_term(tmp_path):
    a = convert(tmp_path, "unranked-consequence.vcf")
    assert a.var["Feature"].iloc[0] == "ENST2"
    assert a.var["Consequence"].iloc[0] == "intron_variant"
    assert a.var["Consequence"].iloc[1] == "invented_variant"
    assert a.var["consequence_rank"].tolist()[1:] == [42, 13, 42]


def test_an_unranked_term_is_neither_coding_nor_lof(tmp_path):
    a = convert(tmp_path, "unranked-consequence.vcf")
    assert not a.var["is_coding"].iloc[1] and not a.var["is_lof"].iloc[1]
    assert a.var["is_coding"].iloc[2]


def test_every_info_shape_and_type_reads_back(tmp_path):
    a = convert(tmp_path, "typed-fields.vcf")
    multiallelic_t, multiallelic_c = "chr1:300:G:A", "chr1:300:G:C"
    assert a.var["AC"].tolist()[:4] == [300, 70000, 1, 2]
    assert pd.isna(a.var["AC"].iloc[4])
    assert str(a.var["AC"].dtype) == "Int64"
    assert a.var["AF"].iloc[0] == 0.5 and np.isnan(a.var.loc[multiallelic_c, "AF"])
    assert a.var.loc[multiallelic_t, "RD_ref"] == 9 and a.var.loc[multiallelic_c, "RD_ref"] == 9
    assert a.var.loc[multiallelic_t, "RD_alt"] == 1 and a.var.loc[multiallelic_c, "RD_alt"] == 2
    assert a.var["QUARTET_4"].iloc[1] == 4 and pd.isna(a.var["QUARTET_1"].iloc[2])
    assert a.var["NOTE"].iloc[0] == "one,two" and pd.isna(a.var["NOTE"].iloc[4])
    assert a.var.loc[multiallelic_t, "LABEL"] == "p" and pd.isna(a.var.loc[multiallelic_c, "LABEL"])
    assert pd.isna(a.var.loc[multiallelic_t, "CODE"]) and a.var.loc[multiallelic_c, "CODE"] == "q"
    assert a.var["SOLO"].tolist() == [True] + [False] * 6
    assert a.var["id"].iloc[0] == "rs7;rs8" and a.var["filter"].iloc[2] == "q10;s50"
    assert a.var.loc[multiallelic_c, "RS_ref"] == "a" and a.var.loc[multiallelic_c, "RS_alt"] == "c"
    assert a.var.loc[multiallelic_t, "PR_1"] == "u" and a.var.loc[multiallelic_t, "PR_2"] == "v"


@pytest.mark.parametrize("name", ["typed-fields.vcf", "typed-fields.bcf"])
def test_a_lone_dot_in_array_text_info_is_missing(tmp_path, name):
    a = convert(tmp_path, name)
    row = a.var.loc["chr1:500:G:C"]
    for column in ["LABEL", "CODE", "RS_ref", "RS_alt", "PR_1", "PR_2"]:
        assert pd.isna(row[column]), column


def test_every_format_shape_and_type_reads_back(tmp_path):
    a = convert(tmp_path, "typed-fields.vcf")
    multiallelic_t, multiallelic_c = (a.var_names.get_loc(n) for n in ["chr1:300:G:A", "chr1:300:G:C"])
    ao, vaf = dense(a.layers["AO"]), dense(a.layers["VAF"])
    assert ao[1, multiallelic_t] == 2 and ao[1, multiallelic_c] == 3
    assert vaf.dtype == np.float32 and np.isnan(vaf[1, multiallelic_t])
    assert vaf[1, multiallelic_c] == pytest.approx(0.6)
    assert dense(a.layers["PAIR_2"])[0, 1] == 1
    assert dense(a.layers["RF_alt"])[0, multiallelic_c] == pytest.approx(0.9)


def test_typed_fields_match_between_vcf_and_bcf(tmp_path):
    vcf, bcf = (convert(tmp_path, n) for n in ["typed-fields.vcf", "typed-fields.bcf"])
    assert_same_tables(vcf, bcf)


@pytest.mark.parametrize("name", ["multiallelic.vcf", "symbolic.vcf", "haploid.vcf"])
def test_a_bcf_made_by_bcftools_matches_its_vcf(tmp_path, name):
    bcf = tmp_path / f"{Path(name).stem}.bcf"
    subprocess.run(
        ["bcftools", "view", "-Ob", "-o", bcf, DATA / name],
        check=True,
        capture_output=True,
    )
    assert_same_tables(convert(tmp_path, name), convert(tmp_path, bcf))


def test_mudata_constructs(tmp_path):
    md.MuData({"variants": convert(tmp_path, "basic.vcf")})


def test_every_input_format_gives_the_same_tables(tmp_path):
    frames = [convert(tmp_path, n) for n in ["basic.vcf", "basic.vcf.gz", "basic.bcf"]]
    vcf = frames[0]
    assert vcf.shape == (3, 5)
    for other in frames[1:]:
        assert_same_tables(vcf, other)


@pytest.mark.parametrize("name", ["basic.vcf", "basic.vcf.gz", "basic.bcf"])
def test_percent_encoding_is_decoded_exactly_once(tmp_path, name):
    a = convert(tmp_path, name)
    assert a.var.loc["chr1:67890:C:T", "ClinVar_CLNDN"] == "Li-Fraumeni,syndrome"
    assert a.var.loc["chr17:7674220:C:T", "ClinVar_CLNDN"] == "Odd%25name"


def test_missing_ids_read_as_an_all_missing_categorical(tmp_path):
    a = convert(tmp_path, "basic.vcf")
    assert isinstance(a.var["id"].dtype, pd.CategoricalDtype)
    assert a.var["id"].isna().all()


def test_uns_records_provenance(tmp_path):
    a = convert(tmp_path, "basic.vcf")
    assert a.uns["annatar_version"] == "0.1.0"
    assert a.uns["source_vcf"] == "basic.vcf"
    assert "--output" in a.uns["command_line"]
    assert [line.split("=")[0] for line in a.uns["vep_header"]] == ["##VEP", "##VEP-command-line"]
    assert a.uns["csq_fields"][0] == "Allele" and len(a.uns["csq_fields"]) == 19
    assert "MANE Select" in a.uns["pick_rule"]
    assert list(a.uns["x_encoding"])[1].startswith("1:")
    assert list(a.uns["layer_encoding"]) == [
        "not stored: cell not covered, or an integer value missing at a covered cell",
        "0: a measured zero",
        "NaN: a float value missing at a covered cell",
    ]
    assert list(a.uns["skipped_fields"]) == []


def test_stdin_records_a_dash_as_the_source(tmp_path):
    output = tmp_path / "out.h5ad"
    subprocess.run(
        [BINARY, "--output", output, "-"],
        input=(DATA / "basic.vcf").read_bytes(),
        check=True,
        capture_output=True,
    )
    assert read(output).uns["source_vcf"] == "-"


def test_genes_record_the_panel_and_its_digest(tmp_path):
    panel = DATA / "panel-genes.txt"
    a = convert(tmp_path, "basic.vcf", "--genes", panel)
    assert list(a.uns["panel_genes"]) == ["GENE3", "TP53"]
    assert a.uns["panel_genes_sha256"] == hashlib.sha256(panel.read_bytes()).hexdigest()
    assert a.var.loc["chr1:67890:C:T", "worst_panel_impact"] == "MODERATE"
    assert a.var.loc["chr17:7674220:C:T", "worst_panel_impact"] == "MODERATE"
    assert pd.isna(a.var.loc["chr1:12345:A:G", "worst_panel_impact"])


def test_mane_table_records_its_digest(tmp_path):
    table = DATA / "mane-summary.txt"
    a = convert(tmp_path, "basic.vcf", "--mane-table", table)
    assert a.uns["mane_table_sha256"] == hashlib.sha256(table.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    ("flag", "name", "key"),
    [
        ("--genes", "panel-genes.txt", "panel_genes_sha256"),
        ("--mane-table", "mane-summary.txt", "mane_table_sha256"),
    ],
)
def test_a_piped_lookup_file_records_the_digest_of_what_was_read(tmp_path, flag, name, key):
    source = DATA / name
    output = tmp_path / "out.h5ad"
    command = f"""'{BINARY}' --output '{output}' {flag} <(cat '{source}') '{DATA / "basic.vcf"}'"""
    subprocess.run(["bash", "-c", command], check=True, capture_output=True)
    assert read(output).uns[key] == hashlib.sha256(source.read_bytes()).hexdigest()


def test_without_the_lookup_flags_none_of_their_output_exists(tmp_path):
    a = convert(tmp_path, "basic.vcf")
    assert not {"panel_genes", "panel_genes_sha256", "mane_table_sha256"} & set(a.uns)
    assert "worst_panel_impact" not in a.var


def test_exclude_field_drops_a_csq_column(tmp_path):
    a = convert(tmp_path, "basic.vcf", "--exclude-field", "CSQ/SYMBOL")
    assert "SYMBOL" not in a.var and "Consequence" in a.var


def test_exclude_field_drops_info_and_format_fields(tmp_path):
    a = convert(tmp_path, "typed-fields.vcf", "--exclude-field", "INFO/AC", "--exclude-field", "FORMAT/DP")
    assert "AC" not in a.var and "AF" in a.var
    assert "DP" not in a.layers and "AD_alt" in a.layers


def test_the_whole_question_selects_one_variant(tmp_path):
    a = convert(tmp_path, "question.vcf")
    hits = a[:, (a.var["SYMBOL"] == "TP53") & a.var["ClinVar_CLNSIG"].isin(["Pathogenic", "Likely_pathogenic"])]
    assert list(hits.var_names) == ["chr17:7674220:C:T"]
    assert hits.X.nnz == 3


def test_the_readme_mudata_example_selects_carriers(tmp_path):
    variants = convert(tmp_path, "question.vcf")
    mu = md.MuData({"variants": variants})
    pathogenic = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
    hits = (variants.var["SYMBOL"] == "TP53") & variants.var["ClinVar_CLNSIG"].isin(pathogenic)
    carriers = np.asarray(variants[:, hits].layers["AD_alt"].sum(axis=1)).ravel() > 0
    assert carriers.shape == (mu["variants"].n_obs,)
    assert carriers.any()


def test_large_output_reads_without_a_hdf5_plugin(tmp_path, monkeypatch):
    import h5py

    cells, records = 220, 210
    alt = {}
    lines = [
        "##fileformat=VCFv4.3",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        "##contig=<ID=chr1,length=1000000>",
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + [f"c{i}" for i in range(cells)]
        ),
    ]
    for j in range(records):
        calls = []
        for i in range(cells):
            if (i + j) % 7 == 0:
                calls.append("./.:.:.")
            else:
                alt[i, j] = 1 + (i * j) % 9
                calls.append(f"0/1:{5},{alt[i, j]}:{5 + alt[i, j]}")
        lines.append(
            "\t".join(["chr1", str(100 + j), ".", "A", "G", "50", "PASS", ".", "GT:AD:DP"] + calls)
        )
    vcf = tmp_path / "large.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    output = tmp_path / "large.h5ad"
    subprocess.run([BINARY, "--output", output, vcf], check=True, capture_output=True)

    offenders = []

    def check(name, node):
        if isinstance(node, h5py.Dataset):
            filters = [f for f in node._filters if f != "gzip"]
            if filters:
                offenders.append(f"{name}: {filters}")

    with h5py.File(output) as handle:
        handle.visititems(check)
    assert not offenders, "datasets with non-gzip filters: " + "; ".join(offenders)

    monkeypatch.delenv("HDF5_PLUGIN_PATH", raising=False)
    adata = read(output)
    assert adata.shape == (cells, records)
    assert adata.X.nnz == len(alt)
    layer = adata.layers["AD_alt"].toarray()
    for (i, j), value in list(alt.items())[::97]:
        assert layer[i, j] == value


def test_floats_read_back_exactly_and_an_unvalued_csq_subfield_is_categorical(tmp_path):
    lines = [
        "##fileformat=VCFv4.3",
        "##contig=<ID=chr1,length=1000000>",
        '##INFO=<ID=SCORE,Number=A,Type=Float,Description="A score">',
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP. Format: Allele|Consequence|SYMBOL|Note">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "c1"]),
        "\t".join(
            ["chr1", "100", ".", "A", "G", "30.1", "PASS", "SCORE=0.1;CSQ=G|missense_variant|TP53|", "GT:AD", "0/1:5,3"]
        ),
    ]
    vcf = tmp_path / "floats.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    output = tmp_path / "floats.h5ad"
    subprocess.run([BINARY, "--output", output, vcf], check=True, capture_output=True)
    a = read(output)
    assert a.var["qual"].iloc[0] == 30.1
    assert a.var["SCORE"].iloc[0] == 0.1
    assert isinstance(a.var["Note"].dtype, pd.CategoricalDtype)
    assert a.var["Note"].isna().all()
    assert not a.var["Note"].str.contains("x").any()
