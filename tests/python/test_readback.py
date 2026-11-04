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
    md.MuData({"snv": a})


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


def test_haploid_call_counts_as_homozygous(tmp_path):
    a = convert(tmp_path, "haploid.vcf")
    assert dense(a.X)[:, 0].tolist() == [3, 1, 0]


def test_no_csq_has_no_csq_columns_and_works(tmp_path):
    a = convert(tmp_path, "no-csq.vcf")
    assert a.shape == (3, 2)
    assert not {"Allele", "Consequence", "SYMBOL"} & set(a.var.columns)
    assert list(a.uns["csq_fields"]) == []
    assert a.var["consequence_rank"].isna().all()
    assert a.var["DP"].tolist() == [18, 4]


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


def test_star_alts_are_columns_without_annotation(tmp_path):
    a = convert(tmp_path, "symbolic.vcf")
    star = a.var.loc["chr1:300:A:*"]
    assert pd.isna(star["SYMBOL"]) and pd.isna(star["pick_tier"])
    assert dense(a.layers["AD_alt"])[1, a.var_names.get_loc("chr1:300:A:*")] == 4


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
    md.MuData({"snv": convert(tmp_path, "basic.vcf")})


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
