import subprocess

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from conftest import BINARY
from conftest import DATA
from test_parity import assert_same_anndata

import annatar


def converted(tmp_path, name, *flags, source=None):
    output = tmp_path / f"{name}.h5ad"
    source = source or DATA / f"{name}.vcf"
    subprocess.run([BINARY, *flags, source, output], check=True, capture_output=True)
    return output


def merged_by_cli(tmp_path, *inputs, name="merged"):
    output = tmp_path / f"{name}.h5ad"
    subprocess.run([BINARY, "merge", "--output", output, *inputs], check=True, capture_output=True)
    return output


def assert_merges_agree(tmp_path, *inputs):
    written = ad.read_h5ad(merged_by_cli(tmp_path, *inputs))
    assert_same_anndata(annatar.concat(inputs), written)


def picked_differently(tmp_path, name, picked):
    source = (DATA / "csq-types-a.vcf").read_text()
    entry = "A|stop_gained|HIGH|GENE1|ENSG1|ENST1|YES|1035|916|306/393|8/11|"
    assert entry in source
    path = tmp_path / f"{name}.vcf"
    path.write_text(source.replace(entry, picked))
    return converted(tmp_path, name, source=path)


@pytest.mark.parametrize(
    "names",
    [
        ("csq-types-a", "csq-types-b"),
        ("csq-types-b", "csq-types-a"),
        ("csq-types-a", "csq-types-fallback"),
        ("csq-types-fallback", "csq-types-a"),
        ("csq-types-a", "csq-types-fallback", "csq-types-b"),
        ("basic", "multiallelic"),
        ("no-csq", "basic"),
        ("typed-fields", "no-csq"),
        ("no-csq", "typed-fields"),
        ("no-records", "basic"),
        ("basic", "no-records"),
        ("partly-annotated", "basic"),
        ("basic",),
        ("no-records",),
    ],
    ids="+".join,
)
def test_merge_equals_concat(names, tmp_path):
    assert_merges_agree(tmp_path, *(converted(tmp_path, name) for name in names))


def test_identical_variants_sort_by_contig_position_and_alleles(tmp_path):
    first = converted(tmp_path, "basic")
    second = tmp_path / "again.h5ad"
    second.write_bytes(first.read_bytes())
    assert list(annatar.concat([first, second]).var_names) == [
        "chr1:12345:A:G",
        "chr1:67890:C:G",
        "chr1:67890:C:T",
        "chr17:7674220:C:T",
        "chr17:7675000:G:A",
    ]
    assert_merges_agree(tmp_path, first, second)


def library(tmp_path, name, contigs, *variants):
    lines = [
        "##fileformat=VCFv4.3",
        *(f"##contig=<ID={contig},length=100000>" for contig in contigs),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "c1"]),
    ]
    for variant in variants:
        chrom, pos, ref, alt = variant.split(":")
        fields = [chrom, pos, ".", ref, alt, "50", "PASS", ".", "GT:AD", "0/1:5,3"]
        lines.append("\t".join(fields))
    path = tmp_path / f"{name}.vcf"
    path.write_text("\n".join(lines) + "\n")
    return path


def merged_names(tmp_path, *libraries):
    inputs = [converted(tmp_path, path.stem, source=path) for path in libraries]
    assert_merges_agree(tmp_path, *inputs)
    return list(annatar.concat(inputs).var_names)


def test_a_header_puts_chr2_before_chr10(tmp_path):
    first = library(tmp_path, "a", ["chr2", "chr10"], "chr10:5:A:G", "chr2:9:A:G")
    second = library(tmp_path, "b", ["chr2", "chr10"], "chr2:3:A:G")
    assert merged_names(tmp_path, first, second) == ["chr2:3:A:G", "chr2:9:A:G", "chr10:5:A:G"]


def test_a_header_order_need_not_be_alphabetical(tmp_path):
    contigs = ["chrY", "chr1", "chrM"]
    first = library(tmp_path, "a", contigs, "chrM:1:A:G", "chr1:2:A:G")
    second = library(tmp_path, "b", contigs, "chrY:3:A:G")
    assert merged_names(tmp_path, first, second) == ["chrY:3:A:G", "chr1:2:A:G", "chrM:1:A:G"]


def test_positions_sort_as_numbers(tmp_path):
    first = library(tmp_path, "a", ["chr1"], "chr1:100:A:G", "chr1:1000:A:G")
    second = library(tmp_path, "b", ["chr1"], "chr1:20:A:G")
    assert merged_names(tmp_path, first, second) == [
        "chr1:20:A:G",
        "chr1:100:A:G",
        "chr1:1000:A:G",
    ]


def test_ref_then_alt_break_ties(tmp_path):
    first = library(tmp_path, "a", ["chr1"], "chr1:5:C:T", "chr1:5:AT:A")
    second = library(tmp_path, "b", ["chr1"], "chr1:5:C:A")
    assert merged_names(tmp_path, first, second) == ["chr1:5:AT:A", "chr1:5:C:A", "chr1:5:C:T"]


def test_later_libraries_append_their_contigs(tmp_path):
    first = library(tmp_path, "a", ["chr2", "chr1"], "chr1:1:A:G", "chr2:2:A:G")
    second = library(tmp_path, "b", ["chr3", "chr1", "chr2"], "chr3:3:A:G", "chr1:4:A:G")
    assert merged_names(tmp_path, first, second) == [
        "chr2:2:A:G",
        "chr1:1:A:G",
        "chr1:4:A:G",
        "chr3:3:A:G",
    ]
    merged = annatar.concat([tmp_path / "a.h5ad", tmp_path / "b.h5ad"])
    assert list(merged.uns["contigs"]) == ["chr2", "chr1", "chr3"]


def test_a_contig_on_no_list_goes_after_the_listed_ones(tmp_path):
    first = annatar.read_vcf(library(tmp_path, "a", ["chr2"], "chr2:1:A:G"))
    second = annatar.read_vcf(library(tmp_path, "b", [], "chrB:1:A:G", "chrA:2:A:G"))
    del second.uns["contigs"]
    merged = annatar.concat([second, first])
    assert list(merged.var_names) == ["chr2:1:A:G", "chrB:1:A:G", "chrA:2:A:G"]


def test_a_merged_file_keeps_its_contig_order(tmp_path):
    first = library(tmp_path, "a", ["chr2", "chr10"], "chr10:5:A:G", "chr2:9:A:G")
    merged = merged_by_cli(tmp_path, converted(tmp_path, "a", source=first), name="again")
    second = library(tmp_path, "b", ["chr10", "chr2"], "chr10:1:A:G")
    names = list(annatar.concat([merged, converted(tmp_path, "b", source=second)]).var_names)
    assert names == ["chr2:9:A:G", "chr10:1:A:G", "chr10:5:A:G"]


def test_disagreeing_picks_agree(tmp_path):
    intronic = picked_differently(
        tmp_path, "intronic", "A|intron_variant|MODIFIER|GENE1|ENSG1|ENST1|YES|||306/393||"
    )
    coding = picked_differently(
        tmp_path, "coding", "A|stop_gained|HIGH|GENE1|ENSG1|ENST1|YES|1035|916||8/11|"
    )
    assert_merges_agree(tmp_path, intronic, coding)
    assert_merges_agree(tmp_path, coding, intronic)


def test_a_layer_absent_from_a_library_agrees(tmp_path):
    full = converted(tmp_path, "typed-fields")
    fewer = tmp_path / "fewer.h5ad"
    source = DATA / "typed-fields.vcf"
    subprocess.run(
        [BINARY, "--exclude-field", "FORMAT/DP", source, fewer], check=True, capture_output=True
    )
    assert_merges_agree(tmp_path, full, fewer)
    assert_merges_agree(tmp_path, fewer, full)
    assert "DP" in ad.read_h5ad(merged_by_cli(tmp_path, fewer, full)).layers


def test_a_layer_float_in_one_library_is_float(tmp_path):
    text = (DATA / "typed-fields.vcf").read_text()
    declared = "ID=AO,Number=A,Type=Integer"
    assert declared in text
    path = tmp_path / "floats.vcf"
    path.write_text(text.replace(declared, "ID=AO,Number=A,Type=Float"))
    floats = converted(tmp_path, "floats", source=path)
    plain = converted(tmp_path, "typed-fields")
    assert_merges_agree(tmp_path, floats, plain)
    assert_merges_agree(tmp_path, plain, floats)
    assert ad.read_h5ad(merged_by_cli(tmp_path, plain, floats)).layers["AO"].dtype == np.float32


def test_a_merge_is_a_valid_input(tmp_path):
    first = merged_by_cli(
        tmp_path, converted(tmp_path, "csq-types-a"), converted(tmp_path, "csq-types-b"), name="ab"
    )
    assert_merges_agree(tmp_path, first, converted(tmp_path, "basic"))
    assert_merges_agree(tmp_path, converted(tmp_path, "basic"), first)


def test_files_written_by_python_are_valid_inputs(tmp_path):
    first = tmp_path / "first.h5ad"
    second = tmp_path / "second.h5ad"
    annatar.concat([converted(tmp_path, "csq-types-a")]).write_h5ad(first)
    annatar.concat([converted(tmp_path, "typed-fields")]).write_h5ad(second)
    assert_merges_agree(tmp_path, first, second)


def test_plain_string_columns_become_text(tmp_path):
    adata = ad.read_h5ad(converted(tmp_path, "no-csq"))
    adata.var["note"] = ["one", "two"]
    path = tmp_path / "strings.h5ad"
    adata.write_h5ad(path, convert_strings_to_categoricals=False)
    assert_merges_agree(tmp_path, path, converted(tmp_path, "basic"))


def test_the_result_reads_in_h5py_and_anndata(tmp_path):
    merged = merged_by_cli(
        tmp_path, converted(tmp_path, "csq-types-a"), converted(tmp_path, "csq-types-b")
    )
    with h5py.File(merged) as h5ad:
        assert h5ad["X"].attrs["encoding-type"] == "csr_matrix"
        assert h5ad["var/n_cells_covered"].dtype == np.int64
        assert h5ad["var/consequence_rank"].attrs["encoding-type"] == "nullable-integer"
        assert h5ad["var/is_lof"].attrs["encoding-type"] == "nullable-boolean"
        assert h5ad["obs/library"].attrs["encoding-type"] == "categorical"
    adata = ad.read_h5ad(merged)
    assert list(adata.obs_names) == [
        "c1-csq-types-a",
        "c2-csq-types-a",
        "c1-csq-types-b",
        "c2-csq-types-b",
    ]
    assert adata.var.loc["chr1:100:G:A", "ref_reads"] == 24
    assert str(adata.var["consequence_rank"].dtype) == "Int64"
    assert str(adata.var["is_lof"].dtype) == "boolean"
    assert isinstance(adata.var["SYMBOL"].dtype, pd.CategoricalDtype)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda adata: adata.layers.pop("AD_ref"), "no sparse integer layer AD_ref"),
        (
            lambda adata: adata.layers.__setitem__("AD_alt", adata.layers["AD_alt"].astype(float)),
            "no sparse integer layer AD_alt",
        ),
        (lambda adata: setattr(adata, "X", adata.X.astype(float)), "no sparse integer X"),
        (lambda adata: setattr(adata, "X", adata.X.toarray()), "no sparse integer X"),
        (lambda adata: adata.var.pop("ref_reads"), "no var column ref_reads"),
        (lambda adata: adata.var.pop("chrom"), "no var column chrom"),
    ],
)
def test_an_h5ad_annatar_did_not_write_stops_with_a_clear_error(change, message, tmp_path):
    adata = ad.read_h5ad(converted(tmp_path, "no-csq"))
    change(adata)
    path = tmp_path / "changed.h5ad"
    adata.write_h5ad(path)
    other = converted(tmp_path, "basic")
    run = subprocess.run(
        [BINARY, "merge", "--output", tmp_path / "merged.h5ad", other, path],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 1
    assert "changed.h5ad is not annatar's h5ad" in run.stderr
    assert message in run.stderr
    assert "Traceback" not in run.stderr
    assert not (tmp_path / "merged.h5ad").exists()


def test_an_extra_obs_column_is_missing_for_libraries_without_it(tmp_path):
    adata = ad.read_h5ad(converted(tmp_path, "no-csq"))
    adata.obs["leiden"] = 1.5
    path = tmp_path / "annotated.h5ad"
    adata.write_h5ad(path)
    merged = ad.read_h5ad(merged_by_cli(tmp_path, path, converted(tmp_path, "basic")))
    assert merged.obs["leiden"].notna().tolist() == [True] * 3 + [False] * 3


def test_a_merged_file_gains_the_library_name_again(tmp_path):
    once = merged_by_cli(tmp_path, converted(tmp_path, "no-csq"), name="once")
    twice = ad.read_h5ad(merged_by_cli(tmp_path, once, name="twice"))
    assert twice.obs_names[0] == "c1-no-csq-once"
    assert set(twice.obs["library"]) == {"once"}
