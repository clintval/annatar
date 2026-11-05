import gzip
import os
import pty
import signal
import subprocess
import sys
import time

import anndata as ad
import h5py
import pytest
from conftest import BINARY
from conftest import DATA

from annatar import _cli

ESCAPE = "\x1b"
NO_CSQ_NOTE = "no INFO/CSQ to read, so the variants have no annotation"


def annatar(*args, stdin=None, cwd=None, no_color="1", **environment):
    env = {**os.environ, **environment}
    env.pop("RUST_LOG", None)
    env.pop("NO_COLOR", None)
    if no_color is not None:
        env["NO_COLOR"] = no_color
    run = subprocess.run(
        [BINARY, *args], input=stdin, capture_output=True, cwd=cwd, env=env, timeout=120
    )
    return subprocess.CompletedProcess(
        run.args, run.returncode, run.stdout.decode(), run.stderr.decode()
    )


def failed(*args, **options):
    run = annatar(*args, **options)
    assert run.returncode != 0
    assert "Traceback" not in run.stderr
    return run


def succeeded(*args, **options):
    run = annatar(*args, **options)
    assert run.returncode == 0, run.stderr
    assert "Traceback" not in run.stderr
    return run


def test_converts_basic(tmp_path):
    succeeded(DATA / "basic.vcf", tmp_path / "out.h5ad")
    assert (tmp_path / "out.h5ad").exists()


def test_reads_stdin(tmp_path):
    succeeded("-", tmp_path / "out.h5ad", stdin=(DATA / "basic.vcf").read_bytes())
    assert ad.read_h5ad(tmp_path / "out.h5ad").shape == (3, 5)


def test_both_paths_are_required(tmp_path):
    output = tmp_path / "out.h5ad"
    run = failed(output, stdin=(DATA / "basic.vcf").read_bytes())
    assert "<H5AD>" in run.stderr
    assert not output.exists()


@pytest.mark.parametrize("flag", ["--output", "-o"])
def test_the_output_flags_are_gone(flag, tmp_path):
    output = tmp_path / "out.h5ad"
    run = failed(flag, output, DATA / "basic.vcf")
    assert f"unexpected argument '{flag}'" in run.stderr
    assert not output.exists()


def test_threads_read_a_bgzipped_input(tmp_path):
    succeeded("-@", "2", DATA / "basic.vcf.gz", tmp_path / "out.h5ad")
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 5


def test_plain_gzip_input_is_refused_with_advice(tmp_path):
    source = tmp_path / "plain.vcf.gz"
    source.write_bytes(gzip.compress((DATA / "basic.vcf").read_bytes()))
    output = tmp_path / "out.h5ad"
    run = failed(source, output)
    assert "plain.vcf.gz is gzip but not BGZF; recompress it with bgzip" in run.stderr
    assert not output.exists()


def test_no_records_writes_empty_matrix(tmp_path):
    succeeded(DATA / "no-records.vcf", tmp_path / "out.h5ad")
    adata = ad.read_h5ad(tmp_path / "out.h5ad")
    assert adata.shape == (3, 0)
    assert list(adata.obs_names) == ["c1", "c2", "c3"]


def test_symbolic_alts_are_skipped(tmp_path):
    run = succeeded(DATA / "symbolic.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("*, <*> and <NON_REF> ALT alleles skipped: 3") == 1
    names = list(ad.read_h5ad(tmp_path / "out.h5ad").var_names)
    assert names == ["chr1:100:AC:AT", "chr1:200:A:G", "chr1:300:A:G"]


def test_a_structural_alt_with_no_entry_keeps_its_row_and_warns_once(tmp_path):
    run = succeeded(DATA / "structural-no-entry.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("symbolic or breakend ALT alleles with no CSQ entry: 2") == 1
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 4


@pytest.mark.parametrize(
    "name",
    [
        "no-csq-entry.vcf",
        "no-csq-entry-allele-number.vcf",
        "split-after-vep.vcf",
        "split-after-vep-allele-number.vcf",
        "split-collide.vcf",
        "pick-then-split.vcf",
        "single-alt-mismatch.vcf",
    ],
)
def test_no_error_suggests_allele_number(name, tmp_path):
    run = failed(DATA / "errors" / name, tmp_path / "out.h5ad")
    assert "--allele_number" not in run.stderr


@pytest.mark.parametrize("name", ["errors/minimal-subfield.vcf", "errors/minimal-command-line.vcf"])
def test_excluding_info_csq_lets_minimal_output_through(name, tmp_path):
    succeeded("--exclude-field", "INFO/CSQ", DATA / name, tmp_path / "out.h5ad")
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 1


def test_excluding_info_csq_lets_a_split_header_through(tmp_path):
    succeeded(
        "--exclude-field", "INFO/CSQ", DATA / "errors/split-norm-line.vcf", tmp_path / "out.h5ad"
    )


def test_a_split_breakend_does_not_take_its_sibling_by_number(tmp_path):
    run = succeeded(DATA / "pick-split-breakend.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("symbolic or breakend ALT alleles with no CSQ entry: 1") == 1


def test_a_path_holding_the_flag_is_not_minimal_output(tmp_path):
    succeeded(DATA / "minimal-in-a-path.vcf", tmp_path / "out.h5ad")


@pytest.mark.parametrize(
    "name",
    [
        "symbolic-csq.vcf",
        "breakend-csq.vcf",
        "breakend-reference-side.vcf",
        "breakend-reference-side-no-number.vcf",
        "allele-number.vcf",
        "allele-spelling.vcf",
        "gvcf-split.vcf",
        "mixed-placement-breakends.vcf",
        "norm-before-vep.vcf",
        "norm-without-split.vcf",
        "norm-join-after-vep.vcf",
    ],
)
def test_alleles_vep_spells_its_own_way_match_without_warning(name, tmp_path):
    run = succeeded(DATA / name, tmp_path / "out.h5ad")
    assert "no CSQ entry" not in run.stderr


SPLIT = "so it was split after VEP; annotate after splitting"
NORM = "the header shows bcftools norm splitting multiallelics after VEP; annotate after splitting"


@pytest.mark.parametrize(
    ("name", "locus", "message"),
    [
        (
            "no-ad.vcf",
            "",
            "the VCF declares no FORMAT/AD; with bcftools mpileup, add --annotate FORMAT/AD",
        ),
        ("ad-number-dot.vcf", "", "FORMAT/AD must be declared Number=R, Type=Integer"),
        ("ad-length.vcf", "chr1:100", "FORMAT/AD has 3 values but this record needs 2"),
        ("ad-length-uncovered.vcf", "chr1:100", "FORMAT/AD has 3 values but this record needs 2"),
        (
            "ad-length-counts-symbolic.vcf",
            "chr1:100",
            "FORMAT/AD has 2 values but this record needs 3",
        ),
        (
            "csq-subfields.vcf",
            "chr1:100",
            "a CSQ entry has 3 subfields where the header declares 7",
        ),
        (
            "no-csq-entry.vcf",
            "chr1:100",
            "no CSQ entry for ALT T among the Allele values G, C; annotate after splitting",
        ),
        (
            "no-csq-entry-allele-number.vcf",
            "chr1:100",
            "no CSQ entry for ALT T among the Allele values G; annotate after splitting",
        ),
        (
            "split-after-vep.vcf",
            "chr1:100",
            f"a record with one ALT names the CSQ alleles G, AT, {SPLIT}",
        ),
        (
            "split-after-vep-allele-number.vcf",
            "chr1:100",
            f"a record with one ALT names the CSQ alleles G, AT, {SPLIT}",
        ),
        (
            "split-after-vep-snv.vcf",
            "chr1:100",
            f"a record with one ALT names the CSQ alleles A, T, {SPLIT}",
        ),
        (
            "split-after-vep-breakend.vcf",
            "chr1:100",
            "a record with one ALT names the CSQ alleles C[chr1:901[, C[chr1:902[, so it was split",
        ),
        (
            "split-collide.vcf",
            "chr1:100",
            f"a record with one ALT names the CSQ alleles CT, T, {SPLIT}",
        ),
        (
            "pick-then-split.vcf",
            "chr1:100",
            "no CSQ entry for ALT T among the Allele values A; annotate after splitting",
        ),
        (
            "split-stale-number.vcf",
            "chr1:100",
            f"ALLELE_NUM 2 is past this record's ALT count of 1, {SPLIT}",
        ),
        ("split-norm-line.vcf", "", NORM),
        ("split-norm-line-interleaved.vcf", "", NORM),
        (
            "single-alt-mismatch.vcf",
            "chr1:100",
            "no CSQ entry for ALT G among the Allele values T; annotate after splitting",
        ),
        ("minimal-subfield.vcf", "", "VEP --minimal output is not supported"),
        ("minimal-command-line.vcf", "", "VEP --minimal output is not supported"),
        ("minimal-breakend.vcf", "chr1:100", "VEP --minimal output is not supported"),
        (
            "layer-collides-with-ad.vcf",
            "",
            "FORMAT/AD and FORMAT/AD_ref both make the layer AD_ref",
        ),
        ("layer-collides-number-r.vcf", "", "FORMAT/X and FORMAT/X_ref both make the layer X_ref"),
        ("layer-collides-fixed.vcf", "", "FORMAT/X and FORMAT/X_1 both make the layer X_1"),
        ("var-column-collides-info.vcf", "", "two fields make the var column X_ref"),
        ("var-column-collides-csq.vcf", "", "two fields make the var column pos"),
        ("csq-no-format.vcf", "", 'INFO/CSQ has no "Format: " list in its Description'),
        (
            "gt-allele-range.vcf",
            "chr1:100",
            "FORMAT/GT calls allele 3 but the last allele here is 1",
        ),
        ("negative-ad-uncovered.vcf", "chr1:100", "FORMAT/AD holds a negative count, -5"),
        ("negative-ad-skipped-alt.vcf", "chr1:100", "FORMAT/AD holds a negative count, -1"),
        ("bad-info-value.vcf", "chr1:100", "invalid digit"),
        ("bad-format-value.vcf", "chr1:100", "invalid digit"),
        ("duplicate-sample.vcf", "", "sample c1 appears more than once"),
        (
            "duplicate-variant.vcf",
            "chr1:100",
            "variant chr1:100:A:G appears more than once; "
            "remove repeats with bcftools norm -d exact",
        ),
        ("info-number-a-length.vcf", "chr1:100", "INFO/AC has 2 values but this record needs 1"),
        ("info-number-r-length.vcf", "chr1:100", "INFO/RD has 1 values but this record needs 2"),
        ("info-fixed-length.vcf", "chr1:100", "INFO/QUARTET has 3 values but this record needs 4"),
        ("negative-ad.vcf", "chr1:100", "FORMAT/AD holds a negative count, -1"),
        (
            "format-number-a-length.vcf",
            "chr1:100",
            "FORMAT/AF has 2 values but this record needs 1",
        ),
    ],
)
def test_an_error_input_fails_naming_the_problem_and_writes_nothing(name, locus, message, tmp_path):
    run = failed(DATA / "errors" / name, tmp_path / "out.h5ad")
    assert run.returncode == 1
    assert locus in run.stderr
    assert message in run.stderr
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("name", "field", "message"),
    [
        ("basic.vcf", "FORMAT/AD", "FORMAT/AD is required"),
        ("basic.vcf", "FORMAT/GT", "FORMAT/GT cannot be excluded"),
        (
            "basic.vcf",
            "INFO/CSQ/symbol",
            "--exclude-field INFO/CSQ/symbol names nothing in this VCF's header",
        ),
        (
            "basic.vcf",
            "INFO/SYMBOL",
            "--exclude-field INFO/SYMBOL names nothing in this VCF's header",
        ),
        ("basic.vcf", "FORMAT/dp", "--exclude-field FORMAT/dp names nothing in this VCF's header"),
        (
            "no-csq.vcf",
            "INFO/SYMBOL",
            "--exclude-field INFO/SYMBOL names nothing in this VCF's header",
        ),
    ],
)
def test_an_exclusion_that_cannot_be_honored_is_an_error(name, field, message, tmp_path):
    run = failed("--exclude-field", field, DATA / name, tmp_path / "out.h5ad")
    assert message in run.stderr
    assert list(tmp_path.iterdir()) == []


def test_a_header_with_no_csq_gets_one_info_line_and_no_warning(tmp_path):
    run = succeeded(DATA / "no-csq.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count(NO_CSQ_NOTE) == 1
    assert "WARN" not in run.stderr


def test_a_header_with_csq_gets_no_note(tmp_path):
    run = succeeded(DATA / "basic.vcf", tmp_path / "out.h5ad")
    assert NO_CSQ_NOTE not in run.stderr


@pytest.mark.parametrize(
    ("flag", "name"), [("--genes", "panel-genes.txt"), ("--mane-table", "mane-summary.txt")]
)
def test_a_lookup_file_warns_when_there_is_no_csq_to_act_on(flag, name, tmp_path):
    run = succeeded(flag, DATA / name, DATA / "no-csq.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count(f"{flag} does nothing: there is no INFO/CSQ to act on") == 1


@pytest.mark.parametrize("field", ["INFO/CSQ/SYMBOL", "INFO/CSQ"])
def test_an_exclusion_of_absent_annotation_warns(field, tmp_path):
    run = succeeded("--exclude-field", field, DATA / "no-csq.vcf", tmp_path / "out.h5ad")
    message = f"--exclude-field {field} does nothing: there is no INFO/CSQ to act on"
    assert run.stderr.count(message) == 1


def test_excluding_info_csq_acts_as_a_header_with_no_csq(tmp_path):
    run = succeeded(
        "--exclude-field",
        "INFO/CSQ",
        "--exclude-field",
        "INFO/CSQ/no_such_subfield",
        "--genes",
        DATA / "panel-genes.txt",
        DATA / "basic.vcf",
        tmp_path / "out.h5ad",
    )
    assert run.stderr.count(NO_CSQ_NOTE) == 1
    assert "--exclude-field INFO/CSQ does nothing" not in run.stderr
    assert run.stderr.count("--exclude-field INFO/CSQ/no_such_subfield does nothing") == 1
    assert run.stderr.count("--genes does nothing") == 1


def test_excluding_info_csq_allows_a_csq_header_with_no_format_list(tmp_path):
    succeeded(
        "--exclude-field", "INFO/CSQ", DATA / "errors/csq-no-format.vcf", tmp_path / "out.h5ad"
    )


def test_a_malformed_exclude_field_is_an_error(tmp_path):
    run = failed("--exclude-field", "SYMBOL", DATA / "basic.vcf", tmp_path / "out.h5ad")
    assert "expected INFO/ID, FORMAT/ID or INFO/CSQ/subfield" in run.stderr


@pytest.mark.parametrize(
    ("flag", "name"), [("--genes", "no-such-genes.txt"), ("--mane-table", "no-such-mane.txt")]
)
def test_a_missing_lookup_file_is_an_error_naming_the_path(flag, name, tmp_path):
    missing = tmp_path / name
    run = failed(flag, missing, DATA / "basic.vcf", tmp_path / "out.h5ad")
    assert str(missing) in run.stderr


def test_an_unwritable_output_is_an_error_naming_the_path(tmp_path):
    output = tmp_path / "no-such-directory" / "out.h5ad"
    run = failed(DATA / "basic.vcf", output)
    assert str(output) in run.stderr


def test_a_failed_run_keeps_an_earlier_output(tmp_path):
    output = tmp_path / "out.h5ad"
    output.write_text("earlier")
    failed(DATA / "errors/ad-length.vcf", output)
    assert output.read_text() == "earlier"


def test_a_successful_run_replaces_an_earlier_output(tmp_path):
    output = tmp_path / "out.h5ad"
    output.write_text("earlier")
    succeeded(DATA / "basic.vcf", output)
    assert ad.read_h5ad(output).n_vars == 5
    assert list(tmp_path.iterdir()) == [output]


def test_a_failed_write_keeps_an_earlier_output(tmp_path):
    output = tmp_path / "out.h5ad"
    output.write_text("earlier")
    (tmp_path / "out.h5ad.partial").mkdir()
    run = failed(DATA / "basic.vcf", output)
    assert str(output) in run.stderr
    assert output.read_text() == "earlier"


def test_a_failed_rename_leaves_no_partial_file(tmp_path):
    output = tmp_path / "out.h5ad"
    output.mkdir()
    run = failed(DATA / "basic.vcf", output)
    assert str(output) in run.stderr
    assert not (tmp_path / "out.h5ad.partial").exists()


def test_output_names_differing_by_suffix_do_not_share_a_partial_file(tmp_path):
    (tmp_path / "cells.h5ad.partial").mkdir()
    output = tmp_path / "cells.1"
    succeeded(DATA / "basic.vcf", output)
    assert ad.read_h5ad(output).n_vars == 5


def test_standard_output_is_refused(tmp_path):
    run = failed(DATA / "basic.vcf", "-", cwd=tmp_path)
    assert "HDF5 cannot stream" in run.stderr
    assert list(tmp_path.iterdir()) == []


def test_a_non_utf8_path_is_an_ordinary_error(tmp_path):
    run = failed(b"in\xff.vcf", tmp_path / "out.h5ad")
    assert run.returncode == 1
    assert "opening" in run.stderr


def test_records_with_no_alt_warn_once(tmp_path):
    run = succeeded(DATA / "no-alt.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("records with no ALT skipped: 1") == 1
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 1


def test_a_mane_table_is_ignored_with_a_warning_when_the_header_declares_mane_fields(tmp_path):
    run = succeeded(
        "--mane-table", DATA / "mane-summary.txt", DATA / "basic.vcf", tmp_path / "out.h5ad"
    )
    assert run.stderr.count("--mane-table is ignored") == 1


def test_no_gt_warns_once(tmp_path):
    run = succeeded(DATA / "no-gt.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("no FORMAT/GT") == 1


def test_missing_tie_fields_warn_once_with_csq(tmp_path):
    run = succeeded(DATA / "multiallelic.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("MANE_SELECT") == 1


def test_missing_tie_fields_stay_quiet_without_csq(tmp_path):
    run = succeeded(DATA / "no-csq.vcf", tmp_path / "out.h5ad")
    assert "MANE_SELECT" not in run.stderr


def test_missing_values_warn_once_per_field(tmp_path):
    run = succeeded(DATA / "missing-ref-count.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("AD_ref: missing in 1 of the cells with reads") == 1


def test_missing_values_count_only_at_cells_with_reads(tmp_path):
    run = succeeded(DATA / "copies.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("DP: missing in 2 of the cells with reads") == 1
    assert "AD_ref: missing" not in run.stderr
    assert "AD_alt: missing" not in run.stderr


def test_a_missing_alt_count_is_counted_when_the_reads_are_on_another_alt(tmp_path):
    run = succeeded(DATA / "missing-alt-reads-elsewhere.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("AD_alt: missing in 1 of the cells with reads") == 1
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 2


def test_a_value_failing_its_csq_type_warns_once_per_column(tmp_path):
    run = succeeded(DATA / "csq-types-fallback.vcf", tmp_path / "out.h5ad")
    for message in (
        "INFO/CSQ/STRAND: plus is not an integer, so the column is text in this file",
        "INFO/CSQ/TSL: 1&2 is not an integer",
        "INFO/CSQ/REVEL_score: nan is not a number",
        "INFO/CSQ/gnomADe_AF: inf is not a number",
        "INFO/CSQ/CADD_PHRED: 1.5abc is not a number",
    ):
        assert run.stderr.count(message) == 1


def test_an_unranked_consequence_warns_once_per_term(tmp_path):
    run = succeeded(DATA / "unranked-consequence.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("VEP 116 does not rank the Consequence term invented_variant") == 1
    assert run.stderr.count("Consequence term other_invented_variant") == 1


def test_another_vep_major_warns_once_and_runs(tmp_path):
    run = succeeded(DATA / "vep-115.vcf", tmp_path / "out.h5ad")
    assert run.stderr.count("annotated by VEP 115") == 1
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 5


def test_a_header_with_no_vep_line_gets_no_version_warning(tmp_path):
    run = succeeded(DATA / "multiallelic.vcf", tmp_path / "out.h5ad")
    assert "annotated by VEP" not in run.stderr


def test_a_clean_run_has_no_warnings(tmp_path):
    run = succeeded(DATA / "basic.vcf", tmp_path / "out.h5ad")
    assert "WARN" not in run.stderr


@pytest.mark.parametrize("no_color", [None, "", "1"])
def test_piped_help_has_no_escapes(no_color):
    run = succeeded("--help", no_color=no_color)
    assert "--genes" in run.stdout
    assert ESCAPE not in run.stdout


@pytest.mark.parametrize("no_color", [None, "", "1"])
def test_a_piped_usage_error_has_no_escapes(no_color):
    run = failed("--bogus", no_color=no_color)
    assert "--bogus" in run.stderr
    assert ESCAPE not in run.stderr


@pytest.mark.parametrize("no_color", [None, "", "1"])
def test_piped_help_for_no_arguments_has_no_escapes(no_color):
    run = failed(no_color=no_color)
    assert "--genes" in run.stderr
    assert ESCAPE not in run.stderr


@pytest.mark.parametrize("no_color", [None, "", "1"])
def test_a_failing_run_has_no_escapes_in_a_pipe(no_color, tmp_path):
    run = failed(DATA / "errors/ad-length.vcf", tmp_path / "out.h5ad", no_color=no_color)
    assert "ERROR" in run.stderr
    assert ESCAPE not in run.stderr


def test_forcing_color_in_the_environment_does_not_color_a_pipe():
    run = failed("--bogus", no_color=None, CLICOLOR_FORCE="1")
    assert ESCAPE not in run.stderr


def test_the_help_lists_every_flag():
    help_text = succeeded("--help").stdout
    for expected in (
        "<VCF>",
        "<H5AD>",
        "--genes",
        "--mane-table",
        "--exclude-field",
        "--threads",
        "MIT License",
        "number of copies of this ALT in the called genotype",
    ):
        assert expected in help_text
    for unexpected in (
        "-o, --output",
        "\n      --output",
        "no copies of this ALT",
        "VEP-annotated multi-sample",
        "`",
    ):
        assert unexpected not in help_text


def on_a_terminal(*args, no_color=None):
    env = {key: value for key, value in os.environ.items() if key not in ("NO_COLOR", "RUST_LOG")}
    if no_color is not None:
        env["NO_COLOR"] = no_color
    parent, child = pty.openpty()
    with subprocess.Popen(
        [BINARY, *args], stdin=child, stdout=child, stderr=child, env=env
    ) as process:
        os.close(child)
        chunks = []
        try:
            while chunk := os.read(parent, 65536):
                chunks.append(chunk)
        except OSError:
            pass
        process.wait(timeout=60)
    os.close(parent)
    return b"".join(chunks).decode()


@pytest.mark.parametrize("args", [["--help"], ["--bogus"], ["merge", "--help"]], ids=" ".join)
def test_a_terminal_gets_color_unless_no_color_says_otherwise(args):
    assert ESCAPE in on_a_terminal(*args)
    assert ESCAPE in on_a_terminal(*args, no_color="")
    assert ESCAPE not in on_a_terminal(*args, no_color="1")


def test_a_terminal_gets_colored_logging_unless_no_color_says_otherwise(tmp_path):
    source, output = DATA / "errors/ad-length.vcf", tmp_path / "out.h5ad"
    assert ESCAPE in on_a_terminal(source, output)
    assert ESCAPE not in on_a_terminal(source, output, no_color="1")


def test_a_closed_pipe_is_quiet():
    process = subprocess.Popen([BINARY, "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.stdout.close()
    _, stderr = process.communicate(timeout=60)
    assert stderr == b""
    assert process.returncode == 0


def test_the_output_carries_the_real_command_line(tmp_path):
    output = tmp_path / "out.h5ad"
    succeeded("--threads", "2", DATA / "basic.vcf", output)
    assert ad.read_h5ad(output).uns["command_line"] == (
        f"{BINARY} --threads 2 {DATA / 'basic.vcf'} {output}"
    )


def test_each_warning_is_one_log_line_and_python_adds_none(tmp_path):
    run = succeeded(DATA / "vep-115.vcf", tmp_path / "out.h5ad")
    lines = run.stderr.splitlines()
    assert len(lines) == 2, run.stderr
    assert " WARN  annatar::convert] annotated by VEP 115, but" in lines[0]
    assert " INFO  annatar] wrote 3 cells and 5 variants to " in lines[1]


def test_help_stays_quick_because_nothing_heavy_is_imported():
    code = "import sys, annatar._cli; print({'anndata', 'pandas', 'numpy'} & set(sys.modules))"
    assert subprocess.check_output([sys.executable, "-c", code], text=True).strip() == "set()"


def converted(tmp_path, name):
    output = tmp_path / f"{name}.h5ad"
    succeeded(DATA / f"{name}.vcf", output)
    return output


def test_merging_needs_an_output(tmp_path):
    run = failed("merge", converted(tmp_path, "basic"))
    assert "--output" in run.stderr


@pytest.mark.parametrize("flag", ["--output", "-o"])
def test_a_merge_output_must_be_a_file(flag, tmp_path):
    run = failed("merge", flag, "-", converted(tmp_path, "basic"))
    assert "HDF5 cannot stream" in run.stderr


def test_merging_needs_an_input(tmp_path):
    run = failed("merge", "--output", tmp_path / "out.h5ad")
    assert "<H5AD>" in run.stderr


def test_one_input_merges_with_its_name_on_every_barcode(tmp_path):
    output = tmp_path / "out.h5ad"
    run = succeeded("merge", "--output", output, converted(tmp_path, "basic"))
    assert "wrote 3 cells and 5 variants" in run.stderr
    adata = ad.read_h5ad(output)
    assert adata.shape == (3, 5)
    assert adata.obs_names[0] == "c1-basic"


def test_a_missing_input_is_named_and_writes_nothing(tmp_path):
    output = tmp_path / "out.h5ad"
    missing = tmp_path / "missing.h5ad"
    run = failed("merge", "--output", output, converted(tmp_path, "basic"), missing)
    assert f"opening {missing}" in run.stderr
    assert not output.exists()


def test_an_output_that_is_an_input_is_refused_and_left_alone(tmp_path):
    first = converted(tmp_path, "basic")
    second = converted(tmp_path, "no-csq")
    before = first.read_bytes()
    run = failed("merge", "--output", tmp_path / "." / "basic.h5ad", first, second)
    assert "is also an input" in run.stderr
    assert first.read_bytes() == before


def test_a_partial_file_that_is_an_input_is_refused(tmp_path):
    output = tmp_path / "out.h5ad"
    partial = tmp_path / "out.h5ad.partial"
    partial.write_bytes(converted(tmp_path, "basic").read_bytes())
    before = partial.read_bytes()
    run = failed("merge", "--output", output, partial)
    assert "is also an input" in run.stderr
    assert partial.read_bytes() == before


def test_two_libraries_with_one_name_are_refused(tmp_path):
    first = converted(tmp_path, "basic")
    (tmp_path / "elsewhere").mkdir()
    second = tmp_path / "elsewhere" / "basic.h5ad"
    second.write_bytes(first.read_bytes())
    output = tmp_path / "out.h5ad"
    run = failed("merge", "--output", output, first, second)
    assert (
        "more than one library is named basic; "
        "rename an input or pass annatar.concat a mapping of names"
    ) in run.stderr
    assert not output.exists()


def test_a_failed_merge_keeps_an_earlier_output(tmp_path):
    first = converted(tmp_path, "basic")
    output = tmp_path / "out.h5ad"
    output.write_text("earlier")
    failed("merge", "--output", output, first, tmp_path / "missing.h5ad")
    assert output.read_text() == "earlier"


def test_a_file_that_is_not_hdf5_stops_with_a_clear_error(tmp_path):
    source = tmp_path / "text.h5ad"
    source.write_text("not hdf5")
    run = failed("merge", "--output", tmp_path / "out.h5ad", source)
    assert run.returncode == 1
    assert "text.h5ad is not annatar's h5ad" in run.stderr


def test_an_h5ad_from_elsewhere_stops_with_a_clear_error(tmp_path):
    source = tmp_path / "empty.h5ad"
    ad.AnnData().write_h5ad(source)
    run = failed("merge", "--output", tmp_path / "out.h5ad", source)
    assert run.returncode == 1
    assert "empty.h5ad is not annatar's h5ad" in run.stderr


def test_an_input_named_merge_works_with_a_path(tmp_path):
    (tmp_path / "merge").write_bytes((DATA / "basic.vcf").read_bytes())
    succeeded("./merge", "out.h5ad", cwd=tmp_path)
    assert ad.read_h5ad(tmp_path / "out.h5ad").n_vars == 5


def test_an_input_named_help_is_still_an_input(tmp_path):
    (tmp_path / "help").write_bytes((DATA / "basic.vcf").read_bytes())
    succeeded("help", "out.h5ad", cwd=tmp_path)


def test_flags_of_the_default_command_do_not_mix_with_merge(tmp_path):
    failed(
        "--threads", "2", "merge", "--output", tmp_path / "out.h5ad", converted(tmp_path, "basic")
    )


def test_help_shows_merge():
    help_text = succeeded("--help").stdout
    assert "merge" in help_text
    assert "annatar merge --output" in help_text


def test_merge_help_is_decorated_and_plain_when_piped():
    help_text = succeeded("merge", "--help").stdout
    assert "MIT License" in help_text
    assert "--output" in help_text
    assert ESCAPE not in help_text


def test_a_bare_merge_prints_its_help():
    assert "MIT License" in failed("merge").stderr


def test_a_merge_usage_error_has_no_escapes_when_piped():
    assert ESCAPE not in failed("merge", "--bogus").stderr


def encodings(path):
    found = set()
    with h5py.File(path) as h5ad:
        h5ad.visititems(lambda _, item: found.add(item.attrs.get("encoding-type")))
    return found


def test_strings_use_the_encoding_anndata_0_10_reads(tmp_path):
    converted_file = converted(tmp_path, "basic")
    merged = tmp_path / "merged.h5ad"
    succeeded("merge", "--output", merged, converted_file)
    for path in (converted_file, merged):
        found = encodings(path)
        assert "nullable-string-array" not in found
        with h5py.File(path) as h5ad:
            for frame in (h5ad["obs"], h5ad["var"]):
                index = frame[frame.attrs["_index"]]
                assert index.attrs["encoding-type"] == "string-array"


class InterruptedWrite:
    n_obs, n_vars = 3, 5

    def write_h5ad(self, path, **_):
        path.write_bytes(b"half a file")
        raise KeyboardInterrupt


def test_an_interrupted_write_leaves_no_partial_file(tmp_path):
    output = tmp_path / "out.h5ad"
    output.write_text("earlier")
    with pytest.raises(KeyboardInterrupt):
        _cli._write(InterruptedWrite(), output)
    assert list(tmp_path.iterdir()) == [output]
    assert output.read_text() == "earlier"


def test_an_interrupted_conversion_stops_at_once(tmp_path):
    fifo = tmp_path / "cells.vcf"
    os.mkfifo(fifo)
    process = subprocess.Popen(
        [BINARY, fifo, tmp_path / "out.h5ad"], stderr=subprocess.PIPE, text=True
    )
    time.sleep(1)
    process.send_signal(signal.SIGINT)
    _, stderr = process.communicate(timeout=10)
    assert process.returncode != 0
    assert "Traceback" not in stderr
    assert not (tmp_path / "out.h5ad").exists()
    assert not (tmp_path / "out.h5ad.partial").exists()
