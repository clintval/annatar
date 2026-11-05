import logging
import os
import subprocess
import sys
import warnings
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pytest
from conftest import DATA
from conftest import logged
from conftest import read_quietly

import annatar

FIXTURES = sorted(
    path for path in DATA.iterdir() if path.name.endswith((".vcf", ".vcf.gz", ".bcf"))
)
REFUSED = sorted((DATA / "errors").iterdir()) + [DATA / "sites-only.vcf"]

RELEASES_THE_GIL = """
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import annatar

fifo, source = map(Path, sys.argv[1:])
with ThreadPoolExecutor(1) as pool:
    conversion = pool.submit(annatar.read_vcf, fifo)
    fifo.write_bytes(source.read_bytes())
    assert conversion.result().n_vars == 2
"""


REIMPORTS_THE_EXTENSION = """
import importlib
import sys
import warnings

import annatar

del sys.modules["annatar._native"]
with warnings.catch_warnings():
    warnings.simplefilter("error", RuntimeWarning)
    importlib.import_module("annatar._native")
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    annatar.read_vcf(sys.argv[1])
assert any("annotated by VEP 115" in str(warning.message) for warning in caught), caught
"""


@pytest.mark.parametrize("path", REFUSED, ids=lambda path: path.name)
def test_errors_carry_the_cli_message(path, run_cli):
    run = run_cli(path)
    assert run.returncode == 1
    with pytest.raises(ValueError) as caught:
        read_quietly(path)
    assert [str(caught.value)] == logged(run, "ERROR")


def test_a_missing_file_raises_with_the_cli_message(tmp_path, run_cli):
    path = tmp_path / "missing.vcf"
    run = run_cli(path)
    with pytest.raises(FileNotFoundError) as caught:
        annatar.read_vcf(path)
    assert [str(caught.value)] == logged(run, "ERROR")


def test_an_invalid_field_name_raises():
    with pytest.raises(ValueError, match="expected INFO/ID, FORMAT/ID or INFO/CSQ/subfield"):
        annatar.read_vcf(DATA / "basic.vcf", exclude_fields=["CSQ"])


def test_a_lone_string_of_fields_raises():
    with pytest.raises(TypeError):
        annatar.read_vcf(DATA / "basic.vcf", exclude_fields="INFO/CSQ")


def test_zero_threads_raises():
    with pytest.raises(ValueError):
        annatar.read_vcf(DATA / "basic.vcf", threads=0)


@pytest.mark.parametrize("path", FIXTURES + REFUSED, ids=lambda path: path.name)
def test_warnings_match_the_cli_log(path, run_cli):
    run = run_cli(path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            annatar.read_vcf(path)
        except ValueError:
            assert run.returncode == 1
    assert [str(warning.message) for warning in caught] == logged(run, "WARN")
    assert all(warning.category is UserWarning for warning in caught)


def test_a_vcf_with_no_csq_logs_an_info_line(caplog):
    with caplog.at_level(logging.INFO, logger="annatar"), warnings.catch_warnings():
        warnings.simplefilter("error")
        annatar.read_vcf(DATA / "no-csq.vcf")
    assert caplog.messages == ["no INFO/CSQ to read, so the variants have no annotation"]


def test_an_annotated_vcf_logs_nothing(caplog):
    with caplog.at_level(logging.INFO, logger="annatar"):
        annatar.read_vcf(DATA / "basic.vcf")
    assert caplog.messages == []


def test_a_warning_points_at_the_caller():
    with pytest.warns(UserWarning, match="annotated by VEP 115") as caught:
        annatar.read_vcf(DATA / "vep-115.vcf")
    assert Path(caught[0].filename) == Path(__file__)


def test_a_warning_filter_can_stop_the_read():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(UserWarning, match="annotated by VEP 115"):
            annatar.read_vcf(DATA / "vep-115.vcf")


def test_nothing_reaches_stderr(capfd):
    read_quietly(DATA / "vep-115.vcf")
    read_quietly(DATA / "no-csq.vcf")
    captured = capfd.readouterr()
    assert (captured.out, captured.err) == ("", "")


def test_reimporting_the_extension_keeps_warnings_without_an_alarm():
    command = [sys.executable, "-c", REIMPORTS_THE_EXTENSION, DATA / "vep-115.vcf"]
    subprocess.run(command, check=True, timeout=60)


def test_the_gil_is_released_during_conversion(tmp_path):
    fifo = tmp_path / "cells.vcf"
    os.mkfifo(fifo)
    command = [sys.executable, "-c", RELEASES_THE_GIL, fifo, DATA / "no-csq.vcf"]
    subprocess.run(command, check=True, timeout=60)


def reversed_records(lines):
    header = [line for line in lines if line.startswith("#")]
    return header + [line for line in lines if not line.startswith("#")][::-1]


def memory_owner(array):
    while isinstance(array.base, np.ndarray):
        array = array.base
    return array


def test_sparse_arrays_arrive_without_a_copy():
    adata = read_quietly(DATA / "basic.vcf")
    for matrix in (adata.X, *adata.layers.values()):
        for array in (matrix.data, matrix.indices, matrix.indptr):
            assert not memory_owner(array).flags.owndata


def test_contigs_come_from_the_header_in_its_order(tmp_path):
    text = (DATA / "no-csq.vcf").read_text()
    path = tmp_path / "reversed.vcf"
    path.write_text(
        text
        .replace("##contig=<ID=chr1,", "##contig=<ID=TEMP,")
        .replace("##contig=<ID=chr17,", "##contig=<ID=chr1,")
        .replace("##contig=<ID=TEMP,", "##contig=<ID=chr17,")
    )
    assert list(read_quietly(path).uns["contigs"]) == ["chr17", "chr1"]
    assert list(read_quietly(DATA / "no-records.vcf").uns["contigs"]) == ["chr1", "chr17"]


def test_contigs_fall_back_to_the_records_without_a_header_list(tmp_path):
    lines = (DATA / "no-csq.vcf").read_text().splitlines(keepends=True)
    records = [line for line in lines if not line.startswith("##contig")]
    path = tmp_path / "unlisted.vcf"
    path.write_text("".join(reversed_records(records)))
    assert list(read_quietly(path).uns["contigs"]) == ["chr17", "chr1"]


def test_the_call_is_the_provenance():
    path = DATA / "no-csq.vcf"
    adata = read_quietly(path, genes=DATA / "panel-genes.txt", threads=2)
    assert adata.uns["command_line"] == (
        f"annatar.read_vcf({str(path)!r}, genes={str(DATA / 'panel-genes.txt')!r}, threads=2)"
    )
    assert read_quietly(path).uns["command_line"] == f"annatar.read_vcf({str(path)!r})"


def test_the_version_is_the_package_version():
    adata = read_quietly(DATA / "basic.vcf")
    assert annatar.__version__ == version("annatar") == adata.uns["annatar_version"]


def test_only_the_api_is_public():
    assert sorted(annatar.__all__) == ["__version__", "concat", "read_vcf"]
