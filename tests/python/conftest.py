import os
import re
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

import annatar

REPO = Path(__file__).resolve().parents[2]
BINARY = Path(sys.executable).with_name("annatar")
DATA = REPO / "tests" / "data"


def logged(run: subprocess.CompletedProcess, level: str) -> list[str]:
    return re.findall(rf"^\[\S+ {level}\s+\S+\] (.*)$", run.stderr, flags=re.MULTILINE)


def read_quietly(path, **options):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return annatar.read_vcf(path, **options)


@pytest.fixture
def run_cli(tmp_path):
    environment = {**os.environ, "NO_COLOR": "1"}
    environment.pop("RUST_LOG", None)

    def run(path, *flags):
        command = [BINARY, *flags, path, tmp_path / "cli.h5ad"]
        return subprocess.run(command, capture_output=True, text=True, env=environment)

    return run


@pytest.fixture
def cli_output(tmp_path) -> Path:
    return tmp_path / "cli.h5ad"
