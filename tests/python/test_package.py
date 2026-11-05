import subprocess
import sys
from pathlib import Path

import pytest
from conftest import BINARY
from conftest import DATA

import annatar._native as native

EXTENSION = Path(native.__file__)


def symbols(path: Path) -> list[str]:
    table = subprocess.run(["nm", "-a", path], capture_output=True, text=True, check=True)
    return [line.split()[-1].lstrip("_") for line in table.stdout.splitlines() if line.strip()]


def libraries(path: Path) -> str:
    tool = ["otool", "-L"] if sys.platform == "darwin" else ["ldd"]
    return subprocess.run([*tool, path], capture_output=True, text=True, check=True).stdout


def test_the_cli_shows_what_hdf5_looks_like():
    assert any(name.startswith("H5") for name in symbols(BINARY))
    assert b"HDF5" in BINARY.read_bytes()


def test_the_extension_has_no_hdf5_symbols():
    assert not [name for name in symbols(EXTENSION) if name.startswith("H5")]


def test_the_extension_holds_no_hdf5_library():
    assert b"HDF5" not in EXTENSION.read_bytes()


def test_the_extension_links_no_hdf5_library():
    assert "hdf5" not in libraries(EXTENSION).lower()


@pytest.mark.parametrize(
    "modules",
    [
        ("annatar", "h5py", "anndata"),
        ("annatar", "anndata", "h5py"),
        ("h5py", "annatar", "anndata"),
        ("h5py", "anndata", "annatar"),
        ("anndata", "annatar", "h5py"),
        ("anndata", "h5py", "annatar"),
    ],
)
def test_imports_work_in_any_order(modules, tmp_path):
    code = "\n".join([
        *(f"import {module}" for module in modules),
        "import sys",
        "adata = annatar.read_vcf(sys.argv[1])",
        "adata.write_h5ad(sys.argv[2])",
        "assert anndata.read_h5ad(sys.argv[2]).n_vars == adata.n_vars == 2",
        "assert h5py.File(sys.argv[2])['X'].attrs['encoding-type'] == 'csr_matrix'",
    ])
    command = [sys.executable, "-c", code, DATA / "no-csq.vcf", tmp_path / "out.h5ad"]
    subprocess.run(command, check=True, timeout=120)
