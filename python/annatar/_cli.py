from __future__ import annotations

import contextlib
import os
import sys
import warnings
from pathlib import Path
from typing import Any

from annatar import _native


def _write(adata: Any, output: Path) -> tuple[int, int]:
    partial = output.with_name(f"{output.name}.partial")
    try:
        if not output.parent.is_dir():
            raise FileNotFoundError(f"{output.parent} is not a directory")
        adata.write_h5ad(partial, compression="gzip", compression_opts=4)
        partial.replace(output)
    except Exception as error:
        raise OSError(f"cannot write {output}: {error}") from error
    finally:
        with contextlib.suppress(OSError):
            partial.unlink()
    return adata.n_obs, adata.n_vars


def _write_conversion(conversion: _native.Conversion, output: Path) -> tuple[int, int]:
    from annatar._vcf import _to_anndata

    return _write(_to_anndata(conversion), output)


def _refuse_an_input(inputs: list[Path], output: Path) -> None:
    for target in (output, output.with_name(f"{output.name}.partial")):
        for path in inputs:
            if target.exists() and path.exists() and os.path.samefile(target, path):
                raise ValueError(f"{target} is also an input")


def _merge(inputs: list[Path], output: Path) -> tuple[int, int]:
    _refuse_an_input(inputs, output)
    from annatar._concat import concat

    return _write(concat(inputs), output)


def main() -> None:
    """Run the `annatar` command."""
    warnings.simplefilter("ignore")
    try:
        code = _native.main(sys.argv, _write_conversion, _merge)
    except KeyboardInterrupt:
        code = 130
    sys.exit(code)
