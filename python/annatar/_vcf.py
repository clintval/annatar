from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from annatar import _native

StrPath = str | os.PathLike[str]


def _array(kind: str, values: npt.NDArray[Any], extra: Any) -> Any:
    if kind == "Int64":
        return pd.arrays.IntegerArray(values, extra)
    if kind == "boolean":
        return pd.arrays.BooleanArray(values, extra)
    if kind == "category":
        return pd.Categorical.from_codes(values, categories=pd.Index(extra))
    return values


# Object strings write as string-array, which anndata 0.10 reads; pandas 3's str dtype does not.
def _text_index(names: Any) -> pd.Index:
    return pd.Index(names, dtype=object, name="index")


def _frame(names: list[str], columns: list[_native.NamedColumn]) -> pd.DataFrame:
    data = {name: _array(kind, values, extra) for name, kind, values, extra in columns}
    return pd.DataFrame(data, index=_text_index(names))


def _matrix(arrays: _native.Sparse, shape: tuple[int, int]) -> csr_matrix:
    return csr_matrix(arrays, shape=shape, copy=False)


def _to_anndata(conversion: _native.Conversion) -> AnnData:
    shape = len(conversion["obs_names"]), len(conversion["var_names"])
    uns = {
        name: np.array(value, dtype=object) if isinstance(value, list) else value
        for name, value in sorted(conversion["uns"])
    }
    return AnnData(
        X=_matrix(conversion["x"], shape),
        obs=_frame(conversion["obs_names"], conversion["obs"]),
        var=_frame(conversion["var_names"], conversion["var"]),
        layers={name: _matrix(arrays, shape) for name, arrays in sorted(conversion["layers"])},
        uns=uns,
    )


def _call(path: StrPath, **options: Any) -> str:
    arguments = [repr(os.fspath(path))]
    arguments += [f"{name}={value!r}" for name, value in options.items() if value]
    return f"annatar.read_vcf({', '.join(arguments)})"


def read_vcf(
    path: StrPath,
    *,
    genes: StrPath | None = None,
    mane_table: StrPath | None = None,
    exclude_fields: Sequence[str] = (),
    threads: int = 1,
) -> AnnData:
    """Convert a VCF or BCF of single cells into an in-memory AnnData.

    Args:
        path: The VCF or BCF, plain or bgzipped.
        genes: Panel genes preferred by the pick, one HGNC symbol or Ensembl gene ID per line.
        mane_table: NCBI MANE summary, plain or gzipped, used when CSQ has no MANE subfields.
        exclude_fields: Fields to leave out, named `INFO/ID`, `FORMAT/ID` or `INFO/CSQ/subfield`.
        threads: BGZF decompression threads.

    Returns:
        The cells by variants AnnData, with the same contents as the h5ad of the `annatar` command.

    Raises:
        ValueError: If the file is malformed or an option is invalid.
        OSError: If a file cannot be read.
    """
    conversion = _native.read_vcf(
        path,
        genes=genes,
        mane_table=mane_table,
        exclude_fields=exclude_fields,
        threads=threads,
        command_line=_call(
            path,
            genes=os.fspath(genes) if genes is not None else None,
            mane_table=os.fspath(mane_table) if mane_table is not None else None,
            exclude_fields=exclude_fields,
            threads=threads if threads != 1 else None,
        ),
    )
    return _to_anndata(conversion)
