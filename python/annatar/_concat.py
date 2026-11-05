from __future__ import annotations

from collections.abc import Iterable
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from annatar._vcf import StrPath

_NULLABLE = {"bool": "boolean", "int64": "Int64"}


def _named(items: Mapping[str, StrPath | AnnData] | Iterable[StrPath | AnnData]) -> dict[str, Any]:
    if isinstance(items, Mapping):
        return dict(items)
    named: dict[str, Any] = {}
    for position, item in enumerate(items):
        name = str(position) if isinstance(item, AnnData) else Path(item).stem
        if name in named:
            raise ValueError(f"more than one library is named {name}; pass a mapping of names")
        named[name] = item
    return named


def _stacked(parts: list[pd.Series]) -> pd.Series:
    filled = [part for part in parts if len(part)] or parts[:1]
    stacked = pd.concat(filled)
    if stacked.dtype == object:
        return pd.concat([part.astype("string") for part in filled])
    return stacked.astype("float64") if stacked.dtype == "Float64" else stacked


def _merged_column(parts: list[pd.Series], var_names: pd.Index) -> pd.Series:
    stacked = _stacked(parts)
    column = stacked[~stacked.index.duplicated()]
    if len(column) < len(var_names):
        column = column.astype(_NULLABLE.get(column.dtype.name, column.dtype))
    column = column.reindex(var_names)
    if column.dtype == "string":
        categories = pd.Index(column.dropna().unique().tolist())
        return column.astype(pd.CategoricalDtype(categories))
    return column


def _merged_var(libraries: list[AnnData], var_names: pd.Index) -> pd.DataFrame:
    names = dict.fromkeys(name for library in libraries for name in library.var.columns)
    columns = {
        name: _merged_column(
            [library.var[name] for library in libraries if name in library.var], var_names
        )
        for name in names
    }
    return pd.DataFrame(columns, index=var_names)


def _per_variant(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=0)).ravel().astype("int64")


def _recount(merged: AnnData) -> None:
    layers: Any = merged.layers
    ref, alt = layers["AD_ref"], layers["AD_alt"]
    merged.var["n_cells_covered"] = _per_variant((ref + alt) > 0)
    merged.var["n_cells_alt"] = _per_variant(alt > 0)
    merged.var["ref_reads"] = _per_variant(ref)
    merged.var["alt_reads"] = _per_variant(alt)


def _restore_dtypes(merged: Any, libraries: list[Any]) -> None:
    merged.X = merged.X.astype(libraries[0].X.dtype, copy=False)
    dtypes = {name: layer.dtype for library in libraries for name, layer in library.layers.items()}
    for name, dtype in dtypes.items():
        merged.layers[name] = merged.layers[name].astype(dtype, copy=False)


def concat(items: Mapping[str, StrPath | AnnData] | Iterable[StrPath | AnnData]) -> AnnData:
    """Merge libraries into one AnnData, keeping every variant's annotation.

    Cells keep their barcodes with the library name appended and name their library in
    `obs["library"]`. Variants are the union of the libraries. The columns `n_cells_covered`,
    `n_cells_alt`, `ref_reads` and `alt_reads` are recounted from the merged layers. Every other
    `var` value on a shared variant, such as `qual` or an INFO field, comes from the first library
    holding it, and later libraries fill only the columns that library lacks entirely.

    Args:
        items: Paths to h5ad files or AnnData objects. A mapping names each library, otherwise a
            file is named after its path and an AnnData after its position.

    Returns:
        The merged AnnData.

    Raises:
        ValueError: If there is nothing to merge or two libraries have one name.
    """
    named = {
        name: item if isinstance(item, AnnData) else ad.read_h5ad(item)
        for name, item in _named(items).items()
    }
    merged = ad.concat(named, join="outer", label="library", index_unique="-", uns_merge="same")
    libraries = list(named.values())
    _restore_dtypes(merged, libraries)
    merged.var = _merged_var(libraries, merged.var_names)
    _recount(merged)
    return merged
