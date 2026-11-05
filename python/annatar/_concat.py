from __future__ import annotations

from collections.abc import Iterable
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse

from annatar import _native
from annatar._vcf import StrPath
from annatar._vcf import _text_index

_NULLABLE = {"bool": "boolean", "int64": "Int64"}
_RECOUNTED = ("n_cells_covered", "n_cells_alt", "ref_reads", "alt_reads")
_IDENTITY = ("chrom", "pos", "ref", "alt")
_TEXT_IDENTITY = {"chrom": object, "ref": object, "alt": object}


def _named(items: Mapping[str, StrPath | AnnData] | Iterable[StrPath | AnnData]) -> dict[str, Any]:
    if isinstance(items, Mapping):
        return dict(items)
    named: dict[str, Any] = {}
    for position, item in enumerate(items):
        name = str(position) if isinstance(item, AnnData) else Path(item).stem
        if name in named:
            raise ValueError(
                f"more than one library is named {name}; "
                "rename an input or pass annatar.concat a mapping of names"
            )
        named[name] = item
    return named


def _read(path: StrPath) -> AnnData:
    if not Path(path).exists():
        raise FileNotFoundError(f"opening {path}: No such file or directory")
    try:
        return ad.read_h5ad(path)
    except Exception as error:
        raise ValueError(f"{path} is not annatar's h5ad: {error}") from error


def _is_sparse_integer(matrix: Any) -> bool:
    return issparse(matrix) and np.issubdtype(matrix.dtype, np.integer)


def _load(name: str, item: StrPath | AnnData) -> AnnData:
    library = item if isinstance(item, AnnData) else _read(item)
    layers: Any = library.layers
    counts = {"X": library.X} | {f"layer {n}": layers.get(n) for n in ("AD_ref", "AD_alt")}
    missing = [
        f"sparse integer {label}" for label, m in counts.items() if not _is_sparse_integer(m)
    ]
    columns = (*_IDENTITY, *_RECOUNTED)
    missing += [f"var column {name}" for name in columns if name not in library.var]
    if missing:
        source = name if isinstance(item, AnnData) else item
        raise ValueError(f"{source} is not annatar's h5ad: it has no {', '.join(missing)}")
    return library


def _contig_ranks(libraries: list[Any]) -> dict[str, int]:
    listed = [contig for library in libraries for contig in library.uns.get("contigs", ())]
    seen = [c for library in libraries for c in library.var["chrom"].astype(object).unique()]
    return {contig: rank for rank, contig in enumerate(dict.fromkeys([*listed, *seen]))}


def _variant_order(libraries: list[Any], ranks: dict[str, int]) -> pd.Index:
    frames = [library.var[list(_IDENTITY)].astype(_TEXT_IDENTITY) for library in libraries]
    filled = [frame for frame in frames if len(frame)]
    if not filled:
        return frames[0].index
    keys = pd.concat(filled)
    keys = keys[~keys.index.duplicated()]
    order = _native.variant_order(
        keys["chrom"].map(ranks).tolist(),
        keys["pos"].tolist(),
        keys["ref"].tolist(),
        keys["alt"].tolist(),
    )
    return keys.index[order]


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
    if column.dtype == "string" or isinstance(column.dtype, pd.CategoricalDtype):
        values = column.astype(object)
        categories = pd.Index(values.dropna().unique().tolist())
        return values.astype(pd.CategoricalDtype(categories))
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


def _recount(var: pd.DataFrame, layers: dict[str, Any]) -> None:
    ref, alt = layers["AD_ref"], layers["AD_alt"]
    var["n_cells_covered"] = _per_variant((ref + alt) > 0)
    var["n_cells_alt"] = _per_variant(alt > 0)
    var["ref_reads"] = _per_variant(ref)
    var["alt_reads"] = _per_variant(alt)


def _layer_dtypes(libraries: list[Any]) -> dict[str, Any]:
    dtypes: dict[str, Any] = {}
    for library in libraries:
        for name, layer in library.layers.items():
            if name is not None and (name not in dtypes or dtypes[name].kind != "f"):
                dtypes[name] = layer.dtype
    return dtypes


def _canonical(matrix: Any, dtype: Any) -> Any:
    matrix = matrix.astype(dtype, copy=False)
    matrix.sort_indices()
    return matrix


def concat(items: Mapping[str, StrPath | AnnData] | Iterable[StrPath | AnnData]) -> AnnData:
    """Merge libraries into one AnnData, keeping every variant's annotation.

    Cells keep their barcodes with the library name appended and name their library in
    `obs["library"]`, so merging a merged file appends a library name again. An `obs` column some
    libraries lack is missing for their cells. Variants are the union of the libraries, sorted by
    contig in the order of `uns["contigs"]`, then `pos`, `ref` and `alt`. A contig on no list goes
    after the listed ones. The columns `n_cells_covered`, `n_cells_alt`, `ref_reads` and
    `alt_reads` are recounted from the merged layers. Every other `var` value on a shared variant,
    such as `qual` or an INFO field, comes from the first library holding it, and later libraries
    fill only the columns that library lacks entirely.

    Args:
        items: Paths to h5ad files or AnnData objects. A mapping names each library, otherwise a
            file is named after its path and an AnnData after its position.

    Returns:
        The merged AnnData.

    Raises:
        ValueError: If there is nothing to merge, two libraries have one name, or a library is not
            annatar's.
        OSError: If a file cannot be read.
    """
    named = {name: _load(name, item) for name, item in _named(items).items()}
    merged: Any = ad.concat(
        named, join="outer", label="library", index_unique="-", uns_merge="same"
    )
    libraries: list[Any] = list(named.values())
    ranks = _contig_ranks(libraries)
    order = _variant_order(libraries, ranks)
    if not merged.var_names.equals(order):
        merged = merged[:, order].copy()
    merged_layers: Any = merged.layers
    layers = {
        name: _canonical(merged_layers[name], dtype)
        for name, dtype in sorted(_layer_dtypes(libraries).items())
    }
    var = _merged_var(libraries, merged.var_names)
    _recount(var, layers)
    return AnnData(
        X=_canonical(merged.X, libraries[0].X.dtype),
        obs=merged.obs.set_axis(_text_index(merged.obs_names)),
        var=var.set_axis(_text_index(var.index)),
        layers=layers,
        uns=dict(sorted({**merged.uns, "contigs": np.array(list(ranks), dtype=object)}.items())),
    )
