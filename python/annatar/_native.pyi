from collections.abc import Callable
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from typing import Any
from typing import TypedDict

import numpy as np
import numpy.typing as npt

Sparse = tuple[npt.NDArray[Any], npt.NDArray[np.integer[Any]], npt.NDArray[np.integer[Any]]]
NamedColumn = tuple[str, str, npt.NDArray[Any], Any]

class Conversion(TypedDict):
    obs_names: list[str]
    obs: list[NamedColumn]
    var_names: list[str]
    var: list[NamedColumn]
    x: Sparse
    layers: list[tuple[str, Sparse]]
    uns: list[tuple[str, str | list[str]]]

__version__: str

def variant_order(
    ranks: list[int],
    starts: list[int],
    ends: list[int],
    references: list[str],
    alternates: list[str],
) -> npt.NDArray[np.uintp]:
    """Return the positions of variants in order of contig rank, start, end, REF and ALT."""

def main(
    args: list[str],
    refuse_an_input: Callable[[list[Path], Path], None],
    write_conversion: Callable[[Conversion, Path], tuple[int, int]],
    merge: Callable[[list[Path], Path], tuple[int, int]],
) -> int:
    """Parse the command line and run it, returning its exit code."""

def read_vcf(
    path: str | PathLike[str],
    *,
    genes: str | PathLike[str] | None,
    mane_table: str | PathLike[str] | None,
    exclude_fields: Sequence[str],
    threads: int,
    command_line: str,
) -> Conversion:
    """Convert a VCF or BCF, handing the result over as NumPy arrays and Python lists."""
