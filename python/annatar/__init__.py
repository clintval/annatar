"""Convert a multi-sample VCF of single cells into AnnData."""

from importlib import import_module
from typing import TYPE_CHECKING
from typing import Any

from annatar._native import __version__

if TYPE_CHECKING:
    from annatar._concat import concat
    from annatar._vcf import read_vcf

__all__ = ["__version__", "concat", "read_vcf"]

# The command imports this package, so anndata waits until something needs it.
_LAZY = {"concat": "annatar._concat", "read_vcf": "annatar._vcf"}


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        return getattr(import_module(_LAZY[name]), name)
    raise AttributeError(f"module 'annatar' has no attribute {name!r}")
