"""Convert a multi-sample VCF of single cells into AnnData."""

from annatar._concat import concat
from annatar._native import __version__
from annatar._vcf import read_vcf

__all__ = ["__version__", "concat", "read_vcf"]
