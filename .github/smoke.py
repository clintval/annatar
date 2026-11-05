"""Check that an installed annatar imports and converts and merges, as a library and a command."""

import subprocess
import sys
from pathlib import Path

import anndata
import h5py

import annatar

COMMAND = Path(sys.executable).with_name("annatar")
HEADER = [
    "##fileformat=VCFv4.3",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
]
COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "c1", "c2"]
RECORD = ["chr1", "100", ".", "A", "G", "50", "PASS", ".", "GT:AD", "0/1:5,3", "0/0:7,0"]

with open("smoke.vcf", "w") as sink:
    sink.write("\n".join([*HEADER, "\t".join(COLUMNS), "\t".join(RECORD)]) + "\n")

adata = annatar.read_vcf("smoke.vcf")
assert adata.shape == (2, 1), adata.shape
assert adata.X.dtype == "int8", adata.X.dtype

shown = subprocess.run([COMMAND, "--help"], check=True, capture_output=True, text=True)
assert "merge" in shown.stdout, shown.stdout

subprocess.run([COMMAND, "smoke.vcf", "smoke.h5ad"], check=True)
with h5py.File("smoke.h5ad") as h5ad:
    assert h5ad["X"].attrs["encoding-type"] == "csr_matrix"
assert anndata.read_h5ad("smoke.h5ad").uns["command_line"].endswith("smoke.vcf smoke.h5ad")

merged = annatar.concat({"a": "smoke.h5ad", "b": adata})
assert merged.var["alt_reads"].tolist() == [6], merged.var

subprocess.run([COMMAND, "merge", "--output", "merged.h5ad", "smoke.h5ad"], check=True)
assert anndata.read_h5ad("merged.h5ad").obs_names[0] == "c1-smoke"
print(f"annatar {annatar.__version__} works")
