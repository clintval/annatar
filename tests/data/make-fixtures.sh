#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bgzip -c basic.vcf > basic.vcf.gz
bcftools view -Ob -o basic.bcf basic.vcf
bcftools view -Ob -o typed-fields.bcf typed-fields.vcf
