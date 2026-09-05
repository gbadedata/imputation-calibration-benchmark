#!/usr/bin/env bash
# Build the tiny fixtures the Nextflow -profile test run uses:
# 10 target samples, 20 panel samples each, region 20:1-2,000,000.
set -euo pipefail
DER=data/derived
OUT=$DER/qc/nf_test
mkdir -p "$OUT"
head -10 "$DER/targets_afr.txt" > "$OUT/t10.txt"
head -20 "$DER/panel_a_afr.txt" > "$OUT/pa20.txt"
head -20 "$DER/panel_b_eur.txt" > "$OUT/pb20.txt"
bcftools view -r 20:1-2000000 -S "$OUT/t10.txt" "$DER/array_targets.qc.vcf.gz" -Oz -o "$OUT/gt_test.vcf.gz"
bcftools view -r 20:1-2000000 -S "$OUT/pa20.txt" "$DER/panelA.vcf.gz" -Oz -o "$OUT/panelA_test.vcf.gz"
bcftools view -r 20:1-2000000 -S "$OUT/pb20.txt" "$DER/panelB.vcf.gz" -Oz -o "$OUT/panelB_test.vcf.gz"
for f in gt_test panelA_test panelB_test; do bcftools index -f -t "$OUT/$f.vcf.gz"; done
echo "nf test fixtures ready in $OUT - run: (cd workflows && nextflow run main.nf -profile test)"
