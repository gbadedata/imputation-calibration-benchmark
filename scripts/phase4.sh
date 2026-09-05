#!/usr/bin/env bash
# Phase 4: extract DS/truth/DR2 tables and run the calibration analysis.
# Run from repo root AFTER phase3.sh gates.
set -euo pipefail

DER=data/derived
RES=results
EX=$DER/qc/phase4
EVAL=$DER/evaluation_sites.rsids.txt
TRUTH=$DER/truth_targets.vcf.gz

for f in "$TRUTH" "$RES/imputed_A.vcf.gz" "$RES/imputed_B.vcf.gz" "$EVAL"; do
  [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
mkdir -p "$EX" evidence/figures

echo "== 4.0 extract truth genotypes at evaluation sites =="
bcftools query -l "$TRUTH" > "$EX/truth_samples.txt"
bcftools query -i "ID=@$EVAL" -f '%ID[\t%GT]\n' "$TRUTH" > "$EX/truth_gt.tsv"
wc -l "$EX/truth_gt.tsv"

for C in A B; do
  echo "== 4.0 extract Condition $C dosages + DR2 at evaluation sites =="
  bcftools query -l "$RES/imputed_$C.vcf.gz" > "$EX/${C}_samples.txt"
  bcftools query -i "ID=@$EVAL" -f '%ID[\t%DS]\n' "$RES/imputed_$C.vcf.gz" > "$EX/${C}_ds.tsv"
  cp "$DER/qc/phase3/${C}_dr2.tsv" "$EX/${C}_dr2.tsv"
done

echo "== 4.1-4.5 accuracy, calibration, abstention, figures, report =="
python -m src.cli calibrate \
  --extract-dir "$EX" \
  --report "$RES/benchmark_report.json" \
  --figures-dir evidence/figures

echo "PHASE 4 COMPLETE - the headline is in results/benchmark_report.json and evidence/figures/"
