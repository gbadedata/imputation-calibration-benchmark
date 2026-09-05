#!/usr/bin/env bash
# Phase 2: kinship removal -> redraw -> chip-cohort QC -> QC report -> BGEN evidence.
# Run from the repo root AFTER a first phase1.sh run. Requires: bcftools, plink2, python venv.
#
# Ordering rationale (BUILD-SPEC 2.1): KING runs across ALL AFR+EUR samples so the
# post-removal redraw is covered; kinship precedes variant QC because it changes the
# target set. Chip QC then runs on the FINAL 100 targets only (BUILD-SPEC 2.2).
set -euo pipefail

DER=data/derived
QC=data/derived/qc
RES=results
mkdir -p "$QC" "$RES"

echo "== 2.1a AFR+EUR candidate list (all of both superpops, psam-style header) =="
awk 'NR>1 && ($3=="AFR" || $3=="EUR"){print $1}' \
  data/raw/integrated_call_samples_v3.20130502.ALL.panel > "$QC/afr_eur.iids.tmp"
{ echo "#IID"; cat "$QC/afr_eur.iids.tmp"; } > "$QC/afr_eur.keep.txt"
rm "$QC/afr_eur.iids.tmp"
N_CAND=$(( $(wc -l < "$QC/afr_eur.keep.txt") - 1 ))
echo "candidates: $N_CAND (expect 1164 = 661 AFR + 503 EUR)"

echo "== 2.1b pgen at array sites for KING =="
plink2 --vcf "$DER/chr20.bisnp.vcf.gz" \
  --extract "$DER/array_sites.rsids.txt" \
  --keep "$QC/afr_eur.keep.txt" \
  --rm-dup force-first \
  --make-pgen --out "$QC/array_afr_eur"

echo "== 2.1c KING kinship (cutoff 0.0884, 2nd degree) =="
plink2 --pfile "$QC/array_afr_eur" --king-cutoff 0.0884 --out "$QC/king"
grep -v '^#' "$QC/king.king.cutoff.out.id" | awk '{print $NF}' > "$QC/king.removals.txt" || true
N_REMOVED=$(wc -l < "$QC/king.removals.txt")
echo "kinship removals: $N_REMOVED"

echo "== 2.1d redraw + rebuild derived data with removals (the one permitted loop) =="
REMOVALS="$QC/king.removals.txt" bash scripts/phase1.sh

echo "== 2.2 chip-cohort variant QC on the final 100 targets =="
plink2 --vcf "$DER/array_targets.vcf.gz" \
  --geno 0.02 --maf 0.01 --hwe 1e-6 \
  --write-snplist --out "$QC/target_qc"
N_KEPT=$(wc -l < "$QC/target_qc.snplist")
echo "array sites kept after chip QC: $N_KEPT"

echo "== 2.2b QC-passed target chip data =="
bcftools view -i "ID=@$QC/target_qc.snplist" "$DER/array_targets.vcf.gz" \
  -Oz -o "$DER/array_targets.qc.vcf.gz"
bcftools index -f -t "$DER/array_targets.qc.vcf.gz"
N_VCF=$(bcftools index -n "$DER/array_targets.qc.vcf.gz")
[ "$N_VCF" -eq "$N_KEPT" ] || { echo "FATAL: qc VCF has $N_VCF variants, snplist $N_KEPT" >&2; exit 1; }

echo "== 2.3 QC report (accounting identity enforced inside) =="
python -m src.cli qc-report \
  --qc-log "$QC/target_qc.log" \
  --snplist "$QC/target_qc.snplist" \
  --king-id "$QC/king.king.cutoff.out.id" \
  --outdir "$RES"

echo "== 2.4 BGEN 1.2 round trip (ref-last: docs/verification.md V3) =="
plink2 --vcf "$DER/array_targets.qc.vcf.gz" \
  --export bgen-1.2 bits=8 --out "$QC/targets_bgen"
plink2 --bgen "$QC/targets_bgen.bgen" ref-last --sample "$QC/targets_bgen.sample" \
  --export vcf bgz --out "$QC/targets_roundtrip"
bcftools query -f '%POS[\t%GT]\n' "$DER/array_targets.qc.vcf.gz" > "$QC/sig_original.txt"
bcftools query -f '%POS[\t%GT]\n' "$QC/targets_roundtrip.vcf.gz" > "$QC/sig_roundtrip.txt"
python -m src.cli bgen-verify \
  --original "$QC/sig_original.txt" \
  --round-tripped "$QC/sig_roundtrip.txt"

echo "PHASE 2 COMPLETE - commit results/qc_report.{json,md}, append derived md5s to PROVENANCE"
