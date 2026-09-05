#!/usr/bin/env bash
# Phase 1 derivation: raw chr20 VCF -> truth/array/panel VCFs + site lists.
# Run from the repo root. Requires: bcftools >= 1.19, tabix, python venv active.
#
# Kinship note (BUILD-SPEC 1.1/2.2): the first run uses no removal list. After
# Phase 2.2 produces plink2's king-cutoff removals, re-run with:
#   REMOVALS=data/derived/king.removals.txt bash scripts/phase1.sh
set -euo pipefail

RAW=data/raw
DER=data/derived
VCF=$RAW/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
PANEL=$RAW/integrated_call_samples_v3.20130502.ALL.panel
REMOVALS=${REMOVALS:-}

# NOTE: the v5b chr20 release carries no rsIDs (ID='.' on every row; docs/verification.md).
# IDs are therefore set to chr:pos and array sites come from deterministic thinning
# (--mode thin). w_hm3.snplist stays in provenance; the hm3 code path is retained unused.

for f in "$VCF" "$PANEL"; do
  [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
mkdir -p "$DER"

echo "== 1/5 biallelic SNP subset + chr:pos IDs (phase-preserving) =="
bcftools view -m2 -M2 -v snps "$VCF" -Ou \
  | bcftools annotate --set-id '%CHROM:%POS' -Oz -o "$DER/chr20.bisnp.vcf.gz"
bcftools index -f -t "$DER/chr20.bisnp.vcf.gz"

echo "== 2/5 site table (ID, POS, AF) =="
bcftools query -f '%ID\t%POS\t%INFO/AF\n' "$DER/chr20.bisnp.vcf.gz" > "$DER/sites_all.tsv"
wc -l "$DER/sites_all.tsv"

echo "== 3/5 sample design =="
if [ -n "$REMOVALS" ]; then
  python -m src.cli design-samples --panel "$PANEL" --outdir "$DER" --removals "$REMOVALS"
else
  python -m src.cli design-samples --panel "$PANEL" --outdir "$DER"
fi

echo "== 4/5 site partition (thin mode) =="
python -m src.cli make-sites --mode thin --site-table "$DER/sites_all.tsv" --outdir "$DER"

echo "== 5/5 derived VCFs =="
bcftools view -S "$DER/targets_afr.txt" "$DER/chr20.bisnp.vcf.gz" \
  -Oz -o "$DER/truth_targets.vcf.gz"
bcftools view -S "$DER/targets_afr.txt" -i "ID=@$DER/array_sites.rsids.txt" \
  "$DER/chr20.bisnp.vcf.gz" -Oz -o "$DER/array_targets.vcf.gz"
bcftools view -S "$DER/panel_a_afr.txt" "$DER/chr20.bisnp.vcf.gz" \
  -Oz -o "$DER/panelA.vcf.gz"
bcftools view -S "$DER/panel_b_eur.txt" "$DER/chr20.bisnp.vcf.gz" \
  -Oz -o "$DER/panelB.vcf.gz"
for f in truth_targets array_targets panelA panelB; do
  bcftools index -f -t "$DER/$f.vcf.gz"
done

echo "== assertions: VCF sample headers exactly match the designed lists =="
diff <(bcftools query -l "$DER/truth_targets.vcf.gz") "$DER/targets_afr.txt"
diff <(bcftools query -l "$DER/array_targets.vcf.gz") "$DER/targets_afr.txt"
diff <(bcftools query -l "$DER/panelA.vcf.gz") "$DER/panel_a_afr.txt"
diff <(bcftools query -l "$DER/panelB.vcf.gz") "$DER/panel_b_eur.txt"

echo "== assertion: array_targets site count matches expectation =="
N_ARRAY_VCF=$(bcftools index -n "$DER/array_targets.vcf.gz")
echo "array_targets variants: $N_ARRAY_VCF (compare against n_array_sites in sites_summary.json)"

echo "PHASE 1 DERIVATION COMPLETE - record counts + md5sums in data/PROVENANCE.md"
