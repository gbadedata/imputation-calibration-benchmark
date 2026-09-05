#!/usr/bin/env bash
# Phase 3: dual-condition Beagle imputation. Run from repo root AFTER phase2.sh gates.
# Flow: full-sample regional dry run (memory gate) -> Condition A (matched AFR panel)
#       -> Condition B (mismatched EUR panel) -> per-condition validation -> run metrics.
set -euo pipefail

DER=data/derived
RAW=data/raw
RES=results
QC=$DER/qc
GT=$DER/array_targets.qc.vcf.gz
PANEL_A=$DER/panelA.vcf.gz
PANEL_B=$DER/panelB.vcf.gz
MAP=$RAW/plink.chr20.GRCh37.map
JAR=$RAW/beagle.27Feb25.75f.jar
EVAL=$DER/evaluation_sites.rsids.txt
TARGETS=$DER/targets_afr.txt
XMX=${ICB_BEAGLE_XMX_GB:-12}
NT=${ICB_BEAGLE_NTHREADS:-6}

for f in "$GT" "$PANEL_A" "$PANEL_B" "$MAP" "$JAR" "$EVAL" "$TARGETS"; do
  [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
mkdir -p "$RES" "$QC/phase3"

echo "== 3.0 memory dry-run gate (all samples, region 20:1-5,000,000) =="
bcftools view -r 20:1-5000000 "$GT" -Oz -o "$QC/phase3/gt_region.vcf.gz"
bcftools view -r 20:1-5000000 "$PANEL_A" -Oz -o "$QC/phase3/panelA_region.vcf.gz"
/usr/bin/time -v java -Xmx${XMX}g -jar "$JAR" \
  gt="$QC/phase3/gt_region.vcf.gz" ref="$QC/phase3/panelA_region.vcf.gz" \
  map="$MAP" out="$QC/phase3/dryrun" nthreads=$NT 2> "$QC/phase3/dryrun.time.log"
DRY_RSS=$(grep "Maximum resident set size" "$QC/phase3/dryrun.time.log" | awk '{print $NF}')
echo "dry-run peak RSS: ${DRY_RSS} kB (gate: must be well under $((XMX*1024*1024)) kB)"

run_condition () {
  local NAME=$1 REF=$2 COND=$3
  echo "== 3.1 Condition $NAME ($COND panel) =="
  /usr/bin/time -v java -Xmx${XMX}g -jar "$JAR" \
    gt="$GT" ref="$REF" map="$MAP" out="$RES/imputed_$NAME" nthreads=$NT \
    2> "$RES/imputed_$NAME.time.log"
  bcftools index -f -t "$RES/imputed_$NAME.vcf.gz"

  echo "== 3.2 extract validation inputs ($NAME) =="
  bcftools view -h "$RES/imputed_$NAME.vcf.gz" > "$QC/phase3/${NAME}_header.txt"
  bcftools query -l "$RES/imputed_$NAME.vcf.gz" > "$QC/phase3/${NAME}_samples.txt"
  bcftools query -f '%ID\n' "$RES/imputed_$NAME.vcf.gz" > "$QC/phase3/${NAME}_ids.txt"
  bcftools query -i "ID=@$EVAL" -f '%ID\t%INFO/DR2\n' "$RES/imputed_$NAME.vcf.gz" \
    > "$QC/phase3/${NAME}_dr2.tsv"

  echo "== 3.3 validate ($NAME, condition=$COND) =="
  python -m src.cli impute-check --condition "$COND" \
    --header "$QC/phase3/${NAME}_header.txt" \
    --samples "$QC/phase3/${NAME}_samples.txt" \
    --targets "$TARGETS" \
    --marker-ids "$QC/phase3/${NAME}_ids.txt" \
    --eval-ids "$EVAL" \
    --dr2-table "$QC/phase3/${NAME}_dr2.tsv" \
    --out "$RES/imputation_$NAME.json"
}

run_condition A "$PANEL_A" matched
run_condition B "$PANEL_B" mismatched

echo "== 3.4 run metrics =="
python - << 'PYEOF'
import json, re
from pathlib import Path

def timelog(path):
    t = Path(path).read_text()
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))
    wall = re.search(r"Elapsed \(wall clock\) time.*: ([\d:.]+)", t).group(1)
    return {"peak_rss_kb": rss, "wall_clock": wall}

metrics = {
    "dry_run": timelog("data/derived/qc/phase3/dryrun.time.log"),
    "condition_A": {**timelog("results/imputed_A.time.log"),
                    **json.loads(Path("results/imputation_A.json").read_text())},
    "condition_B": {**timelog("results/imputed_B.time.log"),
                    **json.loads(Path("results/imputation_B.json").read_text())},
}
Path("results/run_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in
    ("peak_rss_kb", "wall_clock", "eval_coverage", "n_markers_output", "dr2_missing_fraction")}
    for k, v in metrics.items()}, indent=2))
PYEOF

echo "PHASE 3 COMPLETE - commit results/imputation_{A,B}.json + run_metrics.json, append to PROVENANCE"
