#!/usr/bin/env bash
# ==============================================================================
# run_all.sh - Run the complete blind-search pipeline, end to end.
#
# Usage:
#   wsl -d Ubuntu-24.04 -- bash -lc "bash /mnt/c/.../pipeline/run_all.sh"
#
# Every stage is idempotent: if its output already exists and is valid, it is
# skipped. Safe to re-run.
# ==============================================================================
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"
source "$D/00_config.sh"

log "===== 00 preflight self-check ====="
bash "$D/selfcheck.sh" > /dev/null || { log "preflight failed - run: bash pipeline/selfcheck.sh"; exit 1; }

log "===== 01 QC baseline ====="
bash "$D/01_qc_baseline.sh" > /dev/null

log "===== 00b patient HPO terms ====="
python "$D/00b_extract_phenotype.py"

log "===== 02 genome-wide annotation (BLIND) ====="
if [ -f "$WORK/02_annotated.vcf.gz" ] && [ "$(stat -c%s "$WORK/02_annotated.vcf.gz")" -ge 10000000 ]; then
  log "already annotated, skipping"
else
  bash "$D/02_annotate_genomewide.sh"
fi

log "===== 03 inheritance models ====="
python "$D/03_inheritance_models.py"

log "===== 04 population frequency and clinical evidence ====="
python "$D/04_frequency_clinical.py"

log "===== 05 phenotype ranking ====="
python "$D/05_phenotype_rank.py"

log "===== 06 convergence gate ====="
python "$D/06_validate_convergence.py"

log "===== 07 secondary findings ====="
python "$D/07_secondary_findings.py"

log "===== 08 mosaic aneuploidy screen ====="
python "$D/08_mosaic_aneuploidy.py" > /dev/null

log "===== 09 submission file (automated) ====="
python "$D/09_make_submission.py"

log "===== mirroring results to the project folder ====="
bash "$D/sync_results.sh"

log "PIPELINE COMPLETE"
