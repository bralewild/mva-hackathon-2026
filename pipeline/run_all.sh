#!/usr/bin/env bash
# ==============================================================================
# run_all.sh — Ejecuta el pipeline completo de busqueda CIEGA, de punta a punta.
#
# Uso:
#   wsl -d Ubuntu-24.04 -- bash -lc "bash /mnt/c/.../pipeline/run_all.sh"
#
# Cada etapa es idempotente: si su salida ya existe y es valida, se saltea.
# ==============================================================================
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"
source "$D/00_config.sh"

log "===== 01 QC baseline ====="
bash "$D/01_qc_baseline.sh" > /dev/null

log "===== 00b terminos HPO del paciente ====="
python "$D/00b_extract_phenotype.py"

log "===== 02 anotacion genome-wide (CIEGA) ====="
if [ -f "$WORK/02_annotated.vcf.gz" ] && [ "$(stat -c%s "$WORK/02_annotated.vcf.gz")" -ge 10000000 ]; then
  log "ya anotado, se saltea"
else
  bash "$D/02_annotate_genomewide.sh"
fi

log "===== 03 modelos de herencia ====="
python "$D/03_inheritance_models.py"

log "===== 04 frecuencia poblacional y evidencia clinica ====="
python "$D/04_frequency_clinical.py"

log "===== 05 ranking fenotipico ====="
python "$D/05_phenotype_rank.py"

log "===== espejando resultados al proyecto ====="
bash "$D/sync_results.sh"

log "PIPELINE COMPLETO"
