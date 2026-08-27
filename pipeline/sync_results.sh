#!/usr/bin/env bash
# ==============================================================================
# sync_results.sh — Espeja resultados LIVIANOS de WSL al proyecto en Windows
#
# ENTRADA : $RESULTS/*  (en ext4)
# SALIDA  : $PROJECT/results/*  (visible en el proyecto, versionable)
#
# PROPOSITO: los datos pesados se quedan en ext4 por rendimiento, pero los
#            reportes (texto, TSV, CSV, HTML) son de KB y deben verse en el
#            proyecto. Solo copia archivos < 5 MB y nunca VCF/BAM/FASTQ.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail
DST="$PROJECT/results"
mkdir -p "$DST"

n=0
for f in "$RESULTS"/*; do
  [ -f "$f" ] || continue
  case "$f" in *.vcf|*.vcf.gz|*.bam|*.cram|*.fastq*|*.tbi) continue;; esac
  sz=$(stat -c%s "$f")
  [ "$sz" -gt 5242880 ] && { log "omitido (>5MB): $(basename "$f")"; continue; }
  cp -f "$f" "$DST/" && n=$((n+1))
done
log "espejados $n archivos -> $DST"
ls -lh "$DST" | awk 'NR>1 {printf "  %8s  %s\n", $5, $9}'
