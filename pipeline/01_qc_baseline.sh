#!/usr/bin/env bash
# ==============================================================================
# 01_qc_baseline.sh — Caracterizacion del VCF crudo
#
# ENTRADA : $VCF_RAW
# SALIDA  : $RESULTS/01_qc_baseline.txt
# PROPOSITO: establecer build, caller, naming de contigs, campos disponibles y
#            metricas base ANTES de filtrar nada. Sin esto no se puede auditar
#            ninguna decision posterior.
#
# NOTA: este script desactiva 'pipefail' a proposito. Los patrones
#       `... | grep | head` provocan SIGPIPE en el productor y con pipefail+set -e
#       abortarian el script silenciosamente.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail

OUT="$RESULTS/01_qc_baseline.txt"
HDR="$WORK/01_header.txt"
bcftools view -h "$VCF_RAW" > "$HDR"

{
  echo "=============================================================="
  echo " QC BASELINE — $(date -Iseconds)"
  echo "=============================================================="
  echo
  echo "## Muestra"
  bcftools query -l "$VCF_RAW" | sed 's/^/  /'
  echo
  echo "## Referencia declarada en el header"
  grep -E '^##reference' "$HDR" | sed 's/^/  /'
  echo
  echo "## Caller"
  grep -oE 'Version="[0-9.]+"' "$HDR" | sort -u | sed 's/^/  GATK /'
  echo "  lineas GATKCommandLine: $(grep -c '^##GATKCommandLine' "$HDR")"
  echo
  echo "## Naming de contigs (CRITICO para el submission)"
  grep '^##contig' "$HDR" | head -3 | sed 's/^/  /'
  echo "  total contigs: $(grep -c '^##contig' "$HDR")"
  if grep -q '^##contig=<ID=chr' "$HDR"; then
    echo "  ESTILO: UCSC (con prefijo chr)"
  else
    echo "  ESTILO: Ensembl (SIN prefijo chr) -> el submission necesita anteponer '$SUBMIT_CONTIG_PREFIX'"
  fi
  echo
  echo "## Campos FORMAT disponibles"
  grep '^##FORMAT' "$HDR" | sed 's/,Description.*//;s/^/  /'
  echo
  echo "## Filtros definidos"
  grep '^##FILTER' "$HDR" | sed 's/^/  /'
  echo
  echo "## Conteos"
  echo "  totales : $(bcftools view -H "$VCF_RAW" | wc -l)"
  echo "  PASS    : $(bcftools view -H -f PASS "$VCF_RAW" | wc -l)"
} > "$OUT"

cat "$OUT"
log "escrito -> $OUT"
