#!/usr/bin/env bash
# ==============================================================================
# 01_qc_baseline.sh - Characterise the raw VCF
#
# INPUT  : $VCF_RAW
# OUTPUT : $RESULTS/01_qc_baseline.txt
#
# PURPOSE: establish genome build, caller, contig naming, available fields and
#          baseline counts BEFORE filtering anything. Without this, no downstream
#          decision is auditable.
#
# NOTE: this script disables 'pipefail' on purpose. Patterns like
#       `... | grep | head` raise SIGPIPE in the producer, and with
#       pipefail + set -e the script would abort silently.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail

OUT="$RESULTS/01_qc_baseline.txt"
HDR="$WORK/01_header.txt"
bcftools view -h "$VCF_RAW" > "$HDR"

{
  echo "=============================================================="
  echo " QC BASELINE - $(date -Iseconds)"
  echo "=============================================================="
  echo
  echo "## Sample"
  bcftools query -l "$VCF_RAW" | sed 's/^/  /'
  echo
  echo "## Reference declared in the header"
  grep -E '^##reference' "$HDR" | sed 's/^/  /'
  echo
  echo "## Caller"
  grep -oE 'Version="[0-9.]+"' "$HDR" | sort -u | sed 's/^/  GATK /'
  echo "  GATKCommandLine lines: $(grep -c '^##GATKCommandLine' "$HDR")"
  echo
  echo "## Contig naming (CRITICAL for the submission)"
  grep '^##contig' "$HDR" | head -3 | sed 's/^/  /'
  echo "  total contigs: $(grep -c '^##contig' "$HDR")"
  if grep -q '^##contig=<ID=chr' "$HDR"; then
    echo "  STYLE: UCSC (with chr prefix)"
  else
    echo "  STYLE: Ensembl (NO chr prefix) -> the submission must prepend '$SUBMIT_CONTIG_PREFIX'"
  fi
  echo
  echo "## Available FORMAT fields"
  grep '^##FORMAT' "$HDR" | sed 's/,Description.*//;s/^/  /'
  echo
  echo "## Declared filters"
  grep '^##FILTER' "$HDR" | sed 's/^/  /'
  echo
  echo "## Counts"
  echo "  total : $(bcftools view -H "$VCF_RAW" | wc -l)"
  echo "  PASS  : $(bcftools view -H -f PASS "$VCF_RAW" | wc -l)"
} > "$OUT"

cat "$OUT"
log "written -> $OUT"
