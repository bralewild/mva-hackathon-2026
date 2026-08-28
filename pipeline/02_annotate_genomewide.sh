#!/usr/bin/env bash
# ==============================================================================
# 02_annotate_genomewide.sh - Blind genome-wide functional annotation
#
# INPUT  : $VCF_RAW  (5,012,204 variants)
# OUTPUT : $WORK/02_annotated.vcf.gz + .tbi
#          $RESULTS/02_snpeff_summary.html
#
# PURPOSE: annotate EVERY variant with no bias toward candidate genes. This is
#          the step that makes the search blind: the pipeline does not know the
#          disease is MVA, nor that BUB1B / CEP57 / TRIP13 exist.
#
# DATABASE CHOICE: GRCh38.115 (full Ensembl) rather than GRCh38.mane.*
#          MANE covers ~19.3k genes with one transcript each - ideal for
#          REPORTING, but a blind search needs complete coverage. The final
#          report does use MANE transcripts, via VEP in stage 04.
#
# LESSON LEARNED: snpEff exits 0 even when it fails (invalid option, dead
#          download). 'set -e' is not enough - we validate the output.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"

# --- dependency: annotation database ---
bash "$(dirname "$0")/02a_download_snpeff_db.sh"

OUT="$WORK/02_annotated.vcf.gz"
log "annotating 5,012,204 variants with snpEff $SNPEFF_DB (~15-60 min)"

snpEff -Xmx${JAVA_MEM} ann \
  -noLog \
  -stats "$RESULTS/02_snpeff_summary.html" \
  "$SNPEFF_DB" "$VCF_RAW" \
  2> "$LOGS/02_snpeff.log" \
  | bgzip -@ 8 > "$OUT"

# --- VALIDATION: never trust the exit code ---
SIZE=$(stat -c%s "$OUT")
if [ "$SIZE" -lt 10000000 ]; then
  log "ERROR: output is only ${SIZE} bytes - snpEff failed. Tail of the log:"
  tail -20 "$LOGS/02_snpeff.log"
  exit 1
fi
N=$(bcftools view -H "$OUT" | wc -l)
if [ "$N" -lt 5000000 ]; then
  log "ERROR: $N variants annotated, expected ~5,012,204"; exit 1
fi

tabix -f -p vcf "$OUT"
log "OK - $N variants annotated"
log "output  -> $OUT ($(du -h "$OUT" | cut -f1))"
log "summary -> $RESULTS/02_snpeff_summary.html"
