#!/usr/bin/env bash
# ==============================================================================
# 02a_download_snpeff_db.sh - Resumable snpEff database download
#
# OUTPUT: ${SNPEFF_DIR}data/${SNPEFF_DB}/
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# snpEff's built-in downloader (`snpEff download`) implements no retry and no
# timeout: when the socket drops, the Java process sleeps forever with no error
# and no non-zero exit code. It happened here - frozen at 285 MB of 770, for
# 32 minutes, at 0.5 % CPU.
#
# This script uses curl with:
#   -C -                 resume from what is already downloaded
#   --retry 10           retry on transient failures
#   --speed-limit/-time  abort (and retry) if throughput drops below 10 KB/s
#                        for 60 s
# ==============================================================================
source "$(dirname "$0")/00_config.sh"

SNPEFF_DIR=$(ls -d "$HOME"/micromamba/envs/bio/share/snpeff-*/ | head -1)
URL="https://snpeff-public.s3.amazonaws.com/databases/v5_4/snpEff_v5_4_${SNPEFF_DB}.zip"
ZIP="/tmp/snpEff_v5_4_${SNPEFF_DB}.zip"

if [ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ]; then
  log "database already installed: $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"; exit 0
fi

# kill any hung snpEff download
pkill -f 'snpEff.jar download' 2>/dev/null && log "killed a hung snpEff download"

TOTAL=$(curl -sIL "$URL" | tr -d '\r' | grep -i '^content-length' | tail -1 | tr -dc '0-9')
log "remote size: $((TOTAL/1048576)) MB"
[ -f "$ZIP" ] && log "resuming from $(( $(stat -c%s "$ZIP") / 1048576 )) MB"

curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors \
     --speed-limit 10240 --speed-time 60 \
     -C - -o "$ZIP" "$URL"

GOT=$(stat -c%s "$ZIP")
if [ "$GOT" -ne "$TOTAL" ]; then
  log "ERROR: downloaded $GOT bytes, expected $TOTAL"; exit 1
fi
log "download complete and size-verified"

log "extracting into ${SNPEFF_DIR}"
unzip -q -o "$ZIP" -d "$SNPEFF_DIR"
[ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ] || { log "ERROR: data/${SNPEFF_DB} did not appear"; exit 1; }

log "OK - database installed: $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"
rm -f "$ZIP"
