#!/usr/bin/env bash
# ==============================================================================
# sync_results.sh - Mirror lightweight results from WSL to the Windows project
#
# INPUT  : $RESULTS/*        (on ext4)
# OUTPUT : $PROJECT/results/*  (visible in the project, version-controllable)
#
# PURPOSE: heavy data stays on ext4 for performance, but the reports are a few
#          kilobytes and should be visible from the project folder. Only files
#          under 5 MB are copied, and never VCF/BAM/FASTQ.
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
  [ "$sz" -gt 5242880 ] && { log "skipped (>5MB): $(basename "$f")"; continue; }
  cp -f "$f" "$DST/" && n=$((n+1))
done
log "mirrored $n files -> $DST"
ls -lh "$DST" | awk 'NR>1 {printf "  %8s  %s\n", $5, $9}'
