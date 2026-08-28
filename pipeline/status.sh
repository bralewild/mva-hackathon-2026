#!/usr/bin/env bash
# ==============================================================================
# status.sh - Pipeline status: which stages are done, which are pending.
#
# Usage: wsl -d Ubuntu-24.04 -- bash -lc "bash .../pipeline/status.sh [--progress]"
#
#   --progress   also counts annotated variants in the in-flight output.
#                Decompresses ~500 MB, takes ~30 s, so it is opt-in.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail
SNPEFF_DIR=$(ls -d "$HOME"/micromamba/envs/bio/share/snpeff-*/ 2>/dev/null | head -1)

# $3 = minimum size in bytes for a stage to count as complete.
# Needed because an empty bgzip file is 28 bytes and would pass a plain `test -s`.
chk() {
  local min=${3:-1}
  if [ -f "$2" ] && [ "$(stat -c%s "$2")" -ge "$min" ]; then
    printf "  [OK]        %-32s %8s\n" "$1" "$(du -h "$2" | cut -f1)"
  else
    printf "  [pending]   %-32s\n" "$1"
  fi
}

echo "=================================================================="
echo " PIPELINE STATUS - $(date +%H:%M:%S)"
echo "=================================================================="
echo
echo "## Stages"
chk "01 QC baseline"            "$RESULTS/01_qc_baseline.txt"        500
chk "02 genome-wide annotation" "$WORK/02_annotated.vcf.gz"     10000000
chk "03 inheritance models"     "$WORK/03_candidates.tsv"           1000
chk "04 frequency + clinical"   "$WORK/04_rare_candidates.tsv"      1000
chk "05 phenotype ranking"      "$RESULTS/05_ranked_genes.tsv"       500
chk "06 convergence gate"       "$RESULTS/06_convergence_report.txt" 500
chk "07 secondary findings"     "$RESULTS/07_secondary_findings.txt" 300
echo
echo "## snpEff database ($SNPEFF_DB)"
if [ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ]; then
  echo "  INSTALLED   $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"
else
  Z="/tmp/snpEff_v5_4_${SNPEFF_DB}.zip"
  [ -f "$Z" ] && echo "  downloading $(du -h "$Z" | cut -f1)" || echo "  not installed"
fi
echo
echo "## Running processes"
PID=$(pgrep -f snpEff.jar | head -1)
if [ -n "$PID" ]; then
  ps -o etime=,%cpu= -p "$PID" | sed 's/^/  snpEff running - elapsed:/'
else
  echo "  none"
fi
echo
df -h /home | awk 'NR==2 {print "## Disk (ext4): "$3" used, "$4" free"}'

# --- real annotation progress (opt-in: status.sh --progress) ---
# Decompresses ~500 MB, ~30 s, so it does not run by default.
# zcat tolerates the truncated bgzip snpEff is still writing; bgzip -dc does not.
if [ "${1:-}" = "--progress" ] && [ -f "$WORK/02_annotated.vcf.gz" ]; then
  echo
  echo "## Real annotation progress (counting variants)"
  N=$(zcat "$WORK/02_annotated.vcf.gz" 2>/dev/null | grep -vc '^#')
  LAST=$(zcat "$WORK/02_annotated.vcf.gz" 2>/dev/null | tail -1 | cut -f1,2)
  TOTAL=5012204
  echo "  annotated : $N of $TOTAL"
  echo "  position  : chromosome $(echo "$LAST" | cut -f1), pos $(echo "$LAST" | cut -f2)"
  awk -v n="$N" -v t="$TOTAL" 'BEGIN{
    p=n*100/t; b=int(p/4);
    s=""; for(i=0;i<b;i++)s=s"#"; for(i=b;i<25;i++)s=s".";
    printf "  [%s] %.1f %%\n", s, p }'
fi
