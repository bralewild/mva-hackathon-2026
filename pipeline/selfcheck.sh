#!/usr/bin/env bash
# ==============================================================================
# selfcheck.sh - Preflight validation
#
# Run this before the pipeline (or before trusting a checkout). It verifies the
# environment, the inputs and the repository itself, and reports what is wrong
# rather than failing halfway through a one-hour annotation.
#
#   bash pipeline/selfcheck.sh
#
# Exit code 0 if everything required is present, 1 otherwise.
# ==============================================================================
set -uo pipefail
D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$D/00_config.sh"
set +e

PASS=0; WARN=0; FAIL=0
ok()   { printf "  \033[32m[ ok ]\033[0m   %-34s %s\n" "$1" "${2:-}"; PASS=$((PASS+1)); }
warn() { printf "  \033[33m[warn]\033[0m   %-34s %s\n" "$1" "${2:-}"; WARN=$((WARN+1)); }
fail() { printf "  \033[31m[FAIL]\033[0m   %-34s %s\n" "$1" "${2:-}"; FAIL=$((FAIL+1)); }

echo "=================================================================="
echo " PREFLIGHT SELF-CHECK"
echo "=================================================================="
echo
echo "## Resolved paths"
echo "  PROJECT  $PROJECT"
echo "  BASE     $BASE   (override with MVA_BASE)"
echo "  THREADS  $THREADS"
echo

echo "## Required tools"
for t in bcftools samtools tabix bgzip snpEff java python3 curl unzip; do
  p=$(command -v "$t" 2>/dev/null)
  if [ -n "$p" ]; then ok "$t" "$p"; else fail "$t" "not on PATH"; fi
done
echo

echo "## Versions"
v=$(bcftools --version 2>/dev/null | head -1 | awk '{print $2}')
[ -n "$v" ] && ok "bcftools" "$v" || fail "bcftools" "no version"
jv=$(java -version 2>&1 | head -1 | grep -oE '"[0-9]+' | tr -d '"')
if [ -n "$jv" ] && [ "$jv" -ge 21 ] 2>/dev/null; then
  ok "java" "$jv (snpEff 5.4 needs 21+)"
else
  fail "java" "${jv:-unknown} - snpEff 5.4 requires 21+"
fi
pv=$(python3 -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null)
[ -n "$pv" ] && ok "python3" "$pv" || fail "python3" "missing"
echo

echo "## snpEff installation"
if SD=$(snpeff_dir 2>/dev/null); then
  ok "snpEff.jar" "$SD"
  if [ -d "${SD}data/${SNPEFF_DB}" ]; then
    ok "database $SNPEFF_DB" "$(du -sh "${SD}data/${SNPEFF_DB}" 2>/dev/null | cut -f1)"
  else
    warn "database $SNPEFF_DB" "not installed - stage 02a will download ~770 MB"
  fi
else
  fail "snpEff directory" "not found - set SNPEFF_HOME"
fi
echo

echo "## Python stages import cleanly"
for f in "$D"/0*.py "$D"/_paths.py; do
  n=$(basename "$f")
  if python3 -c "import ast,sys; ast.parse(open('$f').read())" 2>/dev/null; then
    ok "$n" "syntax"
  else
    fail "$n" "syntax error"
  fi
done
if (cd "$D" && python3 -c "import _paths" 2>/dev/null); then
  ok "_paths import" "shared path module resolves"
else
  fail "_paths import" "shared path module fails"
fi
echo

echo "## Shell stages parse"
for f in "$D"/*.sh; do
  n=$(basename "$f")
  bash -n "$f" 2>/dev/null && ok "$n" "syntax" || fail "$n" "syntax error"
done
echo

echo "## Input data"
if [ -f "$VCF_RAW" ]; then
  ok "VCF" "$(du -h "$VCF_RAW" | cut -f1)"
  [ -f "${VCF_RAW}.tbi" ] && ok "VCF index" "present" || fail "VCF index" "missing .tbi"
  if bcftools view -h "$VCF_RAW" >/dev/null 2>&1; then
    ok "VCF readable" "$(bcftools query -l "$VCF_RAW" | head -1)"
  else
    fail "VCF readable" "bcftools cannot parse it"
  fi
else
  fail "VCF" "missing at $VCF_RAW - supply your own authorised copy"
fi
[ -f "$RAW/Challenge_Clinical_Phenotype_1.docx" ] \
  && ok "clinical document" "present" \
  || fail "clinical document" "missing - stage 00b needs it"
echo

echo "## Disk"
free_gb=$(df -BG --output=avail "$BASE" 2>/dev/null | tail -1 | tr -dc '0-9')
if [ -n "$free_gb" ] && [ "$free_gb" -ge 40 ]; then
  ok "free space" "${free_gb} GB (>= 40 GB)"
elif [ -n "$free_gb" ]; then
  warn "free space" "${free_gb} GB - 40 GB recommended"
fi
echo

echo "## Network"
for host in rest.ensembl.org snpeff-public.s3.amazonaws.com purl.obolibrary.org; do
  code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 15 "https://$host" 2>/dev/null)
  [ "$code" != "000" ] && ok "$host" "HTTP $code" || fail "$host" "unreachable"
done
echo

echo "## Repository hygiene"
if command -v git >/dev/null && git -C "$PROJECT" rev-parse --git-dir >/dev/null 2>&1; then
  leak=$(git -C "$PROJECT" ls-files | grep -cE '\.(vcf|bam|cram|fastq|docx)$|patient_hpo')
  [ "$leak" -eq 0 ] && ok "no patient data tracked" || fail "patient data tracked" "$leak files"
  # Exclude this file (it contains the pattern itself) and compiled bytecode,
  # which embeds the resolved path but is never committed.
  hard=$(grep -rlI --exclude='selfcheck.sh' --exclude-dir='__pycache__' \
           -e '/mnt/c/Users' -e '/home/[a-z]*/mva"' "$D" 2>/dev/null | wc -l)
  if [ "$hard" -eq 0 ]; then
    ok "no machine-specific paths"
  else
    warn "hardcoded paths" "$hard files"
    grep -rlI --exclude='selfcheck.sh' --exclude-dir='__pycache__' \
      -e '/mnt/c/Users' -e '/home/[a-z]*/mva"' "$D" 2>/dev/null | sed 's/^/           /'
  fi
else
  warn "git" "not a repository - skipping hygiene checks"
fi

echo
echo "=================================================================="
printf " %d passed, %d warnings, %d failures\n" "$PASS" "$WARN" "$FAIL"
echo "=================================================================="
[ "$FAIL" -eq 0 ] && echo " READY" || echo " NOT READY - resolve the failures above"
exit $([ "$FAIL" -eq 0 ] && echo 0 || echo 1)
