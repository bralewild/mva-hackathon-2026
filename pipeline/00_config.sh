#!/usr/bin/env bash
# ==============================================================================
# 00_config.sh - Shared pipeline configuration
#
# Usage:  source "$(dirname "$0")/00_config.sh"
#
# PORTABILITY
# -----------
# Nothing here is tied to a particular machine. The project root is derived from
# this script's own location, so the pipeline works from any checkout, on any
# Linux host (native, WSL, container or cloud VM).
#
# The data root defaults to ~/mva and can be pointed anywhere:
#
#     MVA_BASE=/scratch/mva bash pipeline/run_all.sh
#
# Code and data are kept apart deliberately. On Windows/WSL the code usually
# lives on the NTFS side and the data on ext4, because I/O across the 9p bridge
# to /mnt/c is 5-10x slower. On a native Linux host the split is optional but
# still useful: the data tree is what has to be deleted under the Data Use
# Agreement, and keeping it separate makes that a single command.
# ==============================================================================
set -euo pipefail

# Project root = the directory containing this script's parent.
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Data root. Override with MVA_BASE.
BASE="${MVA_BASE:-$HOME/mva}"
RAW="$BASE/data/raw"
REF="$BASE/data/ref"
ANNOT="$BASE/data/annot"
WORK="$BASE/work"
RESULTS="$BASE/results"
LOGS="$BASE/logs"

VCF_RAW="$RAW/WGS_EX2312012_HGWCNDSX7.vcf.gz"
SAMPLE="WGS_EX2312012"
PROBAND_ID="PROBAND01"          # required by the evaluator (submit_track1.py)
ASSEMBLY="GRCh38"
CONTIG_STYLE="ensembl"          # the VCF uses 1..22,X,Y WITHOUT a chr prefix
SUBMIT_CONTIG_PREFIX="chr"      # the evaluator expects chr15 (see docs/)

SNPEFF_DB="GRCh38.115"
THREADS="${MVA_THREADS:-$(nproc 2>/dev/null || echo 4)}"
JAVA_MEM="${MVA_JAVA_MEM:-12g}"

# Quality thresholds (rationale in docs/01_pipeline_flow.md)
MIN_GQ=20
MIN_DP=10
MAX_AF_RECESSIVE=0.01           # gnomAD, recessive model
MAX_AF_DOMINANT=0.0001          # gnomAD, dominant model

mkdir -p "$WORK" "$RESULTS" "$LOGS" "$REF" "$ANNOT"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# ------------------------------------------------------------------------------
# snpeff_dir - locate the directory holding snpEff.jar and its data/ folder.
#
# Do not assume a particular environment manager. Resolution order:
#   1. $SNPEFF_HOME, if the caller set it
#   2. derived from the `snpEff` wrapper on PATH (conda, micromamba, brew, apt)
#   3. a few conventional locations
#
# Prints the directory (with a trailing slash) and returns 0, or returns 1.
# ------------------------------------------------------------------------------
snpeff_dir() {
  local cand w d
  if [ -n "${SNPEFF_HOME:-}" ] && [ -f "${SNPEFF_HOME%/}/snpEff.jar" ]; then
    printf '%s/\n' "${SNPEFF_HOME%/}"; return 0
  fi
  if w=$(command -v snpEff 2>/dev/null); then
    d=$(dirname "$(readlink -f "$w")")
    for cand in "$d"/../share/snpeff-*/ "$d"/../share/snpEff/ "$d"/; do
      [ -f "${cand}snpEff.jar" ] && { printf '%s\n' "$(cd "$cand" && pwd)/"; return 0; }
    done
  fi
  for cand in "$HOME"/micromamba/envs/*/share/snpeff-*/ \
              "$HOME"/miniconda3/envs/*/share/snpeff-*/ \
              "$HOME"/mambaforge/envs/*/share/snpeff-*/ \
              /usr/share/snpeff/ /opt/snpEff/; do
    [ -f "${cand}snpEff.jar" ] && { printf '%s\n' "$(cd "$cand" && pwd)/"; return 0; }
  done
  return 1
}
