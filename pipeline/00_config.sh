#!/usr/bin/env bash
# ==============================================================================
# 00_config.sh - Shared pipeline configuration
#
# Usage:  source "$(dirname "$0")/00_config.sh"
#
# Code -> /mnt/c/Users/user/Documents/real-kid-mva-hackathon  (version-controlled)
# Data -> ~/mva                                                (ext4, fast I/O)
# ==============================================================================
set -euo pipefail

PROJECT="/mnt/c/Users/user/Documents/real-kid-mva-hackathon"
BASE="$HOME/mva"
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
THREADS=28
JAVA_MEM="12g"

# Quality thresholds (rationale in docs/01_pipeline_flow.md)
MIN_GQ=20
MIN_DP=10
MAX_AF_RECESSIVE=0.01           # gnomAD, recessive model
MAX_AF_DOMINANT=0.0001          # gnomAD, dominant model

mkdir -p "$WORK" "$RESULTS" "$LOGS" "$REF" "$ANNOT"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
