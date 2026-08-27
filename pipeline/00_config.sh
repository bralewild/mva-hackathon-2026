#!/usr/bin/env bash
# ==============================================================================
# 00_config.sh — Configuracion compartida del pipeline
#
# Uso:  source "$(dirname "$0")/00_config.sh"
#
# Codigo  -> /mnt/c/Users/user/Documents/real-kid-mva-hackathon  (versionable)
# Datos   -> ~/mva                                                (ext4, rapido)
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
PROBAND_ID="PROBAND01"          # exigido por el evaluador (submit_track1.py)
ASSEMBLY="GRCh38"
CONTIG_STYLE="ensembl"          # el VCF usa 1..22,X,Y SIN prefijo chr
SUBMIT_CONTIG_PREFIX="chr"      # el evaluador espera chr15 (ver docs/)

SNPEFF_DB="GRCh38.115"
THREADS=28
JAVA_MEM="12g"

# Umbrales de calidad (justificados en docs/02_metodologia.md)
MIN_GQ=20
MIN_DP=10
MAX_AF_RECESSIVE=0.01           # gnomAD, modelo recesivo
MAX_AF_DOMINANT=0.0001          # gnomAD, modelo dominante

mkdir -p "$WORK" "$RESULTS" "$LOGS" "$REF" "$ANNOT"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
