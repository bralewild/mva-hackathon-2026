#!/usr/bin/env bash
# ==============================================================================
# 99_data_inventory.sh — Inventario de TODAS las fuentes de datos
#
# SALIDA: $RESULTS/99_data_inventory.txt
#
# PROPOSITO
# ---------
# 1. Saber exactamente que hay, donde y cuanto pesa.
# 2. Cumplir el Data Use Agreement: distingue datos DEL PACIENTE (borrables al
#    terminar el hackathon) de recursos PUBLICOS (referencias, BDs) que no
#    contienen informacion del paciente y no requieren borrado.
#
# Ejecutar antes de notificar a MVAHackathon2026@synapse.org.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail
OUT="$RESULTS/99_data_inventory.txt"
SNPEFF_DIR=$(ls -d "$HOME"/micromamba/envs/bio/share/snpeff-*/ 2>/dev/null | head -1)

sz() { [ -e "$1" ] && du -sh "$1" 2>/dev/null | cut -f1 || echo "-"; }

{
echo "=============================================================================="
echo " INVENTARIO DE FUENTES DE DATOS — $(date -Iseconds)"
echo "=============================================================================="
echo
echo "###############################################################"
echo "# A) DATOS DEL PACIENTE  — SUJETOS AL DUA, BORRAR AL TERMINAR  #"
echo "###############################################################"
echo
echo "  [$( sz "$RAW" )]  $RAW    (descargado de HuggingFace, gated)"
[ -d "$RAW" ] && ls -lh "$RAW" | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  [$( sz "$WORK" )]  $WORK    (intermedios derivados del VCF)"
[ -d "$WORK" ] && ls -lh "$WORK" 2>/dev/null | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  [$( sz "$RESULTS" )]  $RESULTS    (resultados; contienen coordenadas del paciente)"
[ -d "$RESULTS" ] && ls -lh "$RESULTS" 2>/dev/null | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  ORIGEN REMOTO: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data"
echo "                 (gated; 11 archivos, ~85 GB. Descargamos solo 4, ~318 MB)"
echo "  NO DESCARGADO: los 8 FASTQ (~84 GB) siguen solo en HuggingFace."
echo
echo "###############################################################"
echo "# B) RECURSOS PUBLICOS  — sin datos del paciente, NO borrar    #"
echo "###############################################################"
echo
echo "  LOCALES (en disco):"
echo "    [$( sz "${SNPEFF_DIR}data" )]  ${SNPEFF_DIR}data"
[ -d "${SNPEFF_DIR}data" ] && ls "${SNPEFF_DIR}data" 2>/dev/null | sed 's/^/        BD: /'
echo "    [$( sz "$REF" )]  $REF    (referencia genomica, si se descarga)"
echo "    [$( sz "$ANNOT" )]  $ANNOT    (anotaciones adicionales)"
echo "    [$( sz "$HOME/micromamba" )]  $HOME/micromamba    (entorno de software)"
echo
echo "  REMOTOS (APIs, no se almacena nada):"
echo "    Ensembl REST   https://rest.ensembl.org           coordenadas de genes"
echo "    Ensembl VEP    https://rest.ensembl.org/vep/...   consecuencia, CADD, gnomAD"
echo "    NCBI ClinVar   https://eutils.ncbi.nlm.nih.gov    significancia clinica"
echo "    snpEff DBs     https://snpeff-public.s3.amazonaws.com"
echo
echo "###############################################################"
echo "# C) CODIGO  — sin datos del paciente, versionable             #"
echo "###############################################################"
echo
echo "    $PROJECT"
echo "      .gitignore bloquea: *.vcf* *.bam *.cram *.fastq* *.docx data/ work/"
echo
echo "###############################################################"
echo "# D) COMANDO DE BORRADO (DUA)                                  #"
echo "###############################################################"
echo
echo "    rm -rf $RAW $WORK $RESULTS"
echo "    # luego notificar a MVAHackathon2026@synapse.org"
echo
echo "## Espacio total en ext4"
df -h /home | awk 'NR==2 {print "    usado: "$3"   libre: "$4}'
} > "$OUT"

cat "$OUT"
log "escrito -> $OUT"
