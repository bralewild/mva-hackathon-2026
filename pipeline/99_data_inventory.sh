#!/usr/bin/env bash
# ==============================================================================
# 99_data_inventory.sh - Inventory of EVERY data source
#
# OUTPUT: $RESULTS/99_data_inventory.txt
#
# PURPOSE
# -------
# 1. Know exactly what exists, where, and how large it is.
# 2. Satisfy the Data Use Agreement: distinguishes PATIENT data (to be deleted
#    at hackathon close, including derived datasets) from PUBLIC resources
#    (references, annotation databases) that contain no patient information.
#
# Run this before notifying the organisers that data has been deleted.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail
OUT="$RESULTS/99_data_inventory.txt"
SNPEFF_DIR=$(snpeff_dir 2>/dev/null || echo "")

sz() { [ -e "$1" ] && du -sh "$1" 2>/dev/null | cut -f1 || echo "-"; }

{
echo "=============================================================================="
echo " DATA SOURCE INVENTORY - $(date -Iseconds)"
echo "=============================================================================="
echo
echo "###############################################################"
echo "# A) PATIENT DATA - SUBJECT TO THE DUA, DELETE AT CLOSE        #"
echo "###############################################################"
echo
echo "  [$( sz "$RAW" )]  $RAW    (downloaded from the gated HuggingFace dataset)"
[ -d "$RAW" ] && ls -lh "$RAW" | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  [$( sz "$WORK" )]  $WORK    (intermediates derived from the VCF)"
[ -d "$WORK" ] && ls -lh "$WORK" 2>/dev/null | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  [$( sz "$RESULTS" )]  $RESULTS    (results; contain patient coordinates)"
[ -d "$RESULTS" ] && ls -lh "$RESULTS" 2>/dev/null | awk 'NR>1 {printf "        %8s  %s\n", $5, $9}'
echo
echo "  [$( sz "$PROJECT/results" )]  $PROJECT/results    (mirrored reports - also derived)"
echo
echo "  REMOTE ORIGIN: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data"
echo "                 (gated; 11 files, ~85 GB. Only 4 downloaded, ~318 MB)"
echo "  NOT DOWNLOADED: the 8 FASTQ files (~84 GB) remain only on HuggingFace."
echo
echo "###############################################################"
echo "# B) PUBLIC RESOURCES - no patient data, DO NOT delete         #"
echo "###############################################################"
echo
echo "  LOCAL (on disk):"
echo "    [$( sz "${SNPEFF_DIR}data" )]  ${SNPEFF_DIR}data"
[ -d "${SNPEFF_DIR}data" ] && ls "${SNPEFF_DIR}data" 2>/dev/null | sed 's/^/        DB: /'
echo "    [$( sz "$REF" )]  $REF    (genome reference, if downloaded)"
echo "    [$( sz "$ANNOT" )]  $ANNOT    (HPO ontology and annotations)"
echo "    [$( sz "$(dirname "$(dirname "${SNPEFF_DIR:-/nonexistent}")")" )]  software environment (conda/micromamba prefix)"
echo
echo "  REMOTE (APIs - nothing stored):"
echo "    Ensembl REST   https://rest.ensembl.org           gene coordinates"
echo "    Ensembl VEP    https://rest.ensembl.org/vep/...   consequence, CADD, gnomAD"
echo "    NCBI ClinVar   https://eutils.ncbi.nlm.nih.gov    clinical significance"
echo "    HPO            https://purl.obolibrary.org/obo/   ontology and annotations"
echo "    snpEff DBs     https://snpeff-public.s3.amazonaws.com"
echo
echo "###############################################################"
echo "# C) CODE - no patient data, version-controlled                #"
echo "###############################################################"
echo
echo "    $PROJECT"
echo "      .gitignore blocks: *.vcf* *.bam *.cram *.fastq* *.docx"
echo "                         patient_hpo.tsv data/ work/ submissions/bralewild_*"
echo
echo "###############################################################"
echo "# D) DELETION COMMAND (DUA)                                    #"
echo "###############################################################"
echo
echo "    rm -rf $RAW $WORK $RESULTS"
echo "    rm -rf $PROJECT/results"
echo
echo "    Then notify BOTH official addresses:"
echo "      RarediseaserealkidMVAhackathon2026@synapse.org   (Official Rules)"
echo "      MVAHackathon2026@synapse.org                     (dataset gating form)"
echo
echo "## Total space on ext4"
df -h /home | awk 'NR==2 {print "    used: "$3"   free: "$4}'
} > "$OUT"

cat "$OUT"
log "written -> $OUT"
