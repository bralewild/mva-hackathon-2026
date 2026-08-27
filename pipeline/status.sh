#!/usr/bin/env bash
# ==============================================================================
# status.sh — Estado actual del pipeline: que etapa esta hecha y que falta.
#
# Uso: wsl -d Ubuntu-24.04 -- bash -lc "bash .../pipeline/status.sh"
# ==============================================================================
source "$(dirname "$0")/00_config.sh"
set +o pipefail
SNPEFF_DIR=$(ls -d "$HOME"/micromamba/envs/bio/share/snpeff-*/ 2>/dev/null | head -1)

# $3 = tamano minimo en bytes para dar la etapa por completa.
# Necesario porque un bgzip vacio pesa 28 bytes y pasaria un simple `test -s`.
chk() {
  local min=${3:-1}
  if [ -f "$2" ] && [ "$(stat -c%s "$2")" -ge "$min" ]; then
    printf "  [OK]        %-32s %8s\n" "$1" "$(du -h "$2" | cut -f1)"
  else
    printf "  [pendiente] %-32s\n" "$1"
  fi
}

echo "=================================================================="
echo " ESTADO DEL PIPELINE — $(date +%H:%M:%S)"
echo "=================================================================="
echo
echo "## Etapas"
chk "01 QC baseline"           "$RESULTS/01_qc_baseline.txt"        500
chk "02 anotacion genome-wide" "$WORK/02_annotated.vcf.gz"     10000000
chk "03 modelos de herencia"   "$WORK/03_candidates.tsv"           1000
chk "04 frecuencia + clinica"  "$WORK/04_rare_candidates.tsv"      1000
chk "05 ranking fenotipico"    "$RESULTS/05_ranked_genes.tsv"       500
echo
echo "## Base de datos snpEff ($SNPEFF_DB)"
if [ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ]; then
  echo "  INSTALADA   $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"
else
  Z="/tmp/snpEff_v5_4_${SNPEFF_DB}.zip"
  [ -f "$Z" ] && echo "  descargando $(du -h "$Z" | cut -f1)" || echo "  no instalada"
fi
echo
echo "## Procesos activos"
PID=$(pgrep -f snpEff.jar | head -1)
if [ -n "$PID" ]; then
  ps -o etime=,%cpu= -p "$PID" | sed 's/^/  snpEff corriendo — tiempo:/'
else
  echo "  ninguno"
fi
echo
df -h /home | awk 'NR==2 {print "## Disco ext4: "$3" usado, "$4" libre"}'

# --- progreso real de la anotacion (opcional: status.sh --progress) ---
# Descomprime ~300 MB, tarda ~30 s, por eso no corre por defecto.
# zcat tolera el bgzip truncado que snpEff esta escribiendo; bgzip -dc no.
if [ "${1:-}" = "--progress" ] && [ -f "$WORK/02_annotated.vcf.gz" ]; then
  echo
  echo "## Progreso real de la anotacion (contando variantes)"
  N=$(zcat "$WORK/02_annotated.vcf.gz" 2>/dev/null | grep -vc '^#')
  LAST=$(zcat "$WORK/02_annotated.vcf.gz" 2>/dev/null | tail -1 | cut -f1,2)
  TOTAL=5012204
  echo "  anotadas : $N de $TOTAL"
  echo "  posicion : cromosoma $(echo "$LAST" | cut -f1), pos $(echo "$LAST" | cut -f2)"
  awk -v n="$N" -v t="$TOTAL" 'BEGIN{
    p=n*100/t; b=int(p/4);
    s=""; for(i=0;i<b;i++)s=s"#"; for(i=b;i<25;i++)s=s".";
    printf "  [%s] %.1f %%\n", s, p }'
fi
