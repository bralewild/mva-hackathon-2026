#!/usr/bin/env bash
# ==============================================================================
# 02a_download_snpeff_db.sh — Descarga robusta de la BD de snpEff
#
# SALIDA: ${SNPEFF_DIR}data/${SNPEFF_DB}/
#
# POR QUE EXISTE ESTE SCRIPT
# --------------------------
# El descargador interno de snpEff (`snpEff download`) no implementa reintento
# ni timeout: si se cae el socket, el proceso Java queda dormido para siempre
# sin error ni exit code. Nos paso: se congelo a los 285 MB de 770 y quedo
# 32 minutos sin avanzar al 0.5% de CPU.
#
# Este script usa curl con:
#   -C -                reanuda desde lo ya descargado
#   --retry 10          reintenta ante fallos transitorios
#   --speed-limit/-time aborta si baja de 10 KB/s por 60 s (y reintenta)
# ==============================================================================
source "$(dirname "$0")/00_config.sh"

SNPEFF_DIR=$(ls -d "$HOME"/micromamba/envs/bio/share/snpeff-*/ | head -1)
URL="https://snpeff-public.s3.amazonaws.com/databases/v5_4/snpEff_v5_4_${SNPEFF_DB}.zip"
ZIP="/tmp/snpEff_v5_4_${SNPEFF_DB}.zip"

if [ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ]; then
  log "BD ya instalada: $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"; exit 0
fi

# matar cualquier descarga colgada de snpEff
pkill -f 'snpEff.jar download' 2>/dev/null && log "descarga colgada de snpEff terminada"

TOTAL=$(curl -sIL "$URL" | tr -d '\r' | grep -i '^content-length' | tail -1 | tr -dc '0-9')
log "tamano remoto: $((TOTAL/1048576)) MB"
[ -f "$ZIP" ] && log "reanudando desde $(( $(stat -c%s "$ZIP") / 1048576 )) MB"

curl -L --fail --retry 10 --retry-delay 5 --retry-all-errors \
     --speed-limit 10240 --speed-time 60 \
     -C - -o "$ZIP" "$URL"

GOT=$(stat -c%s "$ZIP")
if [ "$GOT" -ne "$TOTAL" ]; then
  log "ERROR: bajados $GOT bytes de $TOTAL esperados"; exit 1
fi
log "descarga completa y verificada por tamano"

log "descomprimiendo en ${SNPEFF_DIR}"
unzip -q -o "$ZIP" -d "$SNPEFF_DIR"
[ -d "${SNPEFF_DIR}data/${SNPEFF_DB}" ] || { log "ERROR: no aparecio data/${SNPEFF_DB}"; exit 1; }

log "OK — BD instalada: $(du -sh "${SNPEFF_DIR}data/${SNPEFF_DB}" | cut -f1)"
rm -f "$ZIP"
