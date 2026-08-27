#!/usr/bin/env bash
# ==============================================================================
# 02_annotate_genomewide.sh — Anotacion funcional CIEGA de todo el genoma
#
# ENTRADA : $VCF_RAW  (5.012.204 variantes)
# SALIDA  : $WORK/02_annotated.vcf.gz + .tbi
#           $RESULTS/02_snpeff_summary.html
#
# PROPOSITO: anotar TODAS las variantes sin ningun sesgo hacia genes candidatos.
#            Este es el paso que hace CIEGA la busqueda: el pipeline no sabe que
#            la enfermedad es MVA ni que existen BUB1B / CEP57 / TRIP13.
#
# BD ELEGIDA: GRCh38.115 (Ensembl completo) en vez de GRCh38.mane.*
#            MANE cubre ~19.3k genes con un transcrito cada uno — ideal para
#            REPORTAR, pero en una busqueda ciega prima la cobertura completa.
#            El reporte final usa transcritos MANE via VEP (etapa 04).
#
# LECCION APRENDIDA: snpEff sale con codigo 0 aunque falle (opcion invalida,
#            descarga muerta). NO alcanza con 'set -e': validamos la salida.
# ==============================================================================
source "$(dirname "$0")/00_config.sh"

# --- dependencia: base de datos ---
bash "$(dirname "$0")/02a_download_snpeff_db.sh"

OUT="$WORK/02_annotated.vcf.gz"
log "anotando 5.012.204 variantes con snpEff $SNPEFF_DB (~15-40 min)"

snpEff -Xmx${JAVA_MEM} ann \
  -noLog \
  -stats "$RESULTS/02_snpeff_summary.html" \
  "$SNPEFF_DB" "$VCF_RAW" \
  2> "$LOGS/02_snpeff.log" \
  | bgzip -@ 8 > "$OUT"

# --- VALIDACION: no confiar en el exit code ---
SIZE=$(stat -c%s "$OUT")
if [ "$SIZE" -lt 10000000 ]; then
  log "ERROR: salida de ${SIZE} bytes — snpEff fallo. Ultimas lineas del log:"
  tail -20 "$LOGS/02_snpeff.log"
  exit 1
fi
N=$(bcftools view -H "$OUT" | wc -l)
if [ "$N" -lt 5000000 ]; then
  log "ERROR: $N variantes anotadas, se esperaban ~5.012.204"; exit 1
fi

tabix -f -p vcf "$OUT"
log "OK — $N variantes anotadas"
log "salida  -> $OUT ($(du -h "$OUT" | cut -f1))"
log "resumen -> $RESULTS/02_snpeff_summary.html"
