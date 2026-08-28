#!/usr/bin/env python3
"""
==============================================================================
04b_seed_cache.py - Reconstruye el cache de VEP desde una corrida previa

ENTRADA : $WORK/04_annotated_candidates.tsv
SALIDA  : $WORK/04_vep_cache.jsonl

PROPOSITO
---------
Utilidad de recuperacion. La primera version de 04_frequency_clinical.py no
persistia resultados por lote: escribia todo al final. En esa corrida, 3 de 47
lotes fallaron con HTTP 500 de Ensembl y sus 600 variantes quedaron sin anotar,
pero las otras 8.784 si se anotaron bien.

Reejecutar la etapa entera para recuperar 600 variantes seria tirar 40 minutos
de consultas validas. Este script siembra el cache incremental con lo que ya
se anoto correctamente, para que al volver a correr la etapa 04 solo se pidan
a Ensembl las variantes que realmente faltan.

Detecta "sin respuesta de VEP" cuando TODOS los campos anotados vienen vacios:
sin consecuencia, sin rsID, sin transcrito MANE y sin CADD.

Se puede correr las veces que haga falta: es idempotente.
==============================================================================
"""
import csv
import json
import os
import sys

BASE = os.path.expanduser("~/mva")
IN_TSV = BASE + "/work/04_annotated_candidates.tsv"
CACHE = BASE + "/work/04_vep_cache.jsonl"


def sin_respuesta(r):
    """True si VEP no devolvio nada para esta variante."""
    vacio = (".", "", None)
    return (r.get("vep_consequence") in vacio
            and r.get("rsid") in vacio
            and r.get("mane") in vacio
            and not r.get("cadd"))


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("falta " + IN_TSV + " - no hay corrida previa para sembrar")

    ya = set()
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as f:
            for line in f:
                try:
                    ya.add(json.loads(line)["key"])
                except Exception:
                    continue
        print("cache existente: {:,} entradas".format(len(ya)))

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    sembradas = faltantes = repetidas = 0

    with open(CACHE, "a", encoding="utf-8") as out:
        for r in rows:
            key = "{}_{}_{}_{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
            if sin_respuesta(r):
                faltantes += 1
                continue
            if key in ya:
                repetidas += 1
                continue
            af = r.get("gnomad_af")
            cadd = r.get("cadd")
            ann = {
                "gnomad_af": float(af) if af not in ("", None) else None,
                "rsid": r.get("rsid", "."),
                "clinvar": r.get("clinvar", "."),
                "mane": r.get("mane", "."),
                "hgvsc": r.get("vep_hgvsc", "."),
                "hgvsp": r.get("vep_hgvsp", "."),
                "cadd": float(cadd) if cadd not in ("", None) else None,
                "consequence": r.get("vep_consequence", "."),
                "impact": r.get("vep_impact", "."),
                "gene": r.get("gene", "."),
            }
            out.write(json.dumps({"key": key, "ann": ann}) + "\n")
            ya.add(key)
            sembradas += 1

    print("=" * 60)
    print(" 04b - SIEMBRA DEL CACHE DE VEP")
    print("=" * 60)
    print("  filas leidas          {:>8,}".format(len(rows)))
    print("  sembradas al cache    {:>8,}".format(sembradas))
    print("  ya estaban            {:>8,}".format(repetidas))
    print("  SIN anotar (a pedir)  {:>8,}".format(faltantes))
    print()
    print("  -> " + CACHE)
    print()
    print("  Ahora corre de nuevo 04_frequency_clinical.py:")
    print("  solo va a consultar a Ensembl las {:,} que faltan.".format(faltantes))


if __name__ == "__main__":
    main()
