#!/usr/bin/env python3
"""
==============================================================================
06_validate_convergence.py - Compuerta de validacion de la busqueda ciega

ENTRADA : $RESULTS/05_ranked_genes.tsv
          $WORK/04_rare_candidates.tsv
SALIDA  : $RESULTS/06_convergence_report.txt

PROPOSITO
---------
Las etapas 01-05 corren SIN saber cual es la enfermedad. Esta etapa es la unica
que aplica conocimiento externo, y lo hace DESPUES de que el ranking esta
cerrado: verifica si el gen top-1 sostiene un par en heterocigosis compuesta
biologicamente plausible.

Es una COMPUERTA DE VALIDACION, no un filtro de busqueda. La distincion importa:
si el conocimiento de la enfermedad entrara antes del ranking, el resultado
seria circular y no probaria nada sobre el metodo.

CRITERIOS DE CONVERGENCIA
-------------------------
1. Separacion: el score del top-1 supera al del top-2 por un margen claro
2. Modelo:    el top-1 es AR_COMPOUND_HET o AR_HOMOZYGOUS
3. Par:       conserva >= 2 variantes raras tras el filtro de frecuencia
4. Severidad: al menos una variante de impacto HIGH
5. Evidencia: respaldo clinico independiente (ClinVar) o computacional (CADD)

Un fallo en cualquiera de los cinco no invalida el hallazgo, pero debe quedar
declarado en el reporte de metodos.
==============================================================================
"""
import csv
import os
import sys

BASE = os.path.expanduser("~/mva")
RANK = BASE + "/results/05_ranked_genes.tsv"
RARE = BASE + "/work/04_rare_candidates.tsv"
OUT = BASE + "/results/06_convergence_report.txt"

MARGEN_MINIMO = 0.15   # 15% de ventaja del top-1 sobre el top-2
CADD_ALTO = 20.0


def fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main():
    for p in (RANK, RARE):
        if not os.path.exists(p):
            sys.exit("falta " + p)

    ranked = list(csv.DictReader(open(RANK), delimiter="\t"))
    if not ranked:
        sys.exit("el ranking esta vacio")
    rare = list(csv.DictReader(open(RARE), delimiter="\t"))

    top = ranked[0]
    segundo = ranked[1] if len(ranked) > 1 else None
    gen = top["gene"]
    vs = [r for r in rare if r["gene"] == gen]

    s1 = fnum(top["pheno_score"])
    s2 = fnum(segundo["pheno_score"]) if segundo else 0.0
    margen = (s1 - s2) / s1 if s1 else 0.0

    highs = [v for v in vs if v.get("vep_impact") == "HIGH" or v.get("impact") == "HIGH"]
    clinvar = [v for v in vs if v.get("clinvar", ".") not in (".", "", None)
               and "pathogenic" in v.get("clinvar", "").lower()]
    cadds = [v for v in vs if fnum(v.get("cadd")) >= CADD_ALTO]

    criterios = [
        ("Separacion sobre el 2do", margen >= MARGEN_MINIMO,
         "{:.1f}% (minimo {:.0f}%)".format(margen * 100, MARGEN_MINIMO * 100)),
        ("Modelo recesivo", top["model"] in ("AR_COMPOUND_HET", "AR_HOMOZYGOUS"),
         top["model"]),
        ("Par de variantes raras", len(vs) >= 2, "{} variantes raras".format(len(vs))),
        ("Al menos una HIGH", len(highs) >= 1, "{} de impacto HIGH".format(len(highs))),
        ("Evidencia independiente", len(clinvar) >= 1 or len(cadds) >= 1,
         "{} en ClinVar patogenica, {} con CADD>={}".format(len(clinvar), len(cadds), CADD_ALTO)),
    ]
    paso = sum(1 for _, ok, _ in criterios if ok)

    L = []
    L.append("=" * 78)
    L.append(" 06 - COMPUERTA DE VALIDACION DE LA BUSQUEDA CIEGA")
    L.append("=" * 78)
    L.append("")
    L.append("Las etapas 01-05 corrieron sin conocer la enfermedad.")
    L.append("Esta etapa evalua el resultado DESPUES de cerrado el ranking.")
    L.append("")
    L.append("## Gen convergente")
    L.append("")
    L.append("  TOP-1 : {}   score {}".format(gen, top["pheno_score"]))
    if segundo:
        L.append("  TOP-2 : {}   score {}".format(segundo["gene"], segundo["pheno_score"]))
        L.append("  margen: {:.1f}%".format(margen * 100))
    L.append("")
    L.append("## Variantes raras en {}".format(gen))
    L.append("")
    for v in sorted(vs, key=lambda x: int(x["pos"])):
        af = v.get("gnomad_af") or "AUSENTE"
        L.append("  {}:{} {}>{}".format(v["chrom"], v["pos"], v["ref"], v["alt"]))
        L.append("     {}  |  impacto {}  |  CADD {}".format(
            v.get("vep_consequence", "."), v.get("vep_impact", "."), v.get("cadd") or "."))
        L.append("     {}".format(v.get("vep_hgvsp", ".")))
        L.append("     gnomAD {}  |  {}  |  ClinVar: {}".format(
            af, v.get("rsid", "."), v.get("clinvar", ".")))
        L.append("     genotipo {}  DP {}  GQ {}  AD {}".format(
            v.get("gt", "."), v.get("dp", "."), v.get("gq", "."), v.get("ad", ".")))
        L.append("")
    L.append("## Criterios de convergencia")
    L.append("")
    for nombre, ok, detalle in criterios:
        L.append("  [{}] {:<26} {}".format("OK" if ok else "--", nombre, detalle))
    L.append("")
    L.append("  {} de {} criterios cumplidos".format(paso, len(criterios)))
    L.append("")
    if paso == len(criterios):
        L.append("  >>> CONVERGENCIA CONFIRMADA <<<")
        L.append("")
        L.append("  El pipeline, sin conocer la enfermedad, ordeno {} en primer".format(gen))
        L.append("  lugar entre {} genes candidatos, partiendo de 5.012.204 variantes.".format(len(ranked)))
    else:
        L.append("  >>> CONVERGENCIA PARCIAL - declarar en el reporte de metodos <<<")
    L.append("")
    L.append("## Limitaciones que siguen vigentes")
    L.append("")
    L.append("  - Singleton: sin padres no se puede probar la fase por pedigri.")
    L.append("    El par se reporta como PRESUNTO en trans.")
    L.append("  - Solo SNV e indels: el VCF no contiene CNV, SV ni expansiones.")
    L.append("  - El filtro de impacto descarta intrones profundos y reguladoras.")
    L.append("  - Un gen sin anotaciones HPO obtiene score 0 aunque fuera el causal.")

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
