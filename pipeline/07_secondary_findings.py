#!/usr/bin/env python3
"""
==============================================================================
07_secondary_findings.py - Hallazgos secundarios clinicamente accionables

ENTRADA : $WORK/04_rare_candidates.tsv
          $RESULTS/05_ranked_genes.tsv
SALIDA  : $RESULTS/07_secondary_findings.txt
          $WORK/07_secondary_candidates.tsv

PROPOSITO
---------
El FAQ del challenge aclara que los hallazgos secundarios NO afectan el score
automatico y se apartan para revision cualitativa del panel de jueces. Es una
via de evaluacion adicional sin riesgo.

Independientemente del concurso, reportar hallazgos secundarios accionables es
practica clinica estandar: si en el genoma de un paciente aparece una variante
patogenica en un gen donde existe intervencion medica disponible, corresponde
señalarlo aunque no explique el cuadro que motivo el estudio.

CRITERIOS
---------
Se marca como hallazgo secundario candidato una variante que:
  1. Este en un gen de la lista ACMG SF v3.2 (genes accionables consensuados)
     o en un gen con enfermedad tratable de la lista curada de abajo
  2. NO sea el gen causal primario (el top-1 del ranking fenotipico)
  3. Tenga impacto HIGH, o MODERATE con CADD >= 20
  4. Sea rara (ya garantizado por la etapa 04)

IMPORTANTE - LIMITE DECLARADO
-----------------------------
Esto NO es un reporte clinico. Es una lista de candidatos para revision humana.
Un hallazgo secundario real exige confirmacion ortogonal, evaluacion ACMG
completa y consejo genetico. Se reporta como hipotesis, nunca como diagnostico.

FUENTE de la lista: ACMG SF v3.2 (Miller et al., Genet Med 2023).
Se incluye embebida y no por descarga para que el pipeline sea reproducible
sin depender de que una URL siga viva.
==============================================================================
"""
import csv
import os
import sys

BASE = os.path.expanduser("~/mva")
RARE = BASE + "/work/04_rare_candidates.tsv"
RANK = BASE + "/results/05_ranked_genes.tsv"
OUT_TSV = BASE + "/work/07_secondary_candidates.tsv"
OUT_TXT = BASE + "/results/07_secondary_findings.txt"

CADD_MIN = 20.0

# ACMG SF v3.2 - genes con hallazgos secundarios reportables (Miller et al. 2023)
ACMG_SF = {
    # Cancer
    "BRCA1", "BRCA2", "PALB2", "TP53", "STK11", "MLH1", "MSH2", "MSH6", "PMS2",
    "APC", "MUTYH", "BMPR1A", "SMAD4", "GREM1", "VHL", "MEN1", "RET", "PTEN",
    "RB1", "SDHD", "SDHAF2", "SDHC", "SDHB", "MAX", "TMEM127", "NF2", "TSC1",
    "TSC2", "WT1", "MET", "BAP1", "CDH1", "CDKN2A", "DICER1", "FH", "MITF",
    "PRKAR1A", "RUNX1", "SDHA", "TRIM37",
    # Cardiovascular
    "MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL3", "ACTC1", "PRKAG2",
    "GLA", "MYL2", "LMNA", "PKP2", "DSP", "DSC2", "TMEM43", "DSG2", "RYR2",
    "KCNQ1", "KCNH2", "SCN5A", "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2",
    "MYH11", "COL3A1", "BAG3", "DES", "FLNC", "TTN", "TNNC1", "CASQ2", "TRDN",
    "CALM1", "CALM2", "CALM3", "HFE",
    # Metabolico / otros
    "OTC", "ATP7B", "BTD", "GAA", "APOB", "LDLR", "PCSK9", "RYR1", "CACNA1S",
    "ACVRL1", "ENG", "RPE65", "TTR", "GCK", "HNF1A", "HNF1B", "HNF4A",
}

# Genes fuera de ACMG SF pero con enfermedad TRATABLE o manejable.
# Se reportan por separado y con menor peso: no son consenso, son criterio propio.
TRATABLES = {
    "SERPINA1": "deficit de alfa-1 antitripsina - manejo pulmonar/hepatico, terapia de aumento",
    "CBS": "homocistinuria - responde a piridoxina/betaina y dieta restringida en metionina",
    "PAH": "fenilcetonuria - manejo dietario",
    "GALT": "galactosemia - manejo dietario",
    "SLC22A5": "deficit primario de carnitina - suplementacion con L-carnitina",
    "ATM": "ataxia-telangiectasia / riesgo oncologico en portadores - vigilancia",
    "MEFV": "fiebre mediterranea familiar - colchicina",
    "G6PD": "deficit de G6PD - evitar farmacos oxidantes",
    "ACADM": "deficit de MCAD - evitar ayuno prolongado",
    "F5": "factor V Leiden - manejo del riesgo trombotico",
    "F2": "protrombina G20210A - manejo del riesgo trombotico",
}


def fnum(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def main():
    for p in (RARE, RANK):
        if not os.path.exists(p):
            sys.exit("falta " + p)

    ranked = list(csv.DictReader(open(RANK), delimiter="\t"))
    primario = ranked[0]["gene"] if ranked else None

    rows = list(csv.DictReader(open(RARE), delimiter="\t"))

    acmg, trat = [], []
    for r in rows:
        gen = r.get("gene", "")
        if not gen or gen == primario:
            continue
        impacto = r.get("vep_impact") or r.get("impact") or "."
        cadd = fnum(r.get("cadd"))
        if not (impacto == "HIGH" or (impacto == "MODERATE" and cadd >= CADD_MIN)):
            continue
        if gen in ACMG_SF:
            acmg.append(r)
        elif gen in TRATABLES:
            trat.append(r)

    def orden(r):
        return (-(1 if (r.get("vep_impact") or r.get("impact")) == "HIGH" else 0),
                -fnum(r.get("cadd")))

    acmg.sort(key=orden)
    trat.sort(key=orden)

    if acmg or trat:
        cols = list((acmg or trat)[0].keys())
        with open(OUT_TSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in acmg + trat:
                w.writerow(r)

    def bloque(titulo, lista, notas=None):
        L = ["", titulo, "-" * len(titulo), ""]
        if not lista:
            L.append("  (ninguno)")
            return L
        for r in lista:
            gen = r["gene"]
            af = r.get("gnomad_af") or "AUSENTE"
            L.append("  {}  {}:{} {}>{}".format(gen, r["chrom"], r["pos"], r["ref"], r["alt"]))
            L.append("     {}  |  impacto {}  |  CADD {}".format(
                r.get("vep_consequence", "."),
                r.get("vep_impact") or r.get("impact") or ".",
                r.get("cadd") or "."))
            L.append("     {}".format(r.get("vep_hgvsp", ".")))
            L.append("     gnomAD {}  |  {}  |  ClinVar: {}".format(
                af, r.get("rsid", "."), r.get("clinvar", ".")))
            L.append("     genotipo {}  DP {}  GQ {}".format(
                r.get("gt", "."), r.get("dp", "."), r.get("gq", ".")))
            if notas and gen in notas:
                L.append("     relevancia: {}".format(notas[gen]))
            L.append("")
        return L

    L = ["=" * 78,
         " 07 - HALLAZGOS SECUNDARIOS CANDIDATOS",
         "=" * 78,
         "",
         "Gen causal primario (excluido de esta lista): {}".format(primario),
         "Variantes raras evaluadas: {:,}".format(len(rows)),
         "",
         "CRITERIO: impacto HIGH, o MODERATE con CADD >= {}".format(CADD_MIN)]
    L += bloque("A) Genes de la lista ACMG SF v3.2 (consenso internacional)", acmg)
    L += bloque("B) Genes con enfermedad tratable (criterio propio, no consenso)",
                trat, TRATABLES)
    L += ["",
          "=" * 78,
          " LIMITE DECLARADO",
          "=" * 78,
          "",
          "  Esto NO es un reporte clinico. Es una lista de candidatos para",
          "  revision humana. Un hallazgo secundario reportable exige:",
          "",
          "    - confirmacion ortogonal de la variante",
          "    - clasificacion ACMG completa por un profesional",
          "    - correlacion con historia clinica y familiar",
          "    - consejo genetico antes de cualquier comunicacion",
          "",
          "  Se reporta como hipotesis para seguimiento, nunca como diagnostico.",
          "",
          "  Fuente de la lista A: ACMG SF v3.2, Miller et al., Genet Med 2023.",
          "  La lista B es criterio propio y se declara como tal.",
          "",
          "  -> {}".format(OUT_TSV)]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
