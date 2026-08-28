#!/usr/bin/env python3
"""
==============================================================================
07_secondary_findings.py - Clinically actionable secondary findings

INPUT  : $WORK/04_rare_candidates.tsv
         $RESULTS/05_ranked_genes.tsv
OUTPUT : $RESULTS/07_secondary_findings.txt
         $WORK/07_secondary_candidates.tsv

PURPOSE
-------
The challenge FAQ states that secondary findings do not affect the automated
score and are set aside for qualitative review by the judging panel. It is a
risk-free additional evaluation route.

Independently of the competition, reporting actionable secondary findings is
standard clinical practice: if a pathogenic variant appears in a gene with an
available medical intervention, it should be flagged even when it does not
explain the presenting phenotype.

CRITERIA
--------
A variant is flagged as a candidate secondary finding when it:
  1. lies in an ACMG SF v3.2 gene (international consensus list), or in a gene
     from the curated treatable-disease list below
  2. is NOT the primary causal gene (top-1 of the phenotype ranking)
  3. has HIGH impact, or MODERATE impact with CADD >= 20
  4. is rare (already guaranteed by stage 04)

IMPORTANT - DECLARED LIMIT
--------------------------
This is NOT a clinical report. It is a list of candidates for human review. A
reportable secondary finding requires orthogonal confirmation, full ACMG
classification and genetic counselling. Results are hypotheses for follow-up,
never diagnoses.

SOURCE of list A: ACMG SF v3.2 (Miller et al., Genet Med 2023). Embedded rather
than downloaded so the pipeline reproduces without depending on a live URL.
==============================================================================
"""
import csv
import os
import sys

from _paths import BASE  # honours MVA_BASE; see pipeline/_paths.py
RARE = BASE + "/work/04_rare_candidates.tsv"
RANK = BASE + "/results/05_ranked_genes.tsv"
OUT_TSV = BASE + "/work/07_secondary_candidates.tsv"
OUT_TXT = BASE + "/results/07_secondary_findings.txt"

MIN_CADD = 20.0

# ACMG SF v3.2 - genes with reportable secondary findings (Miller et al. 2023)
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
    # Metabolic and other
    "OTC", "ATP7B", "BTD", "GAA", "APOB", "LDLR", "PCSK9", "RYR1", "CACNA1S",
    "ACVRL1", "ENG", "RPE65", "TTR", "GCK", "HNF1A", "HNF1B", "HNF4A",
}

# Genes outside ACMG SF but with a treatable or manageable disease.
# Reported separately and with less weight: this is our own judgement, not
# international consensus, and it is declared as such.
TREATABLE = {
    "SERPINA1": "alpha-1 antitrypsin deficiency - pulmonary/hepatic management, augmentation therapy",
    "CBS": "homocystinuria - responds to pyridoxine/betaine and a methionine-restricted diet",
    "PAH": "phenylketonuria - dietary management",
    "GALT": "galactosemia - dietary management",
    "SLC22A5": "primary carnitine deficiency - L-carnitine supplementation",
    "ATM": "ataxia-telangiectasia / carrier cancer risk - surveillance",
    "MEFV": "familial Mediterranean fever - colchicine",
    "G6PD": "G6PD deficiency - avoid oxidative drugs",
    "ACADM": "MCAD deficiency - avoid prolonged fasting",
    "F5": "factor V Leiden - thrombotic risk management",
    "F2": "prothrombin G20210A - thrombotic risk management",
}


def fnum(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def main():
    for p in (RARE, RANK):
        if not os.path.exists(p):
            sys.exit("missing " + p)

    ranked = list(csv.DictReader(open(RANK), delimiter="\t"))
    primary = ranked[0]["gene"] if ranked else None

    rows = list(csv.DictReader(open(RARE), delimiter="\t"))

    acmg, treat = [], []
    for r in rows:
        gene = r.get("gene", "")
        if not gene or gene == primary:
            continue
        impact = r.get("vep_impact") or r.get("impact") or "."
        cadd = fnum(r.get("cadd"))
        if not (impact == "HIGH" or (impact == "MODERATE" and cadd >= MIN_CADD)):
            continue
        if gene in ACMG_SF:
            acmg.append(r)
        elif gene in TREATABLE:
            treat.append(r)

    def order(r):
        return (-(1 if (r.get("vep_impact") or r.get("impact")) == "HIGH" else 0),
                -fnum(r.get("cadd")))

    acmg.sort(key=order)
    treat.sort(key=order)

    if acmg or treat:
        cols = list((acmg or treat)[0].keys())
        with open(OUT_TSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in acmg + treat:
                w.writerow(r)

    def block(title, items, notes=None):
        L = ["", title, "-" * len(title), ""]
        if not items:
            L.append("  (none)")
            return L
        for r in items:
            gene = r["gene"]
            af = r.get("gnomad_af") or "ABSENT"
            L.append("  {}  {}:{} {}>{}".format(gene, r["chrom"], r["pos"], r["ref"], r["alt"]))
            L.append("     {}  |  impact {}  |  CADD {}".format(
                r.get("vep_consequence", "."),
                r.get("vep_impact") or r.get("impact") or ".",
                r.get("cadd") or "."))
            L.append("     {}".format(r.get("vep_hgvsp", ".")))
            L.append("     gnomAD {}  |  {}  |  ClinVar: {}".format(
                af, r.get("rsid", "."), r.get("clinvar", ".")))
            L.append("     genotype {}  DP {}  GQ {}".format(
                r.get("gt", "."), r.get("dp", "."), r.get("gq", ".")))
            if notes and gene in notes:
                L.append("     relevance: {}".format(notes[gene]))
            L.append("")
        return L

    L = ["=" * 78,
         " 07 - CANDIDATE SECONDARY FINDINGS",
         "=" * 78,
         "",
         "Primary causal gene (excluded from this list): {}".format(primary),
         "Rare variants evaluated: {:,}".format(len(rows)),
         "",
         "CRITERION: HIGH impact, or MODERATE with CADD >= {}".format(MIN_CADD)]
    L += block("A) ACMG SF v3.2 genes (international consensus)", acmg)
    L += block("B) Genes with treatable disease (own judgement, not consensus)",
               treat, TREATABLE)
    L += ["",
          "=" * 78,
          " DECLARED LIMIT",
          "=" * 78,
          "",
          "  This is NOT a clinical report. It is a list of candidates for human",
          "  review. A reportable secondary finding requires:",
          "",
          "    - orthogonal confirmation of the variant",
          "    - full ACMG classification by a professional",
          "    - correlation with clinical and family history",
          "    - genetic counselling before any communication",
          "",
          "  Reported as hypotheses for follow-up, never as a diagnosis.",
          "",
          "  Source of list A: ACMG SF v3.2, Miller et al., Genet Med 2023.",
          "  List B is our own judgement and is declared as such.",
          "",
          "  -> {}".format(OUT_TSV)]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
