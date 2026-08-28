#!/usr/bin/env python3
"""
==============================================================================
06_validate_convergence.py - Validation gate for the blind search

INPUT  : $RESULTS/05_ranked_genes.tsv
         $WORK/04_rare_candidates.tsv
OUTPUT : $RESULTS/06_convergence_report.txt

PURPOSE
-------
Stages 01-05 run WITHOUT knowing the disease. This is the only stage that
applies external knowledge, and it does so AFTER the ranking is closed: it
checks whether the top-ranked gene supports a biologically plausible compound
heterozygous pair.

This is a VALIDATION GATE, not a search filter. The distinction matters: if
disease knowledge entered before the ranking, the result would be circular and
would prove nothing about the method.

CONVERGENCE CRITERIA
--------------------
1. Separation: the top-1 score exceeds top-2 by a clear margin
2. Model:      top-1 is AR_COMPOUND_HET or AR_HOMOZYGOUS
3. Pair:       it retains >= 2 rare variants after the frequency filter
4. Severity:   at least one HIGH-impact variant
5. Evidence:   independent clinical (ClinVar) or computational (CADD) support

Failing any of the five does not invalidate the finding, but must be declared in
the methods report.
==============================================================================
"""
import csv
import os
import sys

BASE = os.path.expanduser("~/mva")
RANK = BASE + "/results/05_ranked_genes.tsv"
RARE = BASE + "/work/04_rare_candidates.tsv"
OUT = BASE + "/results/06_convergence_report.txt"

MIN_MARGIN = 0.15   # top-1 must lead top-2 by 15 %
HIGH_CADD = 20.0


def fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main():
    for p in (RANK, RARE):
        if not os.path.exists(p):
            sys.exit("missing " + p)

    ranked = list(csv.DictReader(open(RANK), delimiter="\t"))
    if not ranked:
        sys.exit("the ranking is empty")
    rare = list(csv.DictReader(open(RARE), delimiter="\t"))

    top = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    gene = top["gene"]
    vs = [r for r in rare if r["gene"] == gene]

    s1 = fnum(top["pheno_score"])
    s2 = fnum(second["pheno_score"]) if second else 0.0
    margin = (s1 - s2) / s1 if s1 else 0.0

    highs = [v for v in vs if v.get("vep_impact") == "HIGH" or v.get("impact") == "HIGH"]
    clinvar = [v for v in vs if v.get("clinvar", ".") not in (".", "", None)
               and "pathogenic" in v.get("clinvar", "").lower()]
    cadds = [v for v in vs if fnum(v.get("cadd")) >= HIGH_CADD]

    criteria = [
        ("Separation over 2nd", margin >= MIN_MARGIN,
         "{:.1f}% (minimum {:.0f}%)".format(margin * 100, MIN_MARGIN * 100)),
        ("Recessive model", top["model"] in ("AR_COMPOUND_HET", "AR_HOMOZYGOUS"),
         top["model"]),
        ("Pair of rare variants", len(vs) >= 2, "{} rare variants".format(len(vs))),
        ("At least one HIGH", len(highs) >= 1, "{} of HIGH impact".format(len(highs))),
        ("Independent evidence", len(clinvar) >= 1 or len(cadds) >= 1,
         "{} pathogenic in ClinVar, {} with CADD>={}".format(len(clinvar), len(cadds), HIGH_CADD)),
    ]
    passed = sum(1 for _, ok, _ in criteria if ok)

    L = []
    L.append("=" * 78)
    L.append(" 06 - VALIDATION GATE FOR THE BLIND SEARCH")
    L.append("=" * 78)
    L.append("")
    L.append("Stages 01-05 ran without knowing the disease.")
    L.append("This stage evaluates the result AFTER the ranking was closed.")
    L.append("")
    L.append("## Convergent gene")
    L.append("")
    L.append("  TOP-1 : {}   score {}".format(gene, top["pheno_score"]))
    if second:
        L.append("  TOP-2 : {}   score {}".format(second["gene"], second["pheno_score"]))
        L.append("  margin: {:.1f}%".format(margin * 100))
    L.append("")
    L.append("## Rare variants in {}".format(gene))
    L.append("")
    for v in sorted(vs, key=lambda x: int(x["pos"])):
        af = v.get("gnomad_af") or "ABSENT"
        L.append("  {}:{} {}>{}".format(v["chrom"], v["pos"], v["ref"], v["alt"]))
        L.append("     {}  |  impact {}  |  CADD {}".format(
            v.get("vep_consequence", "."), v.get("vep_impact", "."), v.get("cadd") or "."))
        L.append("     {}".format(v.get("vep_hgvsp", ".")))
        L.append("     gnomAD {}  |  {}  |  ClinVar: {}".format(
            af, v.get("rsid", "."), v.get("clinvar", ".")))
        L.append("     genotype {}  DP {}  GQ {}  AD {}".format(
            v.get("gt", "."), v.get("dp", "."), v.get("gq", "."), v.get("ad", ".")))
        L.append("")
    L.append("## Convergence criteria")
    L.append("")
    for name, ok, detail in criteria:
        L.append("  [{}] {:<26} {}".format("OK" if ok else "--", name, detail))
    L.append("")
    L.append("  {} of {} criteria met".format(passed, len(criteria)))
    L.append("")
    if passed == len(criteria):
        L.append("  >>> CONVERGENCE CONFIRMED <<<")
        L.append("")
        L.append("  Without knowing the disease, the pipeline ranked {} first".format(gene))
        L.append("  among {} candidate genes, starting from 5,012,204 variants.".format(len(ranked)))
    else:
        L.append("  >>> PARTIAL CONVERGENCE - declare this in the methods report <<<")
    L.append("")
    L.append("## Limitations that still apply")
    L.append("")
    L.append("  - Singleton: without parents, phase cannot be proven from pedigree.")
    L.append("    The pair is reported as PRESUMED in trans.")
    L.append("  - SNVs and indels only: the VCF has no CNVs, SVs or repeat expansions.")
    L.append("  - The impact filter discards deep intronic and regulatory variants.")
    L.append("  - A gene with no HPO annotations scores 0 even if it were causal.")

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
