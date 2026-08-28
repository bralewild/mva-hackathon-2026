#!/usr/bin/env python3
"""
==============================================================================
08_make_submission.py - Generate the Track 1 submission CSV automatically

INPUT  : $RESULTS/05_ranked_genes.tsv
         $WORK/04_rare_candidates.tsv
OUTPUT : $PROJECT/submissions/<username>_<model>.csv

PURPOSE
-------
The methods template asks whether the submission file is the automated output of
the computational approach (preferred) or the product of downstream manual
curation. This stage makes the honest answer "automated": the CSV is derived
deterministically from the pipeline's own results, with no hand editing.

Re-running the pipeline regenerates a byte-identical file.

WHAT IT DOES
------------
1. Takes the top-ranked gene from stage 05.
2. Collects its rare variants from stage 04, ordered by position.
3. Emits one row per gene: a pair when two variants support a compound
   heterozygous model, a single variant otherwise.
4. Assigns EPCR from the evidence, not by hand (see epcr_for below).
5. Adds lower-ranked genes as additional rows up to MAX_ROWS, so the ranking the
   evaluator scores is the ranking the pipeline produced.

CHROMOSOME NAMING
-----------------
The source VCF uses Ensembl contig naming (1..22,X,Y). The evaluator performs no
normalisation and expects the UCSC form (chr15), confirmed by the field
description in submit_track1.py, the fallback in groundtruth.py and the official
submission template. SUBMIT_CONTIG_PREFIX is applied here, in one place.
==============================================================================
"""
import csv
import os
import sys

from _paths import WORK, RESULTS, PROJECT

RANK = os.path.join(RESULTS, "05_ranked_genes.tsv")
RARE = os.path.join(WORK, "04_rare_candidates.tsv")

USERNAME = os.environ.get("MVA_USERNAME", "bralewild")
MODEL = os.environ.get("MVA_MODEL", "blind-wgs-triage")
OUT_DIR = os.path.join(PROJECT, "submissions")
OUT_CSV = os.path.join(OUT_DIR, "{}_{}.csv".format(USERNAME, MODEL))

CONTIG_PREFIX = os.environ.get("MVA_CONTIG_PREFIX", "chr")
MAX_ROWS = 10
# Only genes this close to the leader are worth submitting at all. Below it the
# phenotype signal is noise and extra rows only dilute precision.
MIN_SCORE_FRACTION = 0.60

COLS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
        "chrom_2", "pos_2", "ref_2", "alt_2",
        "epcr", "finding_type", "notes"]
PROBAND_ID = os.environ.get("MVA_PROBAND_ID", "PROBAND01")


def fnum(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def epcr_for(rank_index, gene_row, variants):
    """
    Derive the estimated probability of causal relationship from evidence,
    never by hand. Only the ordering matters for scoring, but the value should
    still be defensible.

    Base by rank, adjusted by the strength of the supporting evidence:
      +  a HIGH-impact (truncating) variant
      +  a ClinVar pathogenic classification
      +  a complete compound heterozygous pair
    """
    base = [0.60, 0.30, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04][min(rank_index, 9)]
    bonus = 0.0
    impacts = [(v.get("vep_impact") or v.get("impact") or "") for v in variants]
    if "HIGH" in impacts:
        bonus += 0.15
    if any("pathogenic" in (v.get("clinvar") or "").lower() for v in variants):
        bonus += 0.15
    if len(variants) >= 2 and gene_row.get("model") == "AR_COMPOUND_HET":
        bonus += 0.05
    return round(min(base + bonus, 0.99), 2)


def note_for(gene_row, variants):
    """A compact, factual rationale. No claim the evidence does not support."""
    bits = ["{} - {} model".format(gene_row["gene"], gene_row.get("model", "?").replace("_", " ").lower())]
    for v in variants:
        hgvsp = (v.get("vep_hgvsp") or ".").split(":")[-1]
        af = v.get("gnomad_af")
        bits.append("{} {} (CADD {}, gnomAD {}{})".format(
            hgvsp,
            (v.get("vep_consequence") or ".").split(",")[0],
            v.get("cadd") or "n/a",
            af if af else "absent",
            ", ClinVar " + v["clinvar"] if v.get("clinvar", ".") not in (".", "") else ""))
    bits.append("phenotype rank {} of {} by HPO semantic similarity".format(
        gene_row["_rank"], gene_row["_total"]))
    if len(variants) >= 2:
        bits.append("presumed in trans; phase not provable in a singleton")
    return ". ".join(bits) + "."


def main():
    for p in (RANK, RARE):
        if not os.path.exists(p):
            sys.exit("missing " + p + " - run the pipeline first")

    ranked = list(csv.DictReader(open(RANK), delimiter="\t"))
    rare = list(csv.DictReader(open(RARE), delimiter="\t"))
    if not ranked:
        sys.exit("the ranking is empty")

    top_score = fnum(ranked[0]["pheno_score"])
    cutoff = top_score * MIN_SCORE_FRACTION

    by_gene = {}
    for r in rare:
        by_gene.setdefault(r["gene"], []).append(r)

    rows, used = [], 0
    for i, g in enumerate(ranked):
        if used >= MAX_ROWS:
            break
        if fnum(g["pheno_score"]) < cutoff:
            break
        vs = sorted(by_gene.get(g["gene"], []), key=lambda x: int(x["pos"]))
        if not vs:
            continue
        # A compound heterozygous claim needs exactly a pair of heterozygotes.
        pair = vs[:2] if (g.get("model") == "AR_COMPOUND_HET" and len(vs) >= 2) else vs[:1]
        g = dict(g, _rank=i + 1, _total=len(ranked))
        row = {
            "proband_id": PROBAND_ID,
            "chrom_1": CONTIG_PREFIX + pair[0]["chrom"], "pos_1": pair[0]["pos"],
            "ref_1": pair[0]["ref"], "alt_1": pair[0]["alt"],
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": epcr_for(i, g, pair),
            "finding_type": "primary",
            "notes": note_for(g, pair),
        }
        if len(pair) == 2:
            row.update({"chrom_2": CONTIG_PREFIX + pair[1]["chrom"], "pos_2": pair[1]["pos"],
                        "ref_2": pair[1]["ref"], "alt_2": pair[1]["alt"]})
        rows.append(row)
        used += 1

    # Required by the challenge: rows sorted by EPCR descending.
    rows.sort(key=lambda r: -float(r["epcr"]))

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)

    print("=" * 78)
    print(" 08 - SUBMISSION FILE (automated output, no manual editing)")
    print("=" * 78)
    print()
    print("  rows written : {} (max {})".format(len(rows), MAX_ROWS))
    print("  score cutoff : {:.4f}  ({:.0f}% of the leader's {:.4f})".format(
        cutoff, MIN_SCORE_FRACTION * 100, top_score))
    print()
    for i, r in enumerate(rows, 1):
        v2 = "  +  {}:{} {}>{}".format(r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"]) if r["chrom_2"] else ""
        print("  {}. epcr {}  {}:{} {}>{}{}".format(
            i, r["epcr"], r["chrom_1"], r["pos_1"], r["ref_1"], r["alt_1"], v2))
    print()
    print("  -> " + OUT_CSV)


if __name__ == "__main__":
    main()
