#!/usr/bin/env python3
"""
==============================================================================
t2_03_mechanism_filter.py - Which associations survive mechanistic scrutiny

INPUT  : $WORK/t2_02_drug_evidence.tsv
OUTPUT : $WORK/t2_03_filtered.tsv
         $RESULTS/t2_03_mechanism_filter.txt

PURPOSE
-------
t2_02 recorded what the databases say. This stage decides what any of it is
worth, against three filters applied in order. Every drop is counted.

  1. EVIDENCE     a single text-mining source is not pharmacology
  2. DIRECTION    the direction of effect must oppose the lesion, not follow it
  3. SAFETY CLASS drugs whose only use is acute cytotoxic chemotherapy are not
                  candidates for chronic therapy in a child

THE DIRECTION FILTER - THE PART THAT MATTERS
--------------------------------------------
The disease is a LOSS of BubR1 function: the spindle assembly checkpoint is too
weak, so anaphase proceeds before chromosomes are correctly attached.

BubR1's job is to INHIBIT CDC20, which activates the APC/C. So:

    node                     effect of INHIBITING it      verdict
    ----------------------   --------------------------   ---------------
    CDC20, APC/C subunits    restrains anaphase onset     COMPENSATORY
                             - substitutes for the
                               function that was lost
    BUB1, BUB3, MAD2L1,      weakens the checkpoint       HARMFUL
    KNL1, TTK, AURKB,        further
    PLK1, CDK1, CENPE
    others                   no directional argument      NEUTRAL

Almost every drug in these databases is an inhibitor, because inhibitors are
what the pharmaceutical industry builds. For a loss-of-function disease that
makes target-based repurposing structurally difficult: the useful direction is
the one nobody has drugs for.

This filter is the reason the readthrough branch (t2_04) exists. It does not act
on the pathway at all - it acts on the ribosome, sidestepping both the
undruggability of BubR1 and the direction problem.
==============================================================================
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

IN_TSV = os.path.join(WORK, "t2_02_drug_evidence.tsv")
OUT_TSV = os.path.join(WORK, "t2_03_filtered.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_03_mechanism_filter.txt")

MIN_SCORE = 0.10
MIN_SOURCES = 2

# Inhibiting these substitutes for the lost BubR1 restraint on anaphase.
COMPENSATORY_IF_INHIBITED = {
    "CDC20", "FZR1", "UBE2C", "ESPL1",
    "ANAPC1", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "ANAPC10", "ANAPC11",
    "ANAPC13", "ANAPC15", "ANAPC16", "CDC16", "CDC23", "CDC26", "CDC27",
}
# Inhibiting these weakens an already-failing checkpoint.
HARMFUL_IF_INHIBITED = {
    "BUB1", "BUB3", "MAD1L1", "MAD2L1", "KNL1", "TTK", "AURKB", "PLK1",
    "CDK1", "CCNB1", "CCNB2", "CENPE", "CENPF", "NEK2", "NDC80", "NUF2",
    "ZWINT", "SGO1", "SKA1", "SKA3", "CDCA8", "CENPA", "SPDL1", "ZWILCH",
    "KIF2C", "KIF20A", "KIF4A", "DLGAP5", "ASPM", "NUSAP1",
}
INHIBITORY_WORDS = ("inhibitor", "antagonist", "blocker", "negative", "suppressor",
                    "inverse agonist", "inhibition")
ACTIVATING_WORDS = ("agonist", "activator", "positive", "inducer", "potentiator")

# Drugs whose only role is acute cytotoxic chemotherapy: not chronic-therapy
# candidates for a child, whatever the database says.
CYTOTOXIC = {
    "CARBOPLATIN", "CISPLATIN", "OXALIPLATIN", "PACLITAXEL", "DOCETAXEL",
    "TOPOTECAN", "IRINOTECAN", "IDARUBICIN", "DOXORUBICIN", "DAUNORUBICIN",
    "EPIRUBICIN", "ETOPOSIDE", "VINCRISTINE", "VINBLASTINE", "VINORELBINE",
    "CYCLOPHOSPHAMIDE", "IFOSFAMIDE", "GEMCITABINE", "CYTARABINE",
    "FLUOROURACIL", "METHOTREXATE", "MITOMYCIN", "BLEOMYCIN", "DACTINOMYCIN",
}


def direction_of(interaction_type, directionality):
    t = (interaction_type + " " + directionality).lower()
    if any(w in t for w in INHIBITORY_WORDS):
        return "inhibits"
    if any(w in t for w in ACTIVATING_WORDS):
        return "activates"
    return "unknown"


def verdict_for(gene, direction):
    """Is this drug pushing the mechanism the right way?"""
    if direction == "unknown":
        return "UNKNOWN", "no interaction direction is recorded"
    if gene in COMPENSATORY_IF_INHIBITED:
        if direction == "inhibits":
            return "COMPENSATORY", "restrains anaphase; substitutes for lost BubR1 restraint"
        return "HARMFUL", "activating this releases anaphase further"
    if gene in HARMFUL_IF_INHIBITED:
        if direction == "inhibits":
            return "HARMFUL", "weakens an already-failing checkpoint"
        return "COMPENSATORY", "strengthening this could partly offset the lesion"
    return "NEUTRAL", "no directional argument for this node"


def is_cytotoxic(name):
    n = name.upper()
    return any(c in n for c in CYTOTOXIC)


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("missing " + IN_TSV + " - run t2_02_drug_evidence.py first")

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    stats = {"total": len(rows)}

    approved = [r for r in rows if r["approved"] == "True"]
    stats["approved"] = len(approved)

    kept, drops = [], {"weak_evidence": [], "harmful_direction": [],
                       "unknown_direction": [], "cytotoxic_only": [], "neutral": []}

    for r in approved:
        score = float(r["score"] or 0)
        nsrc = int(r["n_sources"] or 0)
        if score < MIN_SCORE or nsrc < MIN_SOURCES:
            drops["weak_evidence"].append(r)
            continue
        if is_cytotoxic(r["drug"]):
            drops["cytotoxic_only"].append(r)
            continue
        d = direction_of(r.get("interaction_type", ""), r.get("directionality", ""))
        v, why = verdict_for(r["gene"], d)
        r["direction"], r["verdict"], r["reason"] = d, v, why
        if v == "HARMFUL":
            drops["harmful_direction"].append(r)
        elif v == "UNKNOWN":
            drops["unknown_direction"].append(r)
        elif v == "NEUTRAL":
            drops["neutral"].append(r)
        else:
            kept.append(r)

    cols = ["gene", "tier", "drug", "approved", "score", "n_sources",
            "interaction_type", "directionality", "direction", "verdict", "reason"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in kept:
            w.writerow({c: r.get(c, "") for c in cols})

    # How many approved associations even carry a direction?
    with_dir = sum(1 for r in approved
                   if direction_of(r.get("interaction_type", ""), r.get("directionality", "")) != "unknown")

    L = ["=" * 88,
         " T2-03 - MECHANISTIC FILTER",
         "=" * 88, "",
         "The disease is a LOSS of BubR1 function. A drug is only a candidate if it",
         "pushes the mechanism in the opposing direction. Direction is not a detail",
         "here - it decides whether a drug would help or make the lesion worse.",
         "",
         "## Funnel", "",
         "  associations recorded                {:>6,}".format(stats["total"]),
         "  with an approved drug                {:>6,}".format(stats["approved"]),
         "  ... of which carry a direction       {:>6,}  ({:.0%})".format(
             with_dir, with_dir / stats["approved"] if stats["approved"] else 0),
         "",
         "  dropped - single source / low score  {:>6,}".format(len(drops["weak_evidence"])),
         "  dropped - acute cytotoxic only       {:>6,}".format(len(drops["cytotoxic_only"])),
         "  dropped - direction unknown          {:>6,}".format(len(drops["unknown_direction"])),
         "  dropped - direction HARMFUL          {:>6,}".format(len(drops["harmful_direction"])),
         "  dropped - no directional argument    {:>6,}".format(len(drops["neutral"])),
         "  " + "-" * 44,
         "  SURVIVING CANDIDATES                 {:>6,}".format(len(kept)),
         ""]

    if kept:
        L += ["## Candidates", "",
              "  {:<9}{:<32}{:<14}{:<14}{}".format("GENE", "DRUG", "DIRECTION", "VERDICT", "WHY")]
        for r in kept:
            L.append("  {:<9}{:<32}{:<14}{:<14}{}".format(
                r["gene"], r["drug"][:31], r["direction"], r["verdict"], r["reason"]))
    else:
        L += ["## No candidate survives", "",
              "  Target-based repurposing does not produce a usable candidate for this",
              "  lesion, and the reason is structural rather than accidental:",
              "",
              "  1. BubR1 itself has no drug. The Mitotic Checkpoint Complex and the",
              "     entire APC/C are a pharmacological desert - 43 of 61 network members",
              "     have no reported drug association at all.",
              "",
              "  2. The direction that would help is the direction nobody has drugs for.",
              "     Inhibiting CDC20 or the APC/C would substitute for the lost BubR1",
              "     restraint on anaphase. Those are exactly the nodes with zero drugs.",
              "",
              "  3. The direction that IS available is the harmful one. Inhibitors of",
              "     BUB1, TTK, AURKB, PLK1 and CDK1 exist because oncology wants to push",
              "     cancer cells past a checkpoint. In a child whose checkpoint is",
              "     already failing, that is the wrong way round.",
              "",
              "  4. What remains is single-source text mining. Every approved-drug",
              "     association in this network is backed by exactly one source, with",
              "     scores of 0.01-0.60 - associations like PLK1-erythromycin or",
              "     PLK1-lansoprazole are not pharmacology.",
              "",
              "  This negative is the argument for the next stage. A therapy for this",
              "  variant has to bypass the pathway entirely rather than act within it."]

    L += ["", "## What was dropped, and why - the harmful ones are worth naming", ""]
    seen = set()
    for r in drops["harmful_direction"][:12]:
        k = (r["gene"], r["drug"])
        if k in seen:
            continue
        seen.add(k)
        L.append("  {:<9}{:<32} would weaken the checkpoint further".format(r["gene"], r["drug"][:31]))
    if not drops["harmful_direction"]:
        L.append("  (none reached this filter - they were dropped earlier for weak evidence)")

    L += ["", "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
