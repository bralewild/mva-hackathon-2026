#!/usr/bin/env python3
"""
==============================================================================
t2_04_readthrough_branch.py - The branch that bypasses the pathway

INPUT  : Ensembl REST (BUB1B MANE CDS, to compute the PTC context)
OUTPUT : $WORK/t2_04_candidates.tsv
         $RESULTS/t2_04_readthrough.txt

WHY THIS BRANCH EXISTS
----------------------
t2_03 returned zero candidates, for structural reasons: BubR1 has no drug, the
Mitotic Checkpoint Complex and APC/C have none either (43 of 61 network members
with no association at all), the one direction that would help is the direction
nobody builds drugs for, and 88 % of the approved-drug associations that remain
do not even record whether the drug inhibits or activates.

Target-based repurposing cannot solve a loss-of-function lesion in an
undruggable complex. So this branch does not act on the pathway. It acts on the
ribosome.

THE ARGUMENT
------------
1. One allele carries a premature termination codon (p.Leu737Ter). Translational
   readthrough drugs let the ribosome insert an amino acid at a PTC and produce
   full-length protein.
2. Aneuploidy appears below ~50 % residual BUB1B expression (Hanks et al.). The
   patient has one null allele and one full-length missense allele, so total
   BubR1 sits near that threshold. The therapeutic goal is not cure - it is to
   cross a threshold.
3. The PTC is UGA, the most readthrough-permissive stop codon. The +1 context is
   less favourable. Both are computed here rather than assumed.
4. The transcript is NMD-degraded, so readthrough alone has little substrate.
   The mechanism therefore needs BOTH readthrough induction and NMD inhibition.

That last requirement is what selects the lead candidate: amlexanox is an
FDA-approved drug reported to do both.

WHAT THIS IS NOT
----------------
A treatment recommendation. These are hypotheses for laboratory follow-up. No
readthrough agent has been tested on BUB1B in any system.
==============================================================================
"""
import csv
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

OUT_TSV = os.path.join(WORK, "t2_04_candidates.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_04_readthrough.txt")

TRANSCRIPT = "ENST00000287598"      # BUB1B MANE Select, NM_001211.6
CODON = 737                          # p.Leu737Ter
CDS_POS = 2210                       # c.2210T>G
NEW_BASE = "G"

# Readthrough permissiveness, from the literature: UGA > UAG >> UAA, and the
# nucleotide immediately 3' of the stop matters in the order C > U > A > G.
STOP_RANK = {"TGA": ("UGA", 1, "most permissive - the leakiest stop codon"),
             "TAG": ("UAG", 2, "intermediate fidelity"),
             "TAA": ("UAA", 3, "high fidelity - little intrinsic readthrough")}
CONTEXT_RANK = {"C": (1, "most permissive +1 context"),
                "T": (2, "second"), "A": (3, "third"), "G": (4, "least permissive")}

# Curated candidate set. Approval status stated explicitly, because the
# challenge asks for market-approved medicines and one obvious candidate is not.
CANDIDATES = [
    {
        "drug": "Amlexanox",
        "class": "benzopyranobipyridine carboxylic acid; anti-inflammatory",
        "approval": "FDA-approved (5% oral paste, aphthous ulcers); oral anti-allergic in Japan",
        "mechanism": "DUAL: translational readthrough inducer (TRID) AND NMD inhibitor",
        "why_here": "the lesion needs both - readthrough has no substrate while NMD "
                    "degrades the transcript, and amlexanox addresses both with one molecule",
        "precedent": "COL7A1 in recessive dystrophic epidermolysis bullosa: 8 of 12 PTC "
                     "alleles responded, some reaching >50% of normal full-length protein. "
                     "GDAP1 in patient-derived hiPSC neurons (Charcot-Marie-Tooth).",
        "fit_to_patient": "no known nephrotoxicity - important because the child has "
                          "nephrocalcinosis (HP:0000121), which contraindicates chronic "
                          "aminoglycoside exposure",
        "limitations": "approved formulation is topical/oral-local; systemic exposure for a "
                       "genetic indication would require reformulation and dose-finding. "
                       "Never tested on BUB1B. Readthrough efficiency is allele-specific.",
        "rank": 1,
    },
    {
        "drug": "Gentamicin / Amikacin (aminoglycosides)",
        "class": "aminoglycoside antibiotic",
        "approval": "approved worldwide, decades of clinical use",
        "mechanism": "binds the ribosomal A site, permitting near-cognate tRNA "
                     "incorporation at the PTC",
        "why_here": "the best-characterised readthrough agents; direct precedent in a "
                    "genome-stability gene",
        "precedent": "aminoglycoside-induced readthrough functionally restored BRCA1 "
                     "nonsense alleles - the closest published analogue to BubR1: a large "
                     "nuclear protein in genome maintenance, rescued by an approved drug. "
                     "Also BBS2/ALMS1 ciliopathies (protein AND function restored in "
                     "patient fibroblasts), COL4A5, fucosidosis.",
        "fit_to_patient": "POOR. Nephrotoxic and ototoxic. The child has nephrocalcinosis, "
                          "so chronic systemic exposure is contraindicated. Retained as a "
                          "mechanistic reference and an ex vivo tool compound, not as a "
                          "therapeutic proposal for this patient.",
        "limitations": "cumulative renal and cochlear toxicity; unsuitable for long-term use "
                       "in any patient, and specifically counter-indicated in this one",
        "rank": 2,
    },
    {
        "drug": "Ataluren (PTC124)",
        "class": "oxadiazole; translational readthrough inducer",
        "approval": "NOT currently market-approved. EMA conditional authorisation for "
                    "nonsense-mutation Duchenne was NOT renewed - CHMP negative in Jan 2024, "
                    "re-examined and confirmed Oct 2024, adopted by the European Commission. "
                    "Individual states may still permit named-patient use.",
        "mechanism": "promotes readthrough at the ribosomal A site with a better safety "
                     "profile than aminoglycosides",
        "why_here": "mechanistically the most relevant purpose-built agent, and the reason "
                    "this class exists clinically at all",
        "precedent": "MPS I-H mouse model; extensive DMD trial programme",
        "fit_to_patient": "not applicable while unapproved; efficacy was the ground for "
                          "non-renewal, not safety",
        "limitations": "efficacy not demonstrated to the EMA's satisfaction. Listed for "
                       "completeness and honesty, NOT proposed as a candidate.",
        "rank": 3,
    },
    {
        "drug": "NMD inhibition as an adjunct (e.g. NMDI-1 class)",
        "class": "investigational NMD pathway inhibitors",
        "approval": "investigational - not approved",
        "mechanism": "raises the abundance of the PTC-containing transcript, increasing the "
                     "substrate available for readthrough",
        "why_here": "co-administration of an NMD inhibitor with gentamicin restored "
                    "full-length protein in a Hurler syndrome model; patients with higher "
                    "transcript levels respond better to readthrough drugs",
        "precedent": "Hurler syndrome model (NMDI-1 + gentamicin); choroideremia, where "
                     "NMD efficiency varies and predicts readthrough response",
        "fit_to_patient": "conceptual - included because it explains why a dual-activity "
                          "molecule is preferable to a pure readthrough agent",
        "limitations": "no approved agent; chronic global NMD inhibition has obvious "
                       "safety concerns, since NMD is a housekeeping quality-control pathway",
        "rank": 4,
    },
]

# The patient's HPO terms that carry drug-safety consequences.
PATIENT_SAFETY_FLAGS = {
    "HP:0000121": ("Nephrocalcinosis",
                   "avoid nephrotoxic agents - aminoglycosides, high-dose NSAIDs, "
                   "and anything requiring renal dose adjustment"),
    "HP:0002859": ("Rhabdomyosarcoma",
                   "avoid agents that could promote proliferation or interfere with "
                   "oncological surveillance"),
    "HP:0001508": ("Failure to thrive",
                   "avoid agents with significant gastrointestinal intolerance or "
                   "appetite suppression"),
}


def cds_sequence(tx):
    url = "https://rest.ensembl.org/sequence/id/{}?type=cds;content-type=application/json".format(tx)
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())["seq"]


def main():
    try:
        cds = cds_sequence(TRANSCRIPT)
    except Exception as e:
        sys.exit("could not retrieve the CDS from Ensembl: {}".format(e))

    i = (CODON - 1) * 3
    wt = cds[i:i + 3]
    offset = (CDS_POS - 1) - i            # which base of the codon is mutated
    mut = wt[:offset] + NEW_BASE + wt[offset + 1:]
    plus1 = cds[i + 3:i + 4]
    context = "...{} [{}] {}...".format(cds[i - 6:i], mut, cds[i + 3:i + 9])

    rna, srank, sdesc = STOP_RANK.get(mut, ("?", 9, "not a stop codon"))
    crank, cdesc = CONTEXT_RANK.get(plus1, (9, "unknown"))

    with open(OUT_TSV, "w", newline="") as f:
        cols = ["rank", "drug", "class", "approval", "mechanism", "why_here",
                "precedent", "fit_to_patient", "limitations"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for c in sorted(CANDIDATES, key=lambda x: x["rank"]):
            w.writerow({k: c[k] for k in cols})

    L = ["=" * 88,
         " T2-04 - READTHROUGH BRANCH: BYPASSING AN UNDRUGGABLE PATHWAY",
         "=" * 88, "",
         "t2_03 returned zero candidates. That negative is the argument for this branch:",
         "a loss-of-function lesion in an undruggable complex cannot be addressed from",
         "within the pathway, so this one acts on the ribosome instead.",
         "",
         "## The premature termination codon, computed not assumed", "",
         "  transcript          {} (MANE Select, NM_001211.6)".format(TRANSCRIPT),
         "  CDS length          {:,} nt -> {:,} codons".format(len(cds), len(cds) // 3),
         "  codon {}          {}  (leucine)".format(CODON, wt),
         "  c.{}{}>{}       {} -> {}".format(CDS_POS, wt[offset], NEW_BASE, wt, mut),
         "  stop generated      {} = {}".format(mut, rna),
         "  sequence context    {}".format(context),
         "",
         "  readthrough permissiveness (UGA > UAG >> UAA)",
         "    stop codon        {}  - rank {} of 3 - {}".format(rna, srank, sdesc),
         "    +1 nucleotide     {}  - rank {} of 4 - {}".format(plus1, crank, cdesc),
         "",
         "  MIXED, and worth stating plainly: the stop codon is the most permissive of",
         "  the three, the downstream context is not. Readthrough efficiency at this",
         "  specific allele is an empirical question, and the first thing any follow-up",
         "  should measure.",
         "",
         "## Why a threshold, not a cure", "",
         "  Aneuploidy appears below ~50% residual BUB1B expression. This patient has:",
         "",
         "    allele 1   p.Leu737Ter    NMD-degraded          -> ~0% contribution",
         "    allele 2   p.Asn1002Lys   full-length, kinase   -> reduced function",
         "                              domain substitution",
         "",
         "  Total BubR1 sits near the threshold. The therapeutic target is not to",
         "  restore normal function but to recover enough of allele 1 to cross it -",
         "  a quantitative, measurable goal.",
         "",
         "## Candidates", ""]

    for c in sorted(CANDIDATES, key=lambda x: x["rank"]):
        L += ["  " + "-" * 84,
              "  {}. {}".format(c["rank"], c["drug"]),
              "  " + "-" * 84,
              "     class       {}".format(c["class"]),
              "     APPROVAL    {}".format(c["approval"]),
              "     mechanism   {}".format(c["mechanism"]),
              "     rationale   {}".format(c["why_here"]),
              "     precedent   {}".format(c["precedent"]),
              "     patient fit {}".format(c["fit_to_patient"]),
              "     limits      {}".format(c["limitations"]),
              ""]

    L += ["## Patient-specific safety screen", "",
          "  Candidate drugs are cross-checked against the proband's own phenotype.",
          "  This is where a generic repurposing answer and a patient-specific one",
          "  diverge:", ""]
    for hpo, (label, consequence) in PATIENT_SAFETY_FLAGS.items():
        L.append("    {}  {:<24} {}".format(hpo, label, consequence))
    L += ["",
          "  The nephrocalcinosis flag is decisive. Aminoglycosides are the",
          "  best-evidenced readthrough agents and are approved worldwide - and they",
          "  are contraindicated for chronic use in this specific child. A proposal",
          "  that ignored that would be unusable clinically, whatever the mechanism",
          "  said.",
          "",
          "## What follow-up would look like", "",
          "  1. Quantify BUB1B transcript in patient cells - readthrough response",
          "     correlates with transcript abundance, and NMD efficiency varies.",
          "  2. Measure baseline BubR1 protein against the ~50% threshold.",
          "  3. Test readthrough ex vivo in patient-derived fibroblasts: full-length",
          "     BubR1 by western blot, then premature chromatid separation and",
          "     aneuploidy rate as functional readouts.",
          "  4. Only if protein AND function move together does this become a clinical",
          "     question rather than a laboratory one.",
          "",
          "  Nothing here is a treatment recommendation. No readthrough agent has been",
          "  tested on BUB1B in any system - that gap is the novelty and the risk at",
          "  the same time.",
          "",
          "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
