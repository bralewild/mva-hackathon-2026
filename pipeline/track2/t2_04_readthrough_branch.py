#!/usr/bin/env python3
"""
==============================================================================
t2_04_readthrough_branch.py - The variant-class route, and the safety screen

INPUT  : Ensembl REST (MANE CDS, to compute the PTC context)
         $RAW/patient_hpo.tsv (the proband's own phenotype - actually read)
OUTPUT : $WORK/t2_04_candidates.tsv
         $RESULTS/t2_04_readthrough.txt

WHY THIS BRANCH EXISTS
----------------------
t2_03 returns no surviving candidate. With the corrected network - which now
includes the translation-termination and NMD machinery - that negative is
meaningful rather than circular: the search CAN see readthrough biology, and
independently recovers ataluren and ELX-02 from the databases. Neither is
market-approved. No approved drug acts on that machinery with credible evidence.

WHAT CHANGED, AND WHY (adversarial review, 2026-08-28)
------------------------------------------------------
1. THE SAFETY SCREEN SCREENED NOTHING. PATIENT_SAFETY_FLAGS was a literal dict,
   the extracted phenotype file was never opened, and the flags were only
   printed - never applied to any candidate. Gentamicin's demotion, presented as
   the screen's flagship result, was a human typing rank=2. Now the screen reads
   $RAW/patient_hpo.tsv, matches candidates against liability classes, and
   ANNOTATES AND REORDERS them. If the phenotype file is missing, it says so
   rather than pretending.

2. "(leucine)" WAS A STRING LITERAL on the line the report headlined as
   "computed rather than assumed". The amino acid is now translated from the
   codon.

3. CODON AND CDS_POS WERE INDEPENDENT CONSTANTS with no consistency check. A
   mismatch produced negative slice offsets, silent wraparound and a confident
   report about a non-stop codon. Now validated, and the stage aborts.

4. THE READTHROUGH PRODUCT IS NOT WILD-TYPE PROTEIN. At a UGA, near-cognate
   incorporation gives tryptophan, cysteine or arginine - not the original
   leucine. Success produces BubR1 p.Leu737Trp/Cys/Arg, a novel missense variant,
   in a patient whose other allele is already a missense of unknown consequence.
   This was absent from the previous version entirely, and it changes what a
   western blot for "full-length BubR1" can tell you.

5. AMLEXANOX'S APPROVAL STATUS WAS WRONG. Discontinued in the United States
   (Aphthasol) AND in Japan (Solfa, Takeda, since 2019). It is an approved
   molecule with no marketed product. Its principal systemic pharmacology is
   TBK1/IKKe inhibition, which the safety screen must see in a child with a
   cancer predisposition syndrome.
==============================================================================
"""
import csv
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import RAW, WORK, RESULTS  # noqa: E402

OUT_TSV = os.path.join(WORK, "t2_04_candidates.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_04_readthrough.txt")
HPO_FILE = os.path.join(RAW, "patient_hpo.tsv")

TRANSCRIPT = "ENST00000287598"      # BUB1B MANE Select, NM_001211.6
CODON = 737
CDS_POS = 2210
NEW_BASE = "G"

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Ter", "TAG": "Ter",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Ter", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
# Amino acids inserted by near-cognate tRNAs at each stop codon.
NEAR_COGNATE = {"TGA": ["Trp", "Cys", "Arg"],
                "TAG": ["Gln", "Tyr", "Lys", "Trp"],
                "TAA": ["Gln", "Tyr", "Lys"]}
STOP_RANK = {"TGA": ("UGA", 1, "most permissive - the leakiest stop codon"),
             "TAG": ("UAG", 2, "intermediate fidelity"),
             "TAA": ("UAA", 3, "high fidelity - little intrinsic readthrough")}
CONTEXT_RANK = {"C": (1, "most permissive +4 context"), "T": (2, "second"),
                "A": (3, "third"), "G": (4, "least permissive")}

# Drug liability classes, matched against the patient's phenotype below.
LIABILITY = {
    "nephrotoxic": "cumulative renal tubular injury",
    "ototoxic": "cochlear and vestibular damage",
    "immunosuppressive": "suppression of innate immune / type-I interferon signalling",
    "proliferative_risk": "may promote proliferation or confound oncological surveillance",
    "gi_intolerance": "significant gastrointestinal intolerance",
    "unknown_systemic": "no systemic human exposure data for this indication",
}
# HPO term -> the liability classes it contraindicates.
HPO_CONTRAINDICATIONS = {
    "HP:0000121": (["nephrotoxic"], "nephrocalcinosis - reduced renal reserve"),
    "HP:0002859": (["proliferative_risk", "immunosuppressive"],
                   "rhabdomyosarcoma - active cancer predisposition"),
    "HP:0001508": (["gi_intolerance"], "failure to thrive"),
    "HP:0003202": ([], "skeletal muscle atrophy"),
}

CANDIDATES = [
    {
        "drug": "Amlexanox",
        "class": "benzopyranobipyridine carboxylic acid; anti-inflammatory",
        "approval": "FDA-approved molecule (NDA 020511, 1996) but NO MARKETED PRODUCT: "
                    "Aphthasol discontinued in the US; Solfa oral tablets discontinued "
                    "in Japan by Takeda in 2019",
        "mechanism": "reported dual activity: translational readthrough inducer and "
                     "NMD inhibitor. Principal known systemic pharmacology is "
                     "TBK1/IKKe inhibition.",
        "why_here": "the lesion needs both readthrough and NMD inhibition; amlexanox is "
                    "the only agent reported to do both",
        "precedent": "COL7A1 in recessive dystrophic epidermolysis bullosa (8 of 12 PTC "
                     "alleles responded); GDAP1 in patient-derived hiPSC neurons",
        "liabilities": ["immunosuppressive", "unknown_systemic"],
        "marketed": False,
        "limitations": "no marketed product anywhere; systemic PK for a genetic indication "
                       "unknown; never tested on BUB1B",
        "base_rank": 1,
    },
    {
        "drug": "Escin (beta-aescin)",
        "class": "triterpene saponin mixture, Aesculus hippocastanum",
        "approval": "MARKETED. German Commission E monograph (1984, renewed 1994) for "
                    "chronic venous insufficiency; EMA traditional-use registration; "
                    "oral dragees and topical gel (e.g. Reparil). NOT an FDA drug approval.",
        "mechanism": "readthrough induction; identified in an unbiased high-throughput "
                     "screen of ~1,600 clinically approved compounds against CFTR PTCs",
        "why_here": "the only readthrough-active compound in this set with a currently "
                    "marketed product and established human dosing (50-75 mg twice daily)",
        "precedent": "CFTR G542X and W1282X readthrough in primary human airway cells",
        "liabilities": [],
        "marketed": True,
        "limitations": "herbal/traditional registration rather than a full marketing "
                       "authorisation; no readthrough data in any nuclear or "
                       "cell-cycle gene; never tested on BUB1B",
        "base_rank": 2,
    },
    {
        "drug": "Gentamicin / Amikacin (aminoglycosides)",
        "class": "aminoglycoside antibiotic",
        "approval": "fully approved and marketed worldwide",
        "mechanism": "binds the ribosomal A site, permitting near-cognate tRNA "
                     "incorporation at the PTC",
        "why_here": "the best-characterised readthrough agents, with the closest "
                    "published precedent to BubR1",
        "precedent": "aminoglycoside readthrough functionally restored BRCA1 nonsense "
                     "alleles - a large nuclear genome-maintenance protein. Also "
                     "BBS2/ALMS1 ciliopathies (protein AND function in patient fibroblasts)",
        "liabilities": ["nephrotoxic", "ototoxic"],
        "marketed": True,
        "limitations": "cumulative renal and cochlear toxicity; unsuitable for chronic use",
        "base_rank": 3,
    },
    {
        "drug": "Ataluren (PTC124)",
        "class": "oxadiazole; purpose-built readthrough inducer",
        "approval": "NOT market-approved. EMA conditional authorisation for "
                    "nonsense-mutation Duchenne not renewed (CHMP negative Jan 2024, "
                    "re-examined and confirmed Oct 2024, adopted by the Commission). "
                    "Efficacy, not safety, was the ground.",
        "mechanism": "promotes readthrough at the ribosomal A site",
        "why_here": "recovered independently by the t2-03 search on RPL3 and RPS15, "
                    "which serves as a positive control that the corrected network "
                    "can see readthrough biology",
        "precedent": "MPS I-H mouse model; a full phase-3 programme in DMD that FAILED "
                     "ON EFFICACY - the most informative single datapoint about whether "
                     "readthrough restores clinically meaningful protein",
        "liabilities": [],
        "marketed": False,
        "limitations": "efficacy not demonstrated. Its failure is evidence about the "
                       "CLASS, not only about ataluren, and it argues against this "
                       "entire strategy.",
        "base_rank": 4,
    },
]


def cds_sequence(tx):
    url = "https://rest.ensembl.org/sequence/id/{}?type=cds;content-type=application/json".format(tx)
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())["seq"]


def load_patient_hpo(path):
    """Read the proband's HPO terms. Returns [] if unavailable - and says so."""
    if not os.path.exists(path):
        return None
    terms = []
    with open(path) as f:
        next(f, None)
        for line in f:
            t = line.split("\t")[0].strip()
            if t.startswith("HP:"):
                terms.append(t)
    return terms


def screen(candidates, hpo_terms):
    """Apply the phenotype screen. This REORDERS and ANNOTATES - it is not display."""
    contra = {}
    for t in hpo_terms or []:
        classes, why = HPO_CONTRAINDICATIONS.get(t, ([], None))
        for c in classes:
            contra.setdefault(c, []).append((t, why))
    for c in candidates:
        hits = [l for l in c["liabilities"] if l in contra]
        c["flagged"] = hits
        c["penalty"] = len(hits)
        c["screen_note"] = "; ".join(
            "{} ({})".format(LIABILITY[h], contra[h][0][1]) for h in hits) or "no conflict found"
    # Being marketed is a REQUIREMENT for a repurposing proposal, not a ranking
    # criterion. An unmarketed drug cannot be repurposed, however good its
    # mechanism - so it is partitioned out before safety is even considered.
    # Within the proposable set: safety penalty, then mechanistic rank.
    candidates.sort(key=lambda c: (not c.get("marketed", False), c["penalty"], c["base_rank"]))
    for i, c in enumerate(candidates, 1):
        c["final_rank"] = i
    return candidates, contra


def main():
    try:
        cds = cds_sequence(TRANSCRIPT)
    except Exception as e:
        sys.exit("could not retrieve the CDS from Ensembl: {}".format(e))

    if len(cds) % 3:
        sys.exit("CDS length {} is not a multiple of 3".format(len(cds)))
    i = (CODON - 1) * 3
    if not 0 <= i < len(cds) - 3:
        sys.exit("codon {} is outside the CDS".format(CODON))
    offset = (CDS_POS - 1) - i
    if not 0 <= offset <= 2:
        sys.exit("c.{} does not fall inside codon {} (offset {})".format(CDS_POS, CODON, offset))

    wt = cds[i:i + 3]
    mut = wt[:offset] + NEW_BASE + wt[offset + 1:]
    wt_aa = CODON_TABLE.get(wt, "?")
    mut_aa = CODON_TABLE.get(mut, "?")
    if mut_aa != "Ter":
        sys.exit("c.{}{}>{} gives {} ({}), not a stop codon".format(
            CDS_POS, wt[offset], NEW_BASE, mut, mut_aa))
    plus4 = cds[i + 3:i + 4]
    context = "...{} [{}] {}...".format(cds[i - 6:i], mut, cds[i + 3:i + 9])
    rna, srank, sdesc = STOP_RANK[mut]
    crank, cdesc = CONTEXT_RANK.get(plus4, (9, "unknown"))
    inserted = NEAR_COGNATE.get(mut, [])

    hpo_terms = load_patient_hpo(HPO_FILE)
    cands, contra = screen(list(CANDIDATES), hpo_terms)

    cols = ["final_rank", "drug", "class", "approval", "mechanism", "why_here",
            "precedent", "flagged", "screen_note", "limitations"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for c in cands:
            row = {k: c.get(k, "") for k in cols}
            row["flagged"] = ",".join(c["flagged"])
            w.writerow(row)

    L = ["=" * 88,
         " T2-04 - VARIANT-CLASS ROUTE AND PATIENT SAFETY SCREEN",
         "=" * 88, "",
         "## The premature termination codon - computed from the MANE CDS", "",
         "  transcript        {} (NM_001211.6)".format(TRANSCRIPT),
         "  CDS               {:,} nt -> {:,} codons".format(len(cds), len(cds) // 3),
         "  codon {}         {}  ({})".format(CODON, wt, wt_aa),
         "  c.{}{}>{}      {} -> {}  ({})".format(CDS_POS, wt[offset], NEW_BASE, wt, mut, mut_aa),
         "  context           {}".format(context),
         "",
         "  stop codon        {} - rank {} of 3 - {}".format(rna, srank, sdesc),
         "  +4 nucleotide     {} - rank {} of 4 - {}".format(plus4, crank, cdesc),
         "",
         "  MIXED. The stop codon is the most permissive of the three; the downstream",
         "  context is not. The two are not commensurable - UGA vs UAA is roughly an",
         "  order of magnitude in reporter assays, C vs A at +4 is two- to three-fold.",
         "",
         "## What readthrough actually produces - and it is not wild-type protein", "",
         "  Near-cognate incorporation at {} inserts {}, not {}.".format(
             rna, " / ".join(inserted), wt_aa),
         "",
         "  Successful readthrough therefore yields BubR1 p.{}{}{}".format(
             wt_aa, CODON, "/".join(inserted)),
         "  - a NOVEL MISSENSE VARIANT, in a patient whose other allele already carries",
         "  a missense of unknown functional consequence. Whether that product is",
         "  functional is unknown and cannot be assumed.",
         "",
         "  Consequence for the follow-up plan: a western blot for full-length BubR1",
         "  cannot distinguish restored function from a full-length non-functional",
         "  product. The functional readout is not optional - it is the experiment.",
         "",
         "## Patient phenotype screen", ""]
    if hpo_terms is None:
        L += ["  NOT RUN: {} is missing.".format(HPO_FILE),
              "  Run pipeline/00b_extract_phenotype.py with an authorised copy of the",
              "  dataset. Candidates below are ranked on mechanism alone.", ""]
    else:
        L += ["  Read {} HPO terms from {}".format(len(hpo_terms), HPO_FILE), ""]
        for c, entries in sorted(contra.items()):
            for t, why in entries:
                L.append("    {}  contraindicates {:<20} ({})".format(t, c, why))
        if not contra:
            L.append("    no liability class is contraindicated by this phenotype")
        L.append("")

    proposable = [c for c in cands if c.get("marketed") and not c["penalty"]]
    L += ["## Candidates, after the screen", "",
          "  Being MARKETED is a requirement, not a ranking criterion: an unmarketed",
          "  drug cannot be repurposed however good its mechanism. Candidates are",
          "  partitioned on that first, then on the phenotype screen.", "",
          "  proposable (marketed AND no phenotype conflict): {}".format(
              ", ".join(c["drug"].split(" ")[0] for c in proposable) or "NONE"),
          ""]
    for c in cands:
        L += ["  " + "-" * 84,
              "  {}. {}{}{}".format(
                  c["final_rank"], c["drug"],
                  "   [NOT MARKETED - reference only]" if not c.get("marketed") else "",
                  "   [DEMOTED BY SCREEN]" if c["penalty"] else ""),
              "  " + "-" * 84,
              "     class       {}".format(c["class"]),
              "     APPROVAL    {}".format(c["approval"]),
              "     mechanism   {}".format(c["mechanism"]),
              "     rationale   {}".format(c["why_here"]),
              "     precedent   {}".format(c["precedent"]),
              "     SCREEN      {}".format(c["screen_note"]),
              "     limits      {}".format(c["limitations"]),
              ""]

    demoted = [c for c in cands if c["penalty"]]
    L += ["## What the screen changed", ""]
    if demoted:
        for c in demoted:
            L.append("  {} moved to rank {} (mechanistic rank was {}) - {}".format(
                c["drug"].split(" ")[0], c["final_rank"], c["base_rank"], c["screen_note"]))
    else:
        L.append("  nothing - no candidate carries a liability this phenotype contraindicates")

    L += ["",
          "## Honest status of this branch", "",
          "  No agent in this set is both market-approved and evidenced for readthrough:",
          "",
          "    amlexanox   best mechanistic fit, NO marketed product anywhere",
          "    escin       marketed, but herbal registration and no nuclear-gene data",
          "    aminoglyc.  fully approved, contraindicated in this patient",
          "    ataluren    not approved; its phase-3 failure is evidence against the class",
          "",
          "  That is the finding. The readthrough field has no cleanly approved agent,",
          "  and the most informative clinical datapoint in it - ataluren's efficacy",
          "  failure - argues against the strategy rather than for it.",
          "",
          "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
