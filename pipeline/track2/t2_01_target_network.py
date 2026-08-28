#!/usr/bin/env python3
"""
==============================================================================
t2_01_target_network.py - From the variant to a candidate target set

INPUT  : none (the causal gene is the entry point; see BELOW)
OUTPUT : $WORK/t2_01_targets.tsv
         $RESULTS/t2_01_target_network.txt

WHY THIS STAGE IS NOT "CHEATING" THE WAY TRACK 1 WOULD HAVE BEEN
----------------------------------------------------------------
Track 1 had to be blind: naming the gene up front would have made the result
circular. Track 2 starts *from* the Track 1 answer by design - the challenge
task is to go "from variant/mechanism to candidate medication(s)". The causal
gene is an input here, not a leaked answer.

WHAT CHANGED, AND WHY (adversarial review, 2026-08-28)
------------------------------------------------------
The first version took the top 60 STRING partners at combined score >= 0.9 and
called that the mechanistic neighbourhood. Two defects, both material:

1. THE COMBINED SCORE SELECTS CO-EXPRESSION, NOT MECHANISM. Mitotic genes are
   famously co-expressed and co-cited, so STRING's combined score saturates on
   those channels. Verification against the run showed the network had pulled in
   BRCA2, TOP2A and BIRC5 - and those genes supplied the only well-evidenced
   drugs in the whole set (olaparib, etoposide, trastuzumab). PARP inhibitors and
   topoisomerase poisons have nothing to do with restoring BubR1 dosage. They
   inflated the "43 of 61 targets have no drug" statistic with genes that were
   never candidates.
   FIX: physical-interaction subnetwork only, and the curated core carries the
   mechanism rather than the API.

2. THE SEARCH WAS BLIND TO ITS OWN ANSWER. The proposal that eventually emerged
   acts on the ribosome - translational readthrough plus NMD inhibition. No node
   of the translation-termination or NMD machinery was ever in the network, so
   the search could not have found it. A negative from a search that cannot see
   the class of answer being proposed measures its own blind spot.
   FIX: tier 4 adds the variant-class-actionable machinery. A nonsense allele is
   actionable at the ribosome, so the ribosome belongs in the target set.

Tier is mechanistic DISTANCE from the lesion, not importance. Tier 4 is not
"further away" - it is a different route to the same lesion, reachable because
of the variant CLASS rather than the gene.

SOURCES: curated core from the SAC and translation-termination literature;
STRING v12 physical subnetwork for additional partners.
==============================================================================
"""
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

CAUSAL_GENE = "BUB1B"
VARIANT_CLASS = "nonsense"          # drives tier 4 inclusion
STRING_MIN_SCORE = 900
STRING_NETWORK_TYPE = "physical"    # NOT the combined score - see docstring
STRING_LIMIT = 40

OUT_TSV = os.path.join(WORK, "t2_01_targets.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_01_target_network.txt")

# Tiers 0-2: the lesion and what it mechanistically controls.
CURATED = {
    "BUB1B": (0, "the mutated gene; BubR1, MCC component and CDC20 inhibitor"),
    "BUB1":  (1, "MCC; kinetochore kinase, recruits BubR1"),
    "BUB3":  (1, "MCC; scaffolds BUB1/BUB1B at kinetochores"),
    "MAD2L1": (1, "MCC; co-inhibitor of CDC20 with BubR1"),
    "CDC20": (1, "the direct target of MCC inhibition; APC/C activator"),
    "MAD1L1": (1, "MAD2 template at unattached kinetochores"),
    "KNL1":  (1, "kinetochore scaffold recruiting BUB1/BUB3"),
    "CENPE": (2, "kinetochore motor; BubR1-associated"),
    "TTK":   (2, "MPS1 kinase, apex of checkpoint signalling"),
    "AURKB": (2, "corrects erroneous attachments; checkpoint kinase"),
    "PLK1":  (2, "mitotic kinase, kinetochore-microtubule attachment"),
    "CDK1":  (2, "master mitotic kinase driving entry into mitosis"),
    "CCNB1": (2, "cyclin B1, CDK1 partner, APC/C substrate"),
    "ESPL1": (2, "separase; executes sister chromatid separation"),
    "PTTG1": (2, "securin; inhibits separase until anaphase"),
    "ANAPC1": (2, "APC/C subunit"), "ANAPC2": (2, "APC/C subunit"),
    "ANAPC4": (2, "APC/C subunit"), "ANAPC7": (2, "APC/C subunit"),
    "ANAPC10": (2, "APC/C subunit"), "ANAPC11": (2, "APC/C subunit"),
    "CDC16": (2, "APC/C subunit"), "CDC23": (2, "APC/C subunit"),
    "CDC27": (2, "APC/C subunit"), "CDC26": (2, "APC/C subunit"),
    "FZR1":  (2, "Cdh1; the G1 APC/C activator (not the anaphase one)"),
    "UBE2C": (2, "APC/C E2 enzyme"),
    "PPP2CA": (2, "PP2A catalytic; opposes Aurora B at kinetochores"),
    "PPP2R5A": (2, "PP2A-B56; recruited by the BubR1 KARD motif"),
}

# Tier 4: actionable because of the VARIANT CLASS, not the gene. A premature
# termination codon is addressable at the ribosome and through NMD, entirely
# outside the protein's own pathway. Omitting this tier is what made the first
# version's negative uninformative.
VARIANT_CLASS_TARGETS = {
    "nonsense": {
        "ETF1":  "eRF1; recognises the stop codon - the direct target of readthrough",
        "GSPT1": "eRF3a; termination GTPase, target of molecular-glue degraders",
        "GSPT2": "eRF3b paralogue",
        "UPF1":  "NMD helicase; degrades the PTC transcript that readthrough needs",
        "UPF2":  "NMD core factor",
        "UPF3B": "NMD core factor",
        "SMG1":  "NMD kinase, phosphorylates UPF1",
        "SMG5":  "NMD; UPF1 dephosphorylation",
        "SMG6":  "NMD endonuclease",
        "SMG7":  "NMD; recruits deadenylation machinery",
        "EIF4A3": "exon junction complex; marks the PTC as premature",
        "RPL3":  "large ribosomal subunit; readthrough modulation",
        "RPS15": "small subunit near the decoding site",
    },
}


def string_partners(gene, min_score, network_type, limit):
    """Physical-interaction partners. Returns {} if STRING is unreachable."""
    url = ("https://string-db.org/api/json/network"
           "?identifiers={}&species=9606&required_score={}"
           "&network_type={}&limit={}".format(gene, min_score, network_type, limit))
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.loads(r.read().decode())
    except Exception as e:
        print("  STRING unavailable ({}); curated core only".format(e), file=sys.stderr)
        return {}
    out = {}
    for e in data:
        for side in ("preferredName_A", "preferredName_B"):
            name = e.get(side)
            if name and name != gene:
                out[name] = max(out.get(name, 0.0), float(e.get("score", 0)))
    return out


def main():
    partners = string_partners(CAUSAL_GENE, STRING_MIN_SCORE, STRING_NETWORK_TYPE, STRING_LIMIT)
    truncated = len(partners) >= STRING_LIMIT

    targets = {}
    for g, (tier, why) in CURATED.items():
        targets[g] = {"gene": g, "tier": tier, "rationale": why,
                      "string_score": round(partners.get(g, 0.0), 3),
                      "source": "curated" + ("+STRING" if g in partners else "")}
    for g, sc in partners.items():
        if g not in targets:
            targets[g] = {"gene": g, "tier": 3,
                          "rationale": "physical STRING interactor of {}".format(CAUSAL_GENE),
                          "string_score": round(sc, 3), "source": "STRING-physical"}
    for g, why in VARIANT_CLASS_TARGETS.get(VARIANT_CLASS, {}).items():
        if g not in targets:
            targets[g] = {"gene": g, "tier": 4, "rationale": why,
                          "string_score": round(partners.get(g, 0.0), 3),
                          "source": "variant-class"}

    rows = sorted(targets.values(), key=lambda r: (r["tier"], -r["string_score"], r["gene"]))

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w") as f:
        f.write("gene\ttier\tstring_score\tsource\trationale\n")
        for r in rows:
            f.write("{gene}\t{tier}\t{string_score}\t{source}\t{rationale}\n".format(**r))

    by_tier = {}
    for r in rows:
        by_tier.setdefault(r["tier"], []).append(r)

    L = ["=" * 88,
         " T2-01 - CANDIDATE TARGET NETWORK",
         "=" * 88, "",
         "Entry point   : {} (the Track 1 result - an input, not a leaked answer)".format(CAUSAL_GENE),
         "Variant class : {} - drives tier 4".format(VARIANT_CLASS),
         "STRING        : {} subnetwork, score >= {:.1f}, limit {} -> {} partners{}".format(
             STRING_NETWORK_TYPE, STRING_MIN_SCORE / 1000, STRING_LIMIT, len(partners),
             "  ** TRUNCATED AT LIMIT **" if truncated else ""),
         "Targets       : {}".format(len(rows)),
         "",
         "Tier is mechanistic DISTANCE from the lesion, not importance.",
         "Tier 4 is not further away - it is a different route to the same lesion,",
         "reachable because of the variant CLASS rather than the gene.",
         ""]
    names = {0: "tier 0 - the mutated gene",
             1: "tier 1 - the Mitotic Checkpoint Complex",
             2: "tier 2 - APC/C, separase and checkpoint kinases",
             3: "tier 3 - additional physical interactors",
             4: "tier 4 - actionable by VARIANT CLASS (translation termination and NMD)"}
    for t in sorted(by_tier):
        head = names.get(t, "tier {}".format(t))
        L += ["", head, "-" * len(head), ""]
        for r in by_tier[t]:
            sc = "{:.3f}".format(r["string_score"]) if r["string_score"] else "  -  "
            L.append("  {:<9} STRING {}   {}".format(r["gene"], sc, r["rationale"]))

    L += ["",
          "## Note on network construction", "",
          "  The physical subnetwork is used deliberately. STRING's COMBINED score",
          "  saturates on co-expression and text-mining for mitotic genes, and in the",
          "  first version of this pipeline it pulled BRCA2, TOP2A and BIRC5 into a",
          "  spindle-checkpoint drug search. Those genes then supplied the only",
          "  well-evidenced drugs in the entire result - olaparib, etoposide,",
          "  trastuzumab - none of which has anything to do with BubR1 dosage, and all",
          "  of which distorted the summary statistics.",
          "",
          "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
