#!/usr/bin/env python3
"""
==============================================================================
t2_01_target_network.py - From the variant to a candidate target set

INPUT  : none (the causal gene is the entry point; see BELOW on why that is
         legitimate here while it would not have been in Track 1)
OUTPUT : $WORK/t2_01_targets.tsv
         $RESULTS/t2_01_target_network.txt

WHY THIS STAGE IS NOT "CHEATING" THE WAY TRACK 1 WOULD HAVE BEEN
----------------------------------------------------------------
Track 1 had to be blind: naming the gene up front would have made the result
circular. Track 2 starts *from* the Track 1 answer by design - the challenge
states the task is to go "from variant/mechanism to candidate medication(s)".
The causal gene is an input here, not a leaked answer.

WHAT IT DOES
------------
Builds the neighbourhood in which a drug could plausibly act, with a stated
rationale for every member:

  tier 0  the mutated gene itself
  tier 1  the Mitotic Checkpoint Complex - BubR1's direct functional partners
  tier 2  the APC/C it restrains, and the kinases regulating the checkpoint
  tier 3  high-confidence STRING interactors not already captured

Tier is not importance; it is mechanistic distance from the lesion. A drug
acting further away needs a stronger argument, and the tier records that.

SOURCE: STRING v12 (confidence >= 0.9, physical + functional), plus a curated
core from the SAC literature so the set does not depend on one database being
up on the day.
==============================================================================
"""
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

CAUSAL_GENE = "BUB1B"
STRING_MIN_SCORE = 900          # 0.9 confidence
OUT_TSV = os.path.join(WORK, "t2_01_targets.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_01_target_network.txt")

# Curated core, so the network does not depend on a single API being available.
CURATED = {
    "BUB1B": (0, "the mutated gene; BubR1, MCC component and CDC20 inhibitor"),
    "BUB1":  (1, "MCC; kinetochore kinase, recruits BubR1"),
    "BUB3":  (1, "MCC; scaffolds BUB1/BUB1B at kinetochores"),
    "MAD2L1": (1, "MCC; co-inhibitor of CDC20 with BubR1"),
    "CDC20": (1, "the direct target of MCC inhibition; APC/C activator"),
    "MAD1L1": (1, "MAD2 template at unattached kinetochores"),
    "KNL1":  (1, "kinetochore scaffold recruiting BUB1/BUB3"),
    "CENPE": (2, "kinetochore motor; activates BubR1 kinase activity"),
    "CENPF": (2, "kinetochore component, BubR1-associated"),
    "TTK":   (2, "MPS1 kinase, apex of checkpoint signalling"),
    "AURKB": (2, "corrects erroneous attachments; checkpoint kinase"),
    "PLK1":  (2, "mitotic kinase, kinetochore-microtubule attachment"),
    "CDK1":  (2, "master mitotic kinase driving entry into mitosis"),
    "CCNB1": (2, "cyclin B1, CDK1 partner, APC/C substrate"),
    "NEK2":  (2, "centrosome separation, mitotic fidelity"),
    "ANAPC1": (2, "APC/C subunit"), "ANAPC2": (2, "APC/C subunit"),
    "ANAPC4": (2, "APC/C subunit"), "ANAPC7": (2, "APC/C subunit"),
    "ANAPC10": (2, "APC/C subunit"), "CDC16": (2, "APC/C subunit"),
    "CDC23": (2, "APC/C subunit"), "CDC27": (2, "APC/C subunit"),
}


def string_partners(gene, min_score):
    """High-confidence STRING partners. Returns {} if STRING is unreachable."""
    url = ("https://string-db.org/api/json/network"
           "?identifiers={}&species=9606&required_score={}&limit=60".format(gene, min_score))
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.loads(r.read().decode())
    except Exception as e:
        print("  STRING unavailable ({}); using the curated core only".format(e),
              file=sys.stderr)
        return {}
    out = {}
    for e in data:
        for side in ("preferredName_A", "preferredName_B"):
            name = e.get(side)
            if name and name != gene:
                out[name] = max(out.get(name, 0.0), float(e.get("score", 0)))
    return out


def main():
    partners = string_partners(CAUSAL_GENE, STRING_MIN_SCORE)

    targets = {}
    for g, (tier, why) in CURATED.items():
        targets[g] = {"gene": g, "tier": tier, "rationale": why,
                      "string_score": round(partners.get(g, 0.0), 3),
                      "source": "curated" + ("+STRING" if g in partners else "")}
    for g, sc in partners.items():
        if g not in targets:
            targets[g] = {"gene": g, "tier": 3,
                          "rationale": "high-confidence STRING interactor of {}".format(CAUSAL_GENE),
                          "string_score": round(sc, 3), "source": "STRING"}

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
         "Entry point : {} (the Track 1 result - an input here, not a leaked answer)".format(CAUSAL_GENE),
         "STRING      : confidence >= {:.1f}, {} partners retrieved".format(
             STRING_MIN_SCORE / 1000, len(partners)),
         "Targets     : {}".format(len(rows)),
         "",
         "Tier is mechanistic DISTANCE from the lesion, not importance. A drug",
         "acting further from the lesion needs a stronger argument to be credible.",
         ""]
    names = {0: "tier 0 - the mutated gene",
             1: "tier 1 - the Mitotic Checkpoint Complex",
             2: "tier 2 - APC/C and checkpoint kinases",
             3: "tier 3 - other high-confidence interactors"}
    for t in sorted(by_tier):
        L += ["", names.get(t, "tier {}".format(t)), "-" * len(names.get(t, "")), ""]
        for r in by_tier[t]:
            sc = "{:.3f}".format(r["string_score"]) if r["string_score"] else "  -  "
            L.append("  {:<9} STRING {}   {}".format(r["gene"], sc, r["rationale"]))
    L += ["", "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
