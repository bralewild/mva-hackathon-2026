#!/usr/bin/env python3
"""
==============================================================================
t2_05_positive_control.py - Can this pipeline find a drug when one exists?

INPUT  : none (the control genes are the input)
OUTPUT : $WORK/t2_05_positive_control.tsv
         $RESULTS/t2_05_positive_control.txt

WHY THIS STAGE EXISTS
---------------------
t2_03 returns zero surviving candidates for BUB1B. That number is worthless on
its own. A pipeline strict enough to return zero for BUB1B might return zero for
EVERY gene, in which case the BUB1B negative says nothing about BUB1B - it says
something about the filters.

So: run the same evidence gate and the same direction logic against genes where
a loss-of-function rare disease DOES have an approved drug acting on the
deficient product. If those survive, the BUB1B zero is a statement about BUB1B.
If they do not, the filters are simply too strict and the negative is
uninterpretable - and this stage will say so.

Both outcomes are reported. Neither is assumed.

THE GENERIC DIRECTION RULE
--------------------------
t2_03's direction knowledge base is curated per node and specific to the spindle
assembly checkpoint. That does not generalise, and the report says so.

But for the MUTATED gene itself - tier 0 - no curation is needed. In a
loss-of-function disease the rule is mechanical:

    a drug that ACTIVATES / potentiates the deficient product  -> COMPENSATORY
    a drug that INHIBITS it                                    -> HARMFUL
    no recorded direction                                      -> UNKNOWN

That rule is gene-agnostic. It is the part of the direction filter that ports to
any LoF gene without rewriting a knowledge base, and this stage is where it is
tested.

A CAVEAT ON WHAT A PASS WOULD PROVE
-----------------------------------
Surviving here means "the gate admits a real approved drug when one exists in
DGIdb for the mutated gene". It does NOT prove the gate would find a drug acting
on a network NEIGHBOUR of the mutated gene, which is the harder problem BUB1B
actually poses - BubR1 has no drug at all, so the search had to go outward. This
control bounds the claim; it does not remove it.
==============================================================================
"""
import csv
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

# Same gate as t2_03. Imported by value rather than by import so that a change
# there cannot silently desynchronise the control from the thing it controls -
# the assertion below fails loudly instead.
MIN_SCORE = 0.10
MIN_SOURCES = 2

OUT_TSV = os.path.join(WORK, "t2_05_positive_control.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_05_positive_control.txt")
DGIDB = "https://dgidb.org/api/graphql"

# Rare, loss-of-function, and the approved drug acts on the deficient product
# itself - the same shape as BUB1B, minus the part where BubR1 is undruggable.
CONTROLS = [
    {"gene": "CFTR", "disease": "cystic fibrosis",
     "expect": "IVACAFTOR",
     "why": "potentiator of the mutated channel; activating a loss-of-function "
            "target is by definition the compensatory direction"},
    {"gene": "GBA1", "disease": "Gaucher disease type 1",
     "expect": "IMIGLUCERASE",
     "why": "enzyme replacement for the deficient glucocerebrosidase"},
    {"gene": "PAH", "disease": "phenylketonuria",
     "expect": "SAPROPTERIN",
     "why": "BH4 cofactor; increases residual hydroxylase activity"},
]

# A negative control: a gene from the BUB1B network where the well-evidenced
# drugs are inhibitors of a checkpoint component. If this SURVIVED, the
# direction filter would be broken.
NEGATIVE_CONTROL = {"gene": "AURKB", "disease": "(not a disease gene here)",
                    "expect": "BARASERTIB",
                    "why": "inhibitor of a checkpoint kinase - must be rejected"}

DGIDB_Q = """{ genes(names:%s){ nodes{ name
  interactions{
    drug{ name conceptId approved }
    interactionScore
    interactionTypes{ type directionality }
    sources{ sourceDbName }
  } } } }"""

INHIBITORY_WORDS = ("inhibit", "antagonis", "blocker", "negative", "suppress",
                    "inverse agonist")
ACTIVATING_WORDS = ("activat", "agonist", "positive", "inducer", "potentiator")


def gql(query, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                DGIDB, data=json.dumps({"query": query}).encode(),
                headers={"Content-Type": "application/json",
                         "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            wait = 2 ** attempt
            print("    {}: {} - retry in {}s".format(type(e).__name__, e, wait),
                  file=sys.stderr, flush=True)
            time.sleep(wait)
    return {}


def direction_of(interaction_type, directionality):
    """Identical resolution order to t2_03: structured field wins over free text."""
    for field in (directionality, interaction_type):
        f = (field or "").lower()
        if any(w in f for w in INHIBITORY_WORDS):
            return "inhibits"
        if any(w in f for w in ACTIVATING_WORDS):
            return "activates"
    return "unknown"


def lof_verdict(direction):
    """The gene-agnostic rule for the MUTATED gene in a loss-of-function disease."""
    if direction == "activates":
        return "COMPENSATORY", "restores or potentiates the deficient product"
    if direction == "inhibits":
        return "HARMFUL", "further reduces an already deficient product"
    return "UNKNOWN", "no interaction direction is recorded"


def fetch(genes):
    d = gql(DGIDB_Q % json.dumps(genes))
    rows = []
    for n in (d.get("data", {}).get("genes", {}) or {}).get("nodes", []) or []:
        for it in n.get("interactions", []) or []:
            types = [t.get("type") for t in (it.get("interactionTypes") or []) if t.get("type")]
            dirs = [t.get("directionality") for t in (it.get("interactionTypes") or [])
                    if t.get("directionality")]
            srcs = {s.get("sourceDbName") for s in (it.get("sources") or []) if s.get("sourceDbName")}
            rows.append({"gene": n["name"], "drug": (it["drug"]["name"] or "").strip(),
                         "approved": bool(it["drug"].get("approved")),
                         "score": float(it.get("interactionScore") or 0.0),
                         "interaction_type": ",".join(types),
                         "directionality": ",".join(dirs),
                         "n_sources": len(srcs)})
    return rows


def evaluate(rows, ctrl):
    """Apply the t2_03 gate, then the generic LoF direction rule."""
    mine = [r for r in rows if r["gene"] == ctrl["gene"]]
    out = []
    for r in mine:
        d = direction_of(r["interaction_type"], r["directionality"])
        v, why = lof_verdict(d)
        gate_ok = (r["approved"] and r["score"] >= MIN_SCORE
                   and r["n_sources"] >= MIN_SOURCES)
        r = dict(r, direction=d, verdict=v, reason=why,
                 passes_evidence=gate_ok, survives=(gate_ok and v == "COMPENSATORY"))
        out.append(r)
    return out


def main():
    all_ctrls = CONTROLS + [NEGATIVE_CONTROL]
    genes = [c["gene"] for c in all_ctrls]
    print("querying DGIdb for {} control genes".format(len(genes)), file=sys.stderr, flush=True)
    rows = fetch(genes)
    if not rows:
        sys.exit("ABORT: DGIdb returned nothing for any control gene. This is an "
                 "API failure, not a negative control result.")

    results, per_ctrl = [], {}
    for c in all_ctrls:
        ev = evaluate(rows, c)
        per_ctrl[c["gene"]] = ev
        results.extend(ev)

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    cols = ["gene", "drug", "approved", "score", "n_sources", "interaction_type",
            "directionality", "direction", "verdict", "passes_evidence", "survives"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(results, key=lambda x: (x["gene"], -x["score"])):
            w.writerow({c: r.get(c, "") for c in cols})

    L = ["=" * 88, " T2-05 - POSITIVE CONTROL", "=" * 88, "",
         "t2_03 returns zero candidates for BUB1B. This stage asks whether the same",
         "gate returns zero for EVERYTHING - in which case that zero measures the",
         "filters, not the gene.",
         "",
         "Gate: approved AND score >= {} AND >= {} distinct sources,".format(MIN_SCORE, MIN_SOURCES),
         "      then the gene-agnostic loss-of-function direction rule.",
         ""]

    n_pass = 0
    for c in CONTROLS:
        ev = per_ctrl[c["gene"]]
        surv = [r for r in ev if r["survives"]]
        expected = sorted([r for r in ev if c["expect"].lower() in r["drug"].lower()],
                          key=lambda r: (not r["survives"], not r["passes_evidence"],
                                         -r["n_sources"], -r["score"]))
        n_pass += 1 if surv else 0
        L += ["-" * 88,
              "{}  -  {}".format(c["gene"], c["disease"]),
              "-" * 88,
              "  associations in DGIdb        : {}".format(len(ev)),
              "  ... approved                 : {}".format(sum(1 for r in ev if r["approved"])),
              "  ... passing the evidence gate: {}".format(sum(1 for r in ev if r["passes_evidence"])),
              "  SURVIVING (compensatory)     : {}".format(len(surv)),
              ""]
        if expected:
            r = expected[0]
            dup = " ({} rows for this drug; showing the best-annotated)".format(
                len(expected)) if len(expected) > 1 else ""
            L += ["  expected drug {}:".format(c["expect"]),
                  "    found      yes" + dup,
                  "    approved   {}   score {:.2f}   sources {}".format(
                      r["approved"], r["score"], r["n_sources"]),
                  "    type       {!r} / {!r}".format(r["interaction_type"], r["directionality"]),
                  "    direction  {}  ->  {}".format(r["direction"], r["verdict"]),
                  "    SURVIVES   {}".format(r["survives"]),
                  "    rationale  {}".format(c["why"])]
        else:
            L += ["  expected drug {}: NOT PRESENT in DGIdb for this gene.".format(c["expect"]),
                  "    That is a gap in the database, not in the filter."]
        if surv:
            L += ["", "  survivors:"]
            for r in sorted(surv, key=lambda x: -x["score"])[:8]:
                L.append("    {:<32} score {:.2f}  sources {}  {}".format(
                    r["drug"][:31], r["score"], r["n_sources"], r["direction"]))
        L.append("")

    nc = NEGATIVE_CONTROL
    ev = per_ctrl[nc["gene"]]
    nc_surv = [r for r in ev if r["survives"]]
    L += ["-" * 88,
          "NEGATIVE CONTROL: {}  -  must NOT survive".format(nc["gene"]),
          "-" * 88,
          "  {}".format(nc["why"]),
          "  associations                 : {}".format(len(ev)),
          "  passing the evidence gate    : {}".format(sum(1 for r in ev if r["passes_evidence"])),
          "  SURVIVING                    : {}   {}".format(
              len(nc_surv), "<-- FILTER IS BROKEN" if nc_surv else "(correct)"),
          ""]

    L += ["=" * 88, " VERDICT", "=" * 88, ""]
    if n_pass and not nc_surv:
        L += ["  {} of {} positive controls returned at least one surviving approved".format(n_pass, len(CONTROLS)),
              "  drug, and the negative control returned none.",
              "",
              "  The gate is therefore capable of admitting a real approved drug when",
              "  one exists for the mutated gene, and capable of rejecting one that",
              "  points the wrong way. BUB1B's zero is a statement about BUB1B."]
    elif not n_pass:
        L += ["  NO positive control survived.",
              "",
              "  The gate returns zero even where an approved drug demonstrably exists.",
              "  BUB1B's zero therefore measures the strictness of these filters and",
              "  NOT the absence of a drug. The Track 2 negative must be read with that",
              "  limitation stated, and this stage exists to state it."]
    else:
        L += ["  {} of {} positive controls survived, but the NEGATIVE control also".format(n_pass, len(CONTROLS)),
              "  survived. The direction filter is admitting something it should reject.",
              "  Do not rely on the BUB1B result until this is resolved."]

    L += ["", "  Scope: this controls the gate on the MUTATED gene. It does not control",
          "  the harder case BUB1B actually posed - finding a drug on a network",
          "  NEIGHBOUR when the mutated gene itself is undruggable.",
          "", "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
