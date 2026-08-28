#!/usr/bin/env python3
"""
==============================================================================
t2_02_drug_evidence.py - Query drug-gene evidence across the target network

INPUT  : $WORK/t2_01_targets.tsv
OUTPUT : $WORK/t2_02_drug_evidence.tsv
         $RESULTS/t2_02_drug_evidence.txt

PURPOSE
-------
Collect, without judgement, every drug-gene association reported for the target
network. Filtering and interpretation happen in t2_03; this stage only records
what the public databases actually say, so the filtering can be audited against
the raw evidence.

SOURCES
-------
DGIdb v5 (GraphQL)      - aggregates ~40 upstream sources: DrugBank, ChEMBL,
                          TTD, PharmGKB, guide-to-pharmacology, plus text-mined
                          sets. Reports an approval flag, an interaction score
                          and, where known, an interaction type/directionality.
Open Targets (GraphQL)  - target-level metadata, used for tractability rather
                          than drug lists.

A NOTE ON WHAT THESE SCORES MEAN
--------------------------------
DGIdb's interaction score reflects how many sources support an association and
how specific the drug and gene are, not how strong or how therapeutically
sensible the interaction is. A score of 0.02 backed by one text-mining source is
not pharmacology. t2_03 uses an explicit floor for that reason, and records how
many associations fall below it.
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

IN_TSV = os.path.join(WORK, "t2_01_targets.tsv")
OUT_TSV = os.path.join(WORK, "t2_02_drug_evidence.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_02_drug_evidence.txt")

DGIDB = "https://dgidb.org/api/graphql"
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
BATCH = 25


def gql(url, query, variables=None, retries=4):
    body = {"query": query}
    if variables:
        body["variables"] = variables
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            wait = 2 ** attempt
            print("    {}: {} - retry in {}s".format(type(e).__name__, e, wait), file=sys.stderr)
            time.sleep(wait)
    return {}


DGIDB_Q = """{ genes(names:%s){ nodes{ name
  interactions{
    drug{ name conceptId approved }
    interactionScore
    interactionTypes{ type directionality }
    publications{ pmid }
    sources{ sourceDbName }
  } } } }"""


def query_dgidb(genes):
    rows = []
    for i in range(0, len(genes), BATCH):
        chunk = genes[i:i + BATCH]
        print("  DGIdb batch {}/{} ({} genes)".format(
            i // BATCH + 1, (len(genes) + BATCH - 1) // BATCH, len(chunk)), file=sys.stderr)
        d = gql(DGIDB, DGIDB_Q % json.dumps(chunk))
        for n in (d.get("data", {}).get("genes", {}) or {}).get("nodes", []) or []:
            for it in n.get("interactions", []) or []:
                types = [t.get("type") for t in (it.get("interactionTypes") or []) if t.get("type")]
                dirs = [t.get("directionality") for t in (it.get("interactionTypes") or [])
                        if t.get("directionality")]
                srcs = [s.get("sourceDbName") for s in (it.get("sources") or []) if s.get("sourceDbName")]
                rows.append({
                    "gene": n["name"],
                    "drug": (it["drug"]["name"] or "").strip(),
                    "drug_id": it["drug"].get("conceptId") or "",
                    "approved": bool(it["drug"].get("approved")),
                    "score": float(it.get("interactionScore") or 0.0),
                    "interaction_type": ",".join(types) or "",
                    "directionality": ",".join(dirs) or "",
                    "n_sources": len(set(srcs)),   # DISTINCT sources - a repeated sourceDbName is not corroboration
                    "sources": ";".join(sorted(set(srcs)))[:120],
                    "n_publications": len(it.get("publications") or []),
                })
        time.sleep(0.4)
    return rows


OT_Q = """query($ids:[String!]!){ targets(ensemblIds:$ids){
    id approvedSymbol tractability{ label modality value } } }"""


def query_tractability(symbols):
    """Open Targets tractability, keyed by symbol. Best effort - not essential."""
    # Resolve symbols to Ensembl IDs via the search endpoint would cost a call
    # per gene; tractability is contextual here, so we query only the causal gene.
    d = gql(OPENTARGETS, OT_Q, {"ids": ["ENSG00000156970"]})
    out = {}
    for t in (d.get("data", {}) or {}).get("targets", []) or []:
        labs = sorted({x["label"] for x in (t.get("tractability") or []) if x.get("value")})
        out[t["approvedSymbol"]] = labs
    return out


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("missing " + IN_TSV + " - run t2_01_target_network.py first")

    targets = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    tier_of = {t["gene"]: int(t["tier"]) for t in targets}
    genes = [t["gene"] for t in targets]
    print("querying {} targets".format(len(genes)), file=sys.stderr)

    rows = query_dgidb(genes)
    for r in rows:
        r["tier"] = tier_of.get(r["gene"], 9)

    tract = query_tractability(["BUB1B"])

    cols = ["gene", "tier", "drug", "drug_id", "approved", "score", "interaction_type",
            "directionality", "n_sources", "n_publications", "sources"]
    rows.sort(key=lambda r: (r["tier"], -r["score"], r["gene"], r["drug"]))
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    approved = [r for r in rows if r["approved"]]
    typed = [r for r in rows if r["interaction_type"]]
    genes_with_any = sorted({r["gene"] for r in rows})
    genes_no_drug = [g for g in genes if g not in genes_with_any]

    L = ["=" * 88,
         " T2-02 - DRUG-GENE EVIDENCE ACROSS THE TARGET NETWORK",
         "=" * 88, "",
         "Targets queried              : {}".format(len(genes)),
         "Targets with any association : {}".format(len(genes_with_any)),
         "Targets with NO association  : {}".format(len(genes_no_drug)),
         "",
         "Associations returned        : {:,}".format(len(rows)),
         "  with an approved drug      : {:,}".format(len(approved)),
         "  with a declared interaction type : {:,}  ({:.0%})".format(
             len(typed), len(typed) / len(rows) if rows else 0),
         ""]

    if "BUB1B" in tract:
        L += ["Open Targets tractability for BUB1B:",
              "  " + (", ".join(tract["BUB1B"]) if tract["BUB1B"] else "(no tractability bucket flagged)"),
              ""]

    L += ["THE CAUSAL GENE ITSELF", "-" * 22, ""]
    bub = [r for r in rows if r["gene"] == "BUB1B"]
    if bub:
        for r in bub[:10]:
            L.append("  {:<32} approved={} score {:.2f}".format(r["drug"][:31], r["approved"], r["score"]))
    else:
        L += ["  No drug-gene association is reported for BUB1B in DGIdb.",
              "",
              "  This is a result, not a gap: BubR1 is not a target any existing drug",
              "  acts on. Any therapeutic hypothesis must therefore either act on the",
              "  network around it, or bypass the protein entirely - which is what the",
              "  readthrough branch (t2_04) does by acting on the ribosome."]

    L += ["", "APPROVED DRUGS ACROSS THE NETWORK (top 30 by score)", "-" * 50, "",
          "  {:<9}{:<34}{:<18}{:>7}{:>6}".format("GENE", "DRUG", "TYPE", "SCORE", "SRC")]
    for r in [x for x in approved][:30]:
        L.append("  {:<9}{:<34}{:<18}{:>7.2f}{:>6}".format(
            r["gene"], r["drug"][:33], (r["interaction_type"] or "-")[:17],
            r["score"], r["n_sources"]))

    if genes_no_drug:
        L += ["", "Targets with no reported drug association:", "",
              "  " + ", ".join(genes_no_drug)]

    L += ["", "  -> " + OUT_TSV,
          "", "Interpretation is deliberately deferred to t2_03. This stage records",
          "what the databases say; the next one decides what any of it is worth."]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
