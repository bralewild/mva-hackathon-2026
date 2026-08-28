#!/usr/bin/env python3
"""
==============================================================================
05_phenotype_rank.py - Rank candidate genes by phenotype similarity

INPUT  : $WORK/04_rare_candidates.tsv
         $RAW/patient_hpo.tsv          (the patient's HPO terms)
OUTPUT : $RESULTS/05_ranked_genes.tsv
         $RESULTS/05_ranking_report.txt

PURPOSE
-------
Final step of the blind search. Up to here the pipeline reduced 5,012,204
variants to a set of genes carrying rare, potentially damaging variants, with no
knowledge of the disease. Those genes are now ordered by how closely their known
phenotype resembles the patient's.

METHOD
------
Resnik semantic similarity over the HPO ontology:

  IC(t)      = -log( genes annotated to t or its descendants / total genes )
  sim(a, b)  = max IC over the common ancestors of a and b
  score(gene)= mean over the patient's terms of
               max_b sim(patient_term, b)  for b in the gene's terms

The MEAN (not the sum) is used so that genes with many annotations are not
rewarded for volume. This is the principle behind Phenomizer and the Exomiser
prioritiser, implemented explicitly and auditably.

SOURCES (public, not patient data -> stored in $ANNOT):
  hp.obo                  https://purl.obolibrary.org/obo/hp.obo
  genes_to_phenotype.txt  https://purl.obolibrary.org/obo/hp/hpoa/
==============================================================================
"""
import collections
import csv
import math
import os
import sys
import urllib.request

BASE = os.path.expanduser("~/mva")
ANNOT = BASE + "/data/annot"
IN_TSV = BASE + "/work/04_rare_candidates.tsv"
HPO_PT = BASE + "/data/raw/patient_hpo.tsv"
OUT_TSV = BASE + "/results/05_ranked_genes.tsv"
OUT_TXT = BASE + "/results/05_ranking_report.txt"

SRC = {
    "hp.obo": "https://purl.obolibrary.org/obo/hp.obo",
    "genes_to_phenotype.txt": "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt",
}


def fetch(name, url):
    """Download a public resource into $ANNOT unless it is already there."""
    path = os.path.join(ANNOT, name)
    if os.path.exists(path) and os.path.getsize(path) > 100000:
        return path
    os.makedirs(ANNOT, exist_ok=True)
    print("  downloading " + name + " ...", file=sys.stderr)
    urllib.request.urlretrieve(url, path)
    return path


def parse_obo(path):
    """Return (parents, names). Obsolete terms are discarded."""
    parents = collections.defaultdict(set)
    names = {}
    tid = None
    obsolete = False
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            tid, obsolete = None, False
        elif line.startswith("id: HP:"):
            tid = line[4:].strip()
        elif line.startswith("name:") and tid:
            names[tid] = line[5:].strip()
        elif line.startswith("is_obsolete: true"):
            obsolete = True
            if tid:
                parents.pop(tid, None)
                names.pop(tid, None)
        elif line.startswith("is_a: HP:") and tid and not obsolete:
            parents[tid].add(line[6:].split("!")[0].strip())
    return parents, names


def ancestors(term, parents, cache):
    """Transitive closure of ancestors, including the term itself."""
    if term in cache:
        return cache[term]
    out = {term}
    stack = list(parents.get(term, ()))
    while stack:
        p = stack.pop()
        if p not in out:
            out.add(p)
            stack.extend(parents.get(p, ()))
    cache[term] = out
    return out


def load_gene_annotations(path):
    """gene -> set of HPO terms, detecting columns from the header."""
    gene_terms = collections.defaultdict(set)
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
        low = [c.strip().lower() for c in hdr]
        i_sym = next((i for i, c in enumerate(low) if "symbol" in c), 1)
        i_hpo = next((i for i, c in enumerate(low) if c in ("hpo_id", "hpo-id")), 2)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(i_sym, i_hpo) and p[i_hpo].startswith("HP:"):
                gene_terms[p[i_sym]].add(p[i_hpo])
    return gene_terms


def load_patient_terms(path):
    terms = []
    with open(path) as f:
        next(f, None)
        for line in f:
            t = line.split("\t")[0].strip()
            if t.startswith("HP:"):
                terms.append(t)
    return terms


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("missing " + IN_TSV + " - run 04_frequency_clinical.py first")
    if not os.path.exists(HPO_PT):
        sys.exit("missing " + HPO_PT + " - run 00b_extract_phenotype.py first")

    obo = fetch("hp.obo", SRC["hp.obo"])
    g2p = fetch("genes_to_phenotype.txt", SRC["genes_to_phenotype.txt"])

    parents, names = parse_obo(obo)
    cache = {}

    gene_terms = load_gene_annotations(g2p)
    print("  genes with HPO annotations: {:,}".format(len(gene_terms)), file=sys.stderr)

    # Information content, propagating each annotation to its ancestors
    freq = collections.Counter()
    for terms in gene_terms.values():
        anc = set()
        for t in terms:
            anc |= ancestors(t, parents, cache)
        for a in anc:
            freq[a] += 1
    total_genes = max(len(gene_terms), 1)
    IC = {t: -math.log(c / total_genes) for t, c in freq.items() if c > 0}

    def resnik(a, b):
        common = ancestors(a, parents, cache) & ancestors(b, parents, cache)
        return max((IC.get(c, 0.0) for c in common), default=0.0)

    patient = load_patient_terms(HPO_PT)
    print("  patient terms: {}".format(len(patient)), file=sys.stderr)

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    by_gene = collections.defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    print("  candidate genes: {:,}".format(len(by_gene)), file=sys.stderr)

    def cadd_of(v):
        try:
            return float(v.get("cadd") or 0)
        except (TypeError, ValueError):
            return 0.0

    scored = []
    for gene, variants in by_gene.items():
        gt = gene_terms.get(gene, set())
        if gt:
            sims = [max((resnik(pt, b) for b in gt), default=0.0) for pt in patient]
            score = sum(sims) / len(patient) if patient else 0.0
            matched = sum(1 for s in sims if s > 0)
        else:
            score, matched = 0.0, 0
        best = max(variants, key=lambda v: (v.get("impact") == "HIGH", cadd_of(v)))
        scored.append({
            "gene": gene,
            "pheno_score": round(score, 4),
            "hpo_terms_matched": matched,
            "n_hpo_annotated": len(gt),
            "model": best.get("model", "."),
            "n_variants": len(variants),
            "best_impact": best.get("impact", "."),
            "best_effect": best.get("effect", "."),
            "best_variant": "{}:{}{}>{}".format(best["chrom"], best["pos"], best["ref"], best["alt"]),
            "best_hgvsp": best.get("vep_hgvsp", "."),
            "best_cadd": best.get("cadd", ""),
            "clinvar": best.get("clinvar", "."),
        })

    scored.sort(key=lambda d: (-d["pheno_score"],
                               -(1 if d["best_impact"] == "HIGH" else 0),
                               -(float(d["best_cadd"]) if d["best_cadd"] else 0.0)))

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    if scored:
        cols = list(scored[0].keys())
        with open(OUT_TSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for d in scored:
                w.writerow(d)

    lines = ["=" * 96,
             " 05 - PHENOTYPE RANKING (blind search)",
             "=" * 96,
             "",
             "## Patient HPO terms"]
    for t in patient:
        lines.append("  {}  {}".format(t, names.get(t, "?")))
    lines += ["",
              "## Candidate genes evaluated: {:,}".format(len(scored)),
              "",
              "## TOP 25",
              "",
              "{:<4}{:<12}{:<9}{:<6}{:<18}{:<5}{:<10}{:<26}{:<20}{:<7}{}".format(
                  "#", "GENE", "SCORE", "HPO", "MODEL", "VAR", "IMPACT",
                  "BEST VARIANT", "HGVSp", "CADD", "CLINVAR")]
    for i, d in enumerate(scored[:25], 1):
        lines.append("{:<4}{:<12}{:<9}{:<6}{:<18}{:<5}{:<10}{:<26}{:<20}{:<7}{}".format(
            i,
            str(d["gene"])[:11],
            d["pheno_score"],
            d["hpo_terms_matched"],
            str(d["model"])[:17],
            d["n_variants"],
            str(d["best_impact"])[:9],
            str(d["best_variant"])[:25],
            str(d["best_hgvsp"])[-19:],
            str(d["best_cadd"])[:6],
            str(d["clinvar"])[:20]))
    lines += ["", "-> " + OUT_TSV]

    txt = "\n".join(lines)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
