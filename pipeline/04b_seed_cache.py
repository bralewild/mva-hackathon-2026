#!/usr/bin/env python3
"""
==============================================================================
04b_seed_cache.py - Rebuild the VEP cache from a previous run

INPUT  : $WORK/04_annotated_candidates.tsv
OUTPUT : $WORK/04_vep_cache.jsonl

PURPOSE
-------
Recovery utility. The first version of 04_frequency_clinical.py did not persist
results per batch: it wrote everything at the end. In that run, 3 of 47 batches
failed with HTTP 500 from Ensembl and their 600 variants were left unannotated,
while the other 8,784 were annotated correctly.

Re-running the whole stage to recover 600 variants would throw away 40 minutes
of valid queries. This script seeds the incremental cache with what was already
annotated, so re-running stage 04 asks Ensembl only for what is genuinely
missing.

A variant counts as "no VEP response" when every annotated field is empty: no
consequence, no rsID, no MANE transcript and no CADD.

Idempotent - safe to run as many times as needed.
==============================================================================
"""
import csv
import json
import os
import sys

BASE = os.path.expanduser("~/mva")
IN_TSV = BASE + "/work/04_annotated_candidates.tsv"
CACHE = BASE + "/work/04_vep_cache.jsonl"


def no_response(r):
    """True when VEP returned nothing for this variant."""
    empty = (".", "", None)
    return (r.get("vep_consequence") in empty
            and r.get("rsid") in empty
            and r.get("mane") in empty
            and not r.get("cadd"))


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("missing " + IN_TSV + " - no previous run to seed from")

    already = set()
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as f:
            for line in f:
                try:
                    already.add(json.loads(line)["key"])
                except Exception:
                    continue
        print("existing cache: {:,} entries".format(len(already)))

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    seeded = missing = duplicate = 0

    with open(CACHE, "a", encoding="utf-8") as out:
        for r in rows:
            key = "{}_{}_{}_{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
            if no_response(r):
                missing += 1
                continue
            if key in already:
                duplicate += 1
                continue
            af = r.get("gnomad_af")
            cadd = r.get("cadd")
            ann = {
                "gnomad_af": float(af) if af not in ("", None) else None,
                "rsid": r.get("rsid", "."),
                "clinvar": r.get("clinvar", "."),
                "mane": r.get("mane", "."),
                "hgvsc": r.get("vep_hgvsc", "."),
                "hgvsp": r.get("vep_hgvsp", "."),
                "cadd": float(cadd) if cadd not in ("", None) else None,
                "consequence": r.get("vep_consequence", "."),
                "impact": r.get("vep_impact", "."),
                "gene": r.get("gene", "."),
            }
            out.write(json.dumps({"key": key, "ann": ann}) + "\n")
            already.add(key)
            seeded += 1

    print("=" * 60)
    print(" 04b - SEED THE VEP CACHE")
    print("=" * 60)
    print("  rows read             {:>8,}".format(len(rows)))
    print("  seeded into cache     {:>8,}".format(seeded))
    print("  already present       {:>8,}".format(duplicate))
    print("  UNannotated (to ask)  {:>8,}".format(missing))
    print()
    print("  -> " + CACHE)
    print()
    print("  Now re-run 04_frequency_clinical.py: it will query Ensembl")
    print("  for the {:,} missing variants only.".format(missing))


if __name__ == "__main__":
    main()
