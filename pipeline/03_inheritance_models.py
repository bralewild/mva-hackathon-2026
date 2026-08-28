#!/usr/bin/env python3
"""
==============================================================================
03_inheritance_models.py - Quality, coherence, impact and inheritance model

INPUT  : $WORK/02_annotated.vcf.gz   (snpEff output, genome-wide)
OUTPUT : $WORK/03_candidates.tsv     (one row per candidate variant)
         $RESULTS/03_inheritance_summary.txt

PURPOSE
-------
Reduce 5,012,204 variants to a tractable set WITHOUT using any prior knowledge
about the disease. The pipeline does not know this is MVA, nor that BUB1B /
CEP57 / TRIP13 exist.

Criteria (rationale in docs/01_pipeline_flow.md):
  1. FILTER == PASS
  2. GQ >= MIN_GQ and DP >= MIN_DP           -> genotype confidence
  3. VAF coherent with the called zygosity   -> biological plausibility
  4. snpEff impact HIGH or MODERATE          -> potential functional consequence
  5. Inheritance model:
       AR_COMPOUND_HET : gene with >= 2 heterozygous variants
       AR_HOMOZYGOUS   : gene with >= 1 homozygous-alternate variant
       (no trio, so 'de novo' is not evaluable)

The patient is a singleton: phase cannot be established from pedigree. Pairs are
flagged as "presumed in trans", and the script records whether GATK left
physical phasing (PID/PGT) that could confirm or exclude them.
==============================================================================
"""
import collections
import gzip
import os
import sys

BASE = os.path.expanduser("~/mva")
IN_VCF = BASE + "/work/02_annotated.vcf.gz"
OUT_TSV = BASE + "/work/03_candidates.tsv"
OUT_SUM = BASE + "/results/03_inheritance_summary.txt"

MIN_GQ, MIN_DP = 20, 10

# Acceptable alternate allele fraction (VAF) for a real heterozygous call.
# Without this, mismapping artefacts get through: reads from a paralogue or
# pseudogene align to the real gene and the caller reports them as heterozygous
# with GQ=99 and high DP but VAF ~0.15. Witness case: SERPINA1 (14q32.13, next
# to the SERPINA2 pseudogene) contributed 201 variants at VAF 0.13-0.20,
# including 4 spurious frameshifts.
VAF_HET_MIN, VAF_HET_MAX = 0.25, 0.75
VAF_HOM_MIN = 0.85

KEEP_IMPACT = {"HIGH", "MODERATE"}

# snpEff ANN field: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|
#                   Feature_ID|BioType|Rank|HGVS.c|HGVS.p|...
A_EFFECT, A_IMPACT, A_GENE, A_GENEID, A_FEATID, A_RANK, A_HGVSC, A_HGVSP = 1, 2, 3, 4, 6, 8, 9, 10


def main():
    if not os.path.exists(IN_VCF):
        sys.exit("missing " + IN_VCF + " - run 02_annotate_genomewide.sh first")

    stats = collections.Counter()
    by_gene = collections.defaultdict(list)

    with gzip.open(IN_VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            stats["total"] += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            flt, info, fmt, smp = f[6], f[7], f[8], f[9]

            if flt not in ("PASS", "."):
                stats["dropped_filter"] += 1
                continue
            stats["pass"] += 1

            fd = dict(zip(fmt.split(":"), smp.split(":")))
            gt = fd.get("GT", "./.").replace("|", "/")
            if gt in ("./.", "0/0"):
                stats["dropped_ref_or_nocall"] += 1
                continue
            try:
                gq, dp = int(fd.get("GQ", 0)), int(fd.get("DP", 0))
            except ValueError:
                gq, dp = 0, 0
            if gq < MIN_GQ or dp < MIN_DP:
                stats["dropped_genotype_quality"] += 1
                continue

            # Biological coherence of the genotype. A high GQ means the CALLER is
            # confident, NOT that the variant is real: mismapped paralogue reads
            # yield het calls with GQ=99, DP>50 and VAF ~0.15.
            ad = fd.get("AD", "")
            vaf = None
            if "," in ad:
                try:
                    parts = [int(x) for x in ad.split(",")]
                    tot = sum(parts)
                    if tot > 0:
                        vaf = max(parts[1:]) / tot
                except ValueError:
                    vaf = None
            is_hom = gt in ("1/1", "2/2")
            if vaf is not None:
                if is_hom and vaf < VAF_HOM_MIN:
                    stats["dropped_vaf_incoherent"] += 1
                    continue
                if not is_hom and not (VAF_HET_MIN <= vaf <= VAF_HET_MAX):
                    stats["dropped_vaf_incoherent"] += 1
                    continue
            stats["quality_ok"] += 1

            ann = None
            for kv in info.split(";"):
                if kv.startswith("ANN="):
                    ann = kv[4:]
                    break
            if not ann:
                stats["no_annotation"] += 1
                continue

            best = None
            for a in ann.split(","):
                p = a.split("|")
                if len(p) > A_HGVSP and p[A_IMPACT] in KEEP_IMPACT:
                    best = p
                    break
            if best is None:
                stats["low_impact"] += 1
                continue
            stats["impact_high_or_moderate"] += 1

            zyg = "hom" if is_hom else "het"
            stats["zyg_" + zyg] += 1
            gene = best[A_GENE] or best[A_GENEID]
            by_gene[gene].append(dict(
                chrom=chrom, pos=pos, ref=ref, alt=alt, gt=gt, zyg=zyg,
                gq=gq, dp=dp, ad=fd.get("AD", "."),
                pid=fd.get("PID", "."), pgt=fd.get("PGT", "."),
                effect=best[A_EFFECT], impact=best[A_IMPACT], gene=gene,
                gene_id=best[A_GENEID], feature=best[A_FEATID], rank=best[A_RANK],
                hgvsc=best[A_HGVSC], hgvsp=best[A_HGVSP],
            ))

    # --- inheritance models ---
    rows, models = [], collections.Counter()
    for gene, vs in by_gene.items():
        hets = [v for v in vs if v["zyg"] == "het"]
        homs = [v for v in vs if v["zyg"] == "hom"]
        model = None
        if homs:
            model = "AR_HOMOZYGOUS"
        elif len(hets) >= 2:
            model = "AR_COMPOUND_HET"
        if not model:
            continue
        models[model] += 1
        for v in (homs if homs else hets):
            v["model"] = model
            v["n_het_gene"] = len(hets)
            v["n_hom_gene"] = len(homs)
            rows.append(v)

    cols = ["gene", "model", "chrom", "pos", "ref", "alt", "gt", "zyg", "impact",
            "effect", "hgvsc", "hgvsp", "rank", "gq", "dp", "ad", "pgt", "pid",
            "feature", "gene_id", "n_het_gene", "n_hom_gene"]
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["gene"], x["chrom"], x["pos"])):
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")

    lines = ["=" * 66, " 03 - INHERITANCE MODELS (blind search)", "=" * 66, "",
             "## Reduction funnel"]
    for k in ["total", "dropped_filter", "pass", "dropped_ref_or_nocall",
              "dropped_genotype_quality", "dropped_vaf_incoherent", "quality_ok",
              "no_annotation", "low_impact", "impact_high_or_moderate",
              "zyg_het", "zyg_hom"]:
        lines.append("  {:<28} {:>10,}".format(k, stats[k]))
    lines += ["", "## Genes per model"]
    for m, n in models.most_common():
        lines.append("  {:<20} {:>8,} genes".format(m, n))
    lines += ["", "## Candidate variants written: {:,}".format(len(rows)),
              "   -> " + OUT_TSV]
    txt = "\n".join(lines)
    os.makedirs(os.path.dirname(OUT_SUM), exist_ok=True)
    open(OUT_SUM, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
