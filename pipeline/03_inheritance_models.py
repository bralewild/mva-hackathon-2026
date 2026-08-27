#!/usr/bin/env python3
"""
==============================================================================
03_inheritance_models.py — Filtrado por calidad, impacto y modelo de herencia

ENTRADA : $WORK/02_annotated.vcf.gz   (salida de snpEff, genome-wide)
SALIDA  : $WORK/03_candidates.tsv     (una fila por variante candidata)
          $RESULTS/03_inheritance_summary.txt

PROPOSITO
---------
Reducir 5.012.204 variantes a un conjunto manejable SIN usar ningun
conocimiento previo sobre la enfermedad. El pipeline no sabe que se trata de
MVA ni que existen BUB1B / CEP57 / TRIP13.

Criterios (justificados en docs/02_metodologia.md):
  1. FILTER == PASS
  2. GQ >= MIN_GQ  y  DP >= MIN_DP        -> confianza en el genotipo
  3. Impacto snpEff HIGH o MODERATE       -> potencial consecuencia funcional
  4. Modelo de herencia:
       AR_COMPOUND_HET : gen con >= 2 variantes heterocigotas
       AR_HOMOZYGOUS   : gen con >= 1 variante homocigota alternativa
       (no hay trio, asi que 'de novo' no es evaluable)

El paciente es singleton: no se puede establecer fase por pedigri. Las
parejas se marcan como "presuntas en trans" y se anota si GATK dejo phasing
fisico (PID/PGT) que permita confirmarlas o descartarlas.
==============================================================================
"""
import gzip, os, sys, collections

BASE     = os.path.expanduser("~/mva")
IN_VCF   = f"{BASE}/work/02_annotated.vcf.gz"
OUT_TSV  = f"{BASE}/work/03_candidates.tsv"
OUT_SUM  = f"{BASE}/results/03_inheritance_summary.txt"
MIN_GQ, MIN_DP = 20, 10
KEEP_IMPACT = {"HIGH", "MODERATE"}

# ANN de snpEff: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|
#                Feature_ID|BioType|Rank|HGVS.c|HGVS.p|...
A_EFFECT, A_IMPACT, A_GENE, A_GENEID, A_FEATID, A_RANK, A_HGVSC, A_HGVSP = 1,2,3,4,6,8,9,10

def main():
    if not os.path.exists(IN_VCF):
        sys.exit(f"falta {IN_VCF} — corre antes 02_annotate_genomewide.sh")

    stats = collections.Counter()
    by_gene = collections.defaultdict(list)

    with gzip.open(IN_VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                continue
            stats["total"] += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, flt, info, fmt, smp = f[0], int(f[1]), f[3], f[4], f[6], f[7], f[8], f[9]

            if flt not in ("PASS", "."):
                stats["descartada_filtro"] += 1; continue
            stats["pass"] += 1

            fd = dict(zip(fmt.split(":"), smp.split(":")))
            gt = fd.get("GT", "./.").replace("|", "/")
            if gt in ("./.", "0/0"):
                stats["descartada_ref_o_nocall"] += 1; continue
            try:
                gq, dp = int(fd.get("GQ", 0)), int(fd.get("DP", 0))
            except ValueError:
                gq, dp = 0, 0
            if gq < MIN_GQ or dp < MIN_DP:
                stats["descartada_calidad_gt"] += 1; continue
            stats["calidad_ok"] += 1

            ann = None
            for kv in info.split(";"):
                if kv.startswith("ANN="):
                    ann = kv[4:]; break
            if not ann:
                stats["sin_anotacion"] += 1; continue

            best = None
            for a in ann.split(","):
                p = a.split("|")
                if len(p) > A_HGVSP and p[A_IMPACT] in KEEP_IMPACT:
                    best = p; break
            if best is None:
                stats["impacto_bajo"] += 1; continue
            stats["impacto_alto_o_moderado"] += 1

            zyg = "hom" if gt in ("1/1", "2/2") else "het"
            stats[f"zig_{zyg}"] += 1
            gene = best[A_GENE] or best[A_GENEID]
            by_gene[gene].append(dict(
                chrom=chrom, pos=pos, ref=ref, alt=alt, gt=gt, zyg=zyg,
                gq=gq, dp=dp, ad=fd.get("AD", "."), pid=fd.get("PID", "."), pgt=fd.get("PGT", "."),
                effect=best[A_EFFECT], impact=best[A_IMPACT], gene=gene,
                gene_id=best[A_GENEID], feature=best[A_FEATID], rank=best[A_RANK],
                hgvsc=best[A_HGVSC], hgvsp=best[A_HGVSP],
            ))

    # --- modelos de herencia ---
    rows, models = [], collections.Counter()
    for gene, vs in by_gene.items():
        hets = [v for v in vs if v["zyg"] == "het"]
        homs = [v for v in vs if v["zyg"] == "hom"]
        model = None
        if homs:                model = "AR_HOMOZYGOUS"
        elif len(hets) >= 2:    model = "AR_COMPOUND_HET"
        if not model:
            continue
        models[model] += 1
        for v in (homs if homs else hets):
            v["model"] = model
            v["n_het_gen"] = len(hets)
            v["n_hom_gen"] = len(homs)
            rows.append(v)

    cols = ["gene","model","chrom","pos","ref","alt","gt","zyg","impact","effect",
            "hgvsc","hgvsp","rank","gq","dp","ad","pgt","pid","feature","gene_id",
            "n_het_gen","n_hom_gen"]
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["gene"], x["chrom"], x["pos"])):
            out.write("\t".join(str(r.get(c, ".")) for c in cols) + "\n")

    lines = ["="*66, " 03 — MODELOS DE HERENCIA (busqueda ciega)", "="*66, "",
             "## Embudo"]
    for k in ["total","descartada_filtro","pass","descartada_ref_o_nocall",
              "descartada_calidad_gt","calidad_ok","sin_anotacion","impacto_bajo",
              "impacto_alto_o_moderado","zig_het","zig_hom"]:
        lines.append(f"  {k:<28} {stats[k]:>10,}")
    lines += ["", "## Genes por modelo"]
    for m, n in models.most_common():
        lines.append(f"  {m:<20} {n:>8,} genes")
    lines += ["", f"## Variantes candidatas escritas: {len(rows):,}", f"   -> {OUT_TSV}"]
    txt = "\n".join(lines)
    os.makedirs(os.path.dirname(OUT_SUM), exist_ok=True)
    open(OUT_SUM, "w").write(txt + "\n")
    print(txt)

if __name__ == "__main__":
    main()
