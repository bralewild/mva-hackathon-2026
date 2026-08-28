#!/usr/bin/env python3
"""
==============================================================================
08_mosaic_aneuploidy.py - Screen for mosaic aneuploidy from the VCF alone

INPUT  : the raw VCF (heterozygous SNVs, streamed with bcftools)
OUTPUT : $RESULTS/08_mosaic_aneuploidy.txt
         $WORK/08_baf_windows.tsv      per-window measurements
         $WORK/08_baf_by_chrom.tsv     per-chromosome summary

PURPOSE
-------
The disease is Mosaic Variegated ANEUPLOIDY. Stages 01-08 identify its cause
(biallelic BUB1B); this stage asks whether its effect - chromosomal instability
itself - is measurable in the same patient, from the same file, with no extra
sequencing. It closes the chain genotype -> mechanism -> cellular phenotype.

THE PRINCIPLE
-------------
At a heterozygous SNV in a diploid cell, reads carrying each allele arrive 1:1,
so the B-allele fraction (BAF) centres on 0.50. If a fraction f of cells carries
an extra copy of that chromosome, one allele is over-represented:

    BAF = (1 + f) / (2 + f)   and its mirror   1 / (2 + f)
    deviation d = |BAF - 0.5| = f / (2 (2 + f))     =>     f = 4d / (1 - 2d)

Depth gives an orthogonal signal: a trisomic fraction f raises mean coverage by
(2 + f) / 2.

WHY WINDOWS, AND WHY THAT MATTERS
---------------------------------
An earlier version summarised each chromosome with the median |BAF-0.5| over its
SNVs and estimated the noise from the spread across 22 chromosomes. That was
wrong twice. At ~44x depth, BAF is a ratio of small integers, so the median lands
on a discrete value: fourteen chromosomes tied at exactly 0.0556, the MAD
collapsed to zero, and the z-scores ran to 10^10 - obvious nonsense.

Worse, individual SNVs are not independent units of evidence for a
whole-chromosome event: a single mismapped region contributes tens of thousands
of correlated sites.

This version bins the genome into fixed windows and treats the WINDOW as the
unit. The mean deviation is continuous, and the genome-wide spread of window
values gives an honest noise scale.

THE DISCRIMINATOR
-----------------
Elevated BAF deviation on a chromosome has two very different explanations:

    real mosaic aneuploidy  ->  the shift is UNIFORM along the chromosome:
                                every window is affected, because every cell in
                                the aneuploid fraction carries the whole extra
                                chromosome

    technical artefact      ->  the shift CONCENTRATES in particular windows:
                                centromeric and pericentromeric repeats,
                                segmental duplications, extreme-GC blocks

So the fraction of a chromosome's windows that are elevated separates the two.
A chromosome-wide event should light up nearly all of its windows.

A GC confounder check is also reported: chromosome GC content is correlated with
the observed excess. GC-rich chromosomes (19, 22, 17, 16, 20) are exactly where
library-preparation bias inflates both coverage and allelic imbalance, so a
strong correlation is evidence against a biological reading.

LIMITS - READ BEFORE INTERPRETING
---------------------------------
* Blood only. Aneuploid fractions vary widely by tissue in MVA; a negative here
  does not exclude mosaicism elsewhere.
* BAF cannot separate trisomy from copy-neutral loss of heterozygosity. Depth
  breaks the tie only for genuine gains.
* No matched normal samples were available, so no external control for
  chromosome-specific technical bias exists. The uniformity test and the GC
  correlation are internal substitutes, not replacements.
* A screening statistic, not a karyotype. Confirmation requires cytogenetics or
  a dedicated mosaicism caller.
==============================================================================
"""
import collections
import math
import os
import statistics
import subprocess
import sys

from _paths import RAW, WORK, RESULTS

VCF = os.path.join(RAW, "WGS_EX2312012_HGWCNDSX7.vcf.gz")
OUT_WIN = os.path.join(WORK, "08_baf_windows.tsv")
OUT_CHR = os.path.join(WORK, "08_baf_by_chrom.tsv")
OUT_TXT = os.path.join(RESULTS, "08_mosaic_aneuploidy.txt")

MIN_DP = 20
MIN_GQ = 20
WINDOW = 10_000_000          # 10 Mb
MIN_SITES_WINDOW = 2000      # a window needs enough sites to be meaningful
MIN_WINDOWS_CHROM = 4
AUTOSOMES = [str(i) for i in range(1, 23)]

# GRCh38 GC content per chromosome (%), for the confounder check.
GC = {"1": 41.7, "2": 40.2, "3": 39.7, "4": 38.2, "5": 39.5, "6": 39.6,
      "7": 40.7, "8": 40.2, "9": 42.3, "10": 41.6, "11": 41.6, "12": 40.8,
      "13": 38.5, "14": 40.9, "15": 42.2, "16": 44.8, "17": 45.5, "18": 39.8,
      "19": 48.4, "20": 44.1, "21": 40.9, "22": 47.7}


def stream_sites():
    """Heterozygous biallelic PASS SNVs, streamed with bcftools."""
    cmd = ["bcftools", "query",
           "-i", 'TYPE="snp" && N_ALT=1 && FILTER="PASS" && GT="het" '
                 '&& FMT/DP>={} && FMT/GQ>={}'.format(MIN_DP, MIN_GQ),
           "-f", "%CHROM\t%POS[\t%AD\t%DP]\n", VCF]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[0] not in AUTOSOMES or "," not in f[2]:
            continue
        try:
            ref, alt = (int(x) for x in f[2].split(",")[:2])
            dp, pos = int(f[3]), int(f[1])
        except ValueError:
            continue
        tot = ref + alt
        if tot >= MIN_DP:
            yield f[0], pos, alt / tot, dp
    p.wait()
    if p.returncode != 0:
        sys.exit("bcftools failed: " + p.stderr.read()[:300])


def mosaic_fraction(d):
    if d <= 0:
        return 0.0
    denom = 1 - 2 * d
    return (4 * d / denom) if denom > 0 else float("inf")


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else 0.0


def sparkline(values, lo, hi, bins=20):
    counts = [0] * bins
    for v in values:
        if lo <= v <= hi:
            counts[min(int((v - lo) / (hi - lo) * bins), bins - 1)] += 1
    top = max(counts) or 1
    blocks = " .:-=+*#%@"
    return "".join(blocks[min(int(c / top * (len(blocks) - 1)), len(blocks) - 1)] for c in counts)


def main():
    if not os.path.exists(VCF):
        sys.exit("missing " + VCF)

    win = collections.defaultdict(lambda: {"dev": [], "dp": []})
    total = 0
    for chrom, pos, baf, dp in stream_sites():
        w = (chrom, pos // WINDOW)
        win[w]["dev"].append(abs(baf - 0.5))
        win[w]["dp"].append(dp)
        total += 1
    if not total:
        sys.exit("no heterozygous sites passed the thresholds")

    # ---- per window ----
    windows = []
    for (chrom, wi), d in sorted(win.items(), key=lambda kv: (int(kv[0][0]), kv[0][1])):
        if len(d["dev"]) < MIN_SITES_WINDOW:
            continue
        windows.append({
            "chrom": chrom, "start": wi * WINDOW, "n_sites": len(d["dev"]),
            "mean_dev": statistics.fmean(d["dev"]),
            "median_dp": statistics.median(d["dp"]),
        })

    all_dev = [w["mean_dev"] for w in windows]
    all_dp = [w["median_dp"] for w in windows]
    base_dev = statistics.median(all_dev)
    base_dp = statistics.median(all_dp)
    mad = statistics.median(abs(v - base_dev) for v in all_dev)
    sigma = mad * 1.4826
    if sigma <= 0:                      # never trust a zero scale again
        sigma = statistics.pstdev(all_dev) or 1e-6

    with open(OUT_WIN, "w") as f:
        f.write("chrom\tstart\tn_sites\tmean_dev\tmedian_dp\texcess\tz\n")
        for w in windows:
            e = w["mean_dev"] - base_dev
            f.write("{}\t{}\t{}\t{:.6f}\t{}\t{:.6f}\t{:.2f}\n".format(
                w["chrom"], w["start"], w["n_sites"], w["mean_dev"],
                w["median_dp"], e, e / sigma))

    # ---- per chromosome, windows as the unit of evidence ----
    by_chrom = collections.defaultdict(list)
    for w in windows:
        by_chrom[w["chrom"]].append(w)

    rows = []
    for c in AUTOSOMES:
        ws = by_chrom.get(c, [])
        if len(ws) < MIN_WINDOWS_CHROM:
            continue
        devs = [w["mean_dev"] for w in ws]
        m = statistics.fmean(devs)
        sem = statistics.stdev(devs) / math.sqrt(len(devs)) if len(devs) > 1 else 0.0
        excess = m - base_dev
        elevated = sum(1 for v in devs if v - base_dev > sigma)
        rows.append({
            "chrom": c, "n_windows": len(ws), "n_sites": sum(w["n_sites"] for w in ws),
            "mean_dev": m, "sem": sem, "excess": excess,
            "z": excess / sigma,
            "frac_windows_elevated": elevated / len(ws),
            "depth_ratio": statistics.fmean([w["median_dp"] for w in ws]) / base_dp,
            "est_mosaic_fraction": mosaic_fraction(excess),
            "gc": GC.get(c, 0.0),
            "_devs": devs,
        })
    rows.sort(key=lambda r: -r["excess"])

    with open(OUT_CHR, "w") as f:
        cols = [k for k in rows[0] if not k.startswith("_")]
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("{:.6g}".format(r[c]) if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")

    # ---- GC confounder ----
    r_gc = pearson([r["gc"] for r in rows], [r["excess"] for r in rows])
    r_gc_dp = pearson([r["gc"] for r in rows], [r["depth_ratio"] for r in rows])

    Z_FLAG, UNIFORM_MIN = 3.0, 0.80
    flagged = [r for r in rows if r["z"] >= Z_FLAG]
    uniform = [r for r in flagged if r["frac_windows_elevated"] >= UNIFORM_MIN]
    patchy = [r for r in flagged if r["frac_windows_elevated"] < UNIFORM_MIN]

    L = ["=" * 100,
         " 09 - MOSAIC ANEUPLOIDY SCREEN (B-allele fraction, from the VCF alone)",
         "=" * 100, "",
         "Heterozygous SNVs : {:,}   (PASS, biallelic, DP>={}, GQ>={})".format(total, MIN_DP, MIN_GQ),
         "Windows           : {:,} of {:,} bp, >= {:,} sites each".format(len(windows), WINDOW, MIN_SITES_WINDOW),
         "Autosomes         : {}".format(len(rows)),
         "",
         "Noise scale is taken from the spread of WINDOW values, not from SNV counts:",
         "a single mismapped region contributes tens of thousands of correlated sites,",
         "so SNVs are not independent evidence for a whole-chromosome event.",
         "",
         "  baseline mean |BAF-0.5| : {:.5f}".format(base_dev),
         "  robust sigma (window)   : {:.5f}".format(sigma),
         "  baseline depth          : {:.0f}x".format(base_dp),
         "",
         "{:<6}{:>8}{:>11}{:>11}{:>8}{:>9}{:>9}{:>8}{:>8}  {}".format(
             "chr", "wins", "mean dev", "excess", "z", "uniform", "depth", "mosaic", "GC%",
             "window deviations"),
         "-" * 100]
    lo = base_dev - 3 * sigma
    hi = base_dev + 6 * sigma
    for r in rows:
        flag = "  <<<" if r["z"] >= Z_FLAG else ""
        L.append("{:<6}{:>8}{:>11.5f}{:>11.5f}{:>8.1f}{:>8.0%}{:>8.2f}x{:>7.1%}{:>8.1f}  {}{}".format(
            r["chrom"], r["n_windows"], r["mean_dev"], r["excess"], r["z"],
            r["frac_windows_elevated"], r["depth_ratio"],
            min(r["est_mosaic_fraction"], 9.99), r["gc"],
            sparkline(r["_devs"], lo, hi), flag))

    L += ["", "=" * 100, " CONFOUNDER CHECK: GC CONTENT", "=" * 100, "",
          "  correlation( chromosome GC% , BAF excess )   r = {:+.3f}".format(r_gc),
          "  correlation( chromosome GC% , depth ratio )  r = {:+.3f}".format(r_gc_dp),
          ""]
    if abs(r_gc) >= 0.6 or abs(r_gc_dp) >= 0.6:
        L.append("  STRONG correlation with GC content. Library-preparation GC bias")
        L.append("  inflates both coverage and allelic imbalance on GC-rich chromosomes")
        L.append("  and is a sufficient explanation for the pattern on its own.")
        L.append("  The signal below CANNOT be attributed to biology from this data.")
    elif abs(r_gc) >= 0.3:
        L.append("  MODERATE correlation with GC content. GC bias contributes and cannot")
        L.append("  be excluded as at least a partial explanation.")
    else:
        L.append("  WEAK correlation with GC content. GC bias alone does not explain")
        L.append("  the observed pattern.")

    L += ["", "=" * 100, " INTERPRETATION", "=" * 100, ""]
    if not flagged:
        L.append("  No autosome exceeds z = {} on the window scale.".format(Z_FLAG))
        L.append("  No whole-chromosome mosaicism is detected above the floor below.")
    else:
        if uniform:
            L.append("  UNIFORM elevation (>= {:.0%} of windows) - the pattern a real".format(UNIFORM_MIN))
            L.append("  whole-chromosome event produces:")
            L.append("")
            for r in uniform:
                supp = "supported" if r["depth_ratio"] > 1.02 else "NOT supported"
                L.append("    chr{:<4} excess {:.5f}  ~{:.0%} of cells   depth {:.2f}x ({} by coverage)".format(
                    r["chrom"], r["excess"], min(r["est_mosaic_fraction"], 1.0),
                    r["depth_ratio"], supp))
            L.append("")
        if patchy:
            L.append("  PATCHY elevation (< {:.0%} of windows) - concentrated in part of the".format(UNIFORM_MIN))
            L.append("  chromosome, which is what mapping artefacts look like, not aneuploidy:")
            L.append("")
            for r in patchy:
                L.append("    chr{:<4} excess {:.5f}   only {:.0%} of windows elevated".format(
                    r["chrom"], r["excess"], r["frac_windows_elevated"]))
            L.append("")

    floor = mosaic_fraction(Z_FLAG * sigma)
    L += ["## Detection floor", "",
          "  A window-scale z of 3 corresponds to an excess of {:.5f},".format(Z_FLAG * sigma),
          "  i.e. a mosaic fraction of roughly {:.0%}. Below that, nothing is called.".format(min(floor, 1.0)),
          "",
          "## Limits", "",
          "  - Blood only. Aneuploid fractions vary widely by tissue in MVA.",
          "  - BAF cannot separate trisomy from copy-neutral LOH; depth breaks the tie",
          "    only for genuine gains.",
          "  - No matched normal samples were available, so there is no external control",
          "    for chromosome-specific technical bias. The uniformity test and the GC",
          "    correlation are internal substitutes, not replacements.",
          "  - Whole-chromosome screen; segmental events are not resolved.",
          "  - A screening statistic, not a karyotype.",
          "",
          "  -> {}".format(OUT_CHR),
          "  -> {}".format(OUT_WIN)]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
