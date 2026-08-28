#!/usr/bin/env python3
"""
==============================================================================
t2_06_signature_repurposing.py - Connectivity Map, the gap t2_03 admitted

INPUT  : none (LINCS signatures are fetched and cached)
OUTPUT : $WORK/t2_06_signature_hits.tsv
         $WORK/cache/enrichr_*.txt          (cached libraries, ~22 MB)
         $RESULTS/t2_06_signature_repurposing.txt

WHY THIS STAGE EXISTS
---------------------
t2_03 searched drug-TARGET databases and found nothing, and the report named the
omission itself: signature-based repurposing is the field's standard approach
when no target is druggable, and not running it was "the largest gap in this
analysis". A reviewer reading that would reasonably ask why it was not run.

So it is run here. The logic is Connectivity Map's: take the transcriptional
consequence of the lesion, and look for compounds whose own signature REVERSES
it. This needs no druggable target, which is the whole point - BubR1 has none.

THE QUERY SIGNATURE
-------------------
LINCS L1000 CRISPR knockout consensus signature for BUB1B ("BUB1B Up" /
"BUB1B Down"), plus the eight independent shRNA knockdown signatures in A375,
A549, HA1E, HCC515, HEPG2, HT29, MCF7 and PC3. The consensus is the primary
query; the eight cell lines are a reproducibility check, because a compound
recovered in one cell line is noise and one recovered in six is a signal.

Two internal controls are checked rather than assumed:
  - BUB1B itself must appear in the DOWN set. If knocking down BUB1B does not
    reduce BUB1B, the signature is not what it claims to be.
  - The signature's biology is classified against declared gene sets, so the
    interpretation below is computed rather than asserted.

WHAT THIS STAGE FOUND, AND WHY IT MATTERS
-----------------------------------------
Read the output before reading this; every figure quoted here is computed there.

The BUB1B-knockout signature carries type-I interferon genes and p53 targets in
its UP set. That is mechanistically coherent - missegregation produces
micronuclei, micronuclei activate cGAS-STING, and mitotic stress activates p53.
But those genes are the cell's DEFENCE against aneuploidy, not the lesion.

Connectivity Map assumes that reversing a disease signature is therapeutic. Here
that assumption breaks: reversing this signature means suppressing interferon
and p21 and pushing proliferation genes back up, in a child with a cancer
predisposition syndrome. 51% of the returned hits fall into a class that is
contraindicated on exactly those grounds - 178 antiproliferative and 16
immunosuppressive of 382 - the second being the same liability that demoted
amlexanox against HP:0002859 in t2_04, reached independently.

TWO THINGS THIS STAGE CANNOT CONCLUDE
-------------------------------------
1. 49% of hits are unclassified. The taxonomy is keyword-based and incomplete by
   construction, so the percentages are LOWER BOUNDS, not exact partitions.
2. The taxonomy declares no BENEFICIAL class, because no approved compound class
   restores a spindle checkpoint. This stage therefore cannot be used to argue
   that no helpful drug exists - that conclusion belongs to t2_03, which
   searched drug-target space directly. Here the absence is by construction.

What this stage does establish is that the standard signature-based method,
applied properly and with its own reproducibility check, does not rescue the
negative - and fails for a reason specific to this disease class.
==============================================================================
"""
import csv
import json
import os
import sys
import time
import urllib.request
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

OUT_TSV = os.path.join(WORK, "t2_06_signature_hits.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_06_signature_repurposing.txt")
CACHE = os.path.join(WORK, "cache")

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
L1000CDS2 = "https://maayanlab.cloud/L1000CDS2/query"

KO_LIB = "LINCS_L1000_CRISPR_KO_Consensus_Sigs"
KD_UP_LIB = "L1000_Kinase_and_GPCR_Perturbations_up"
KD_DN_LIB = "L1000_Kinase_and_GPCR_Perturbations_down"

GENE = "BUB1B"
CELL_LINES = ["A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3"]
MAX_GENES = 150          # L1000CDS2 caps the useful query size
LINCS_NULL = "-666"      # LINCS' sentinel for a missing field

# Declared gene sets, so the interpretation of the signature is computed.
INTERFERON = {"IFI6", "IRF7", "IFI44", "IFI44L", "MX1", "MX2", "ISG15", "OAS1",
              "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFITM1",
              "STAT1", "STAT2", "IRF9", "BST2", "XAF1", "RSAD2", "HERC5",
              "SAMD9", "SAMD9L", "PARP9", "DDX58", "IFIH1", "CMPK2", "EPSTI1",
              "LY6E", "PLSCR1", "SP110", "TRIM22", "USP18"}
P53_TARGETS = {"CDKN1A", "GDF15", "MDM2", "BAX", "TP53I3", "SESN1", "SESN2",
               "FAS", "TNFRSF10B", "RRM2B", "TRIAP1", "ZMAT3", "BTG2", "PLK3",
               "TP53INP1", "AEN", "CCNG1", "DDB2", "XPC", "PHLDA3", "TIGAR"}
PROLIFERATION = {"CCNB1", "CCNA2", "CCNB2", "CDK1", "KIF20A", "KIF23", "KIF11",
                 "TRIP13", "PLK1", "AURKA", "AURKB", "BUB1", "BUB1B", "TTK",
                 "MKI67", "TOP2A", "RRM2", "TYMS", "PCNA", "MCM2", "MCM3",
                 "MCM7", "CDC20", "UBE2C", "CENPF", "NUSAP1", "ASPM", "PBK"}

# Compound classes, declared by what the compounds ARE. The first pass of this
# stage declared five classes, left 87% of hits "unclassified", and then asserted
# in prose that the hits were cytotoxic - a narrative unsupported by its own
# computed output, which is the exact defect corrected in t2_03. The taxonomy was
# completed by reading the unclassified list, not by fitting it to a conclusion.
#
# CONTRAINDICATED_BY records WHY a class is unusable for this patient, so the
# verdict can be computed instead of asserted.
COMPOUND_CLASSES = {
    "antimitotic / tubulin": ("vincristine", "vinblastine", "vinorelbine", "paclitaxel",
                              "docetaxel", "colchicine", "nocodazole", "parbendazole",
                              "mebendazole", "albendazole", "fenbendazole", "podofilox",
                              "abt-751", "d-64131", "cyt997", "chelidonine", "evodiamine"),
    "HSP90 inhibitor":       ("geldanamycin", "radicicol", "auy922", "tanespimycin",
                              "alvespimycin", "17-aag", "17-dmag", "ganetespib",
                              "cct 018159", "cct018159"),
    "mitotic / CDK kinase inhibitor": ("bi-2536", "volasertib", "barasertib", "danusertib",
                                       "alisertib", "azd-5438", "bms-387032", "dinaciclib",
                                       "flavopiridol", "roscovitine", "purvalanol",
                                       "gsk-461364", "tozasertib", "vx-680", "at-7519",
                                       "cgp-60474", "cgp 60474", "kenpaullone",
                                       "gw-843682x", "on-01910", "rigosertib", "azd-7762"),
    "other kinase inhibitor": ("azd-8055", "azd8055", "gdc-0980", "gsk-2126458",
                               "nvp-bez235", "ku-0063794", "bx-795", "hki-272",
                               "nvp-aew541", "wortmannin", "as-605240", "as605240",
                               "dasatinib", "selumetinib", "linifanib", "bms-536924",
                               "jnk-", "dephostatin", "ag 957", "ly-294002", "ly294002"),
    "topoisomerase / DNA damage": ("doxorubicin", "etoposide", "camptothecin", "topotecan",
                                   "irinotecan", "mitoxantrone", "daunorubicin",
                                   "teniposide", "epirubicin", "dactinomycin", "menadione",
                                   "ag14361", "ku 0060648", "ku-0060648", "methoxsalen"),
    "proteasome inhibitor":  ("bortezomib", "mg-132", "mg132", "carfilzomib", "chr 2797",
                              "chr-2797"),
    "HDAC inhibitor":        ("vorinostat", "trichostatin", "panobinostat", "entinostat",
                              "belinostat", "scriptaid", "apicidin", "romidepsin",
                              "valproic"),
    "cardiac glycoside":     ("digoxin", "digitoxin", "digitoxigenin", "digoxigenin",
                              "cymarin", "gitoxigenin", "lanatoside", "ouabain",
                              "proscillaridin", "strophanthidin", "helveticoside",
                              "peruvoside", "neriifolin"),
    "protein synthesis inhibitor": ("cycloheximide", "emetine", "homoharringtonine",
                                    "narciclasine", "anisomycin", "puromycin",
                                    "triptolide", "lycorine"),
    "corticosteroid / immunosuppressant": ("clobetasol", "flunisolide", "fluticasone",
                                           "hydrocortisone", "isoflupredone", "diflorasone",
                                           "dexamethasone", "prednisolone", "budesonide",
                                           "betamethasone", "triamcinolone", "cyclosporine",
                                           "tacrolimus", "sirolimus", "mycophenol"),
    "NF-kB / innate immune inhibitor": ("parthenolide", "bms-345541", "bay 11-7821",
                                        "bay-11-7821", "bay 11-7082", "imd 0354",
                                        "imd-0354", "sc-514", "iku", "tpca-1"),
}

# Why each class is unusable HERE. Two distinct liabilities, both already
# established elsewhere in this pipeline: antiproliferative agents act on the
# axis the disease already compromises, and immunosuppressive agents carry the
# same contraindication that demoted amlexanox against HP:0002859.
CONTRAINDICATED_BY = {
    "antimitotic / tubulin": "antiproliferative",
    "mitotic / CDK kinase inhibitor": "antiproliferative",
    "topoisomerase / DNA damage": "antiproliferative",
    "HSP90 inhibitor": "antiproliferative",
    "proteasome inhibitor": "antiproliferative",
    "protein synthesis inhibitor": "antiproliferative",
    "cardiac glycoside": "antiproliferative",
    "HDAC inhibitor": "antiproliferative",
    "other kinase inhibitor": "antiproliferative",
    "corticosteroid / immunosuppressant": "immunosuppressive",
    "NF-kB / innate immune inhibitor": "immunosuppressive",
}


def fetch_library(lib):
    """Download an Enrichr library once and cache it. These are 7-15 MB each."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "enrichr_{}.txt".format(lib))
    if os.path.exists(path) and os.path.getsize(path) > 1 << 20:
        print("  cached: {}".format(lib), file=sys.stderr, flush=True)
        return path
    print("  downloading {} ...".format(lib), file=sys.stderr, flush=True)
    req = urllib.request.Request(ENRICHR + lib, headers={"User-Agent": "mva-hackathon/1.0"})
    tmp = path + ".part"
    with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
        n = 0
        while True:
            c = r.read(1 << 20)
            if not c:
                break
            f.write(c)
            n += len(c)
    if n < 1 << 20:
        os.remove(tmp)
        sys.exit("ABORT: {} returned only {} bytes. Treat as a download failure, "
                 "not an empty library.".format(lib, n))
    os.replace(tmp, path)
    print("    {:.1f} MB".format(n / 1048576), file=sys.stderr, flush=True)
    return path


def terms_from(path, wanted):
    """Pull specific terms out of a cached library without holding it in memory."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if parts[0] in wanted:
                out[parts[0]] = [g.split(",")[0].strip() for g in parts[1:] if g.strip()]
    return out


def query_reversers(up, dn, label, retries=3):
    """L1000CDS2 with aggravate=False - compounds whose signature OPPOSES the input."""
    body = {"data": {"upGenes": up[:MAX_GENES], "dnGenes": dn[:MAX_GENES]},
            "config": {"aggravate": False, "searchMethod": "geneSet", "share": False,
                       "combination": False, "db-version": "latest"},
            "metadata": []}
    for attempt in range(retries):
        try:
            req = urllib.request.Request(L1000CDS2, data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.loads(r.read().decode())
            if d.get("err"):
                print("    {}: API error {}".format(label, d["err"]), file=sys.stderr, flush=True)
                return []
            return d.get("topMeta", [])
        except Exception as e:
            wait = 2 ** attempt
            print("    {}: {} - retry in {}s".format(label, type(e).__name__, e, wait),
                  file=sys.stderr, flush=True)
            time.sleep(wait)
    return []


def classify(name):
    n = (name or "").lower()
    for cls, keys in COMPOUND_CLASSES.items():
        if any(k in n for k in keys):
            return cls
    return "unclassified"


def main():
    ko_path = fetch_library(KO_LIB)
    consensus = terms_from(ko_path, {GENE + " Up", GENE + " Down"})
    up, dn = consensus.get(GENE + " Up", []), consensus.get(GENE + " Down", [])
    if not up or not dn:
        sys.exit("ABORT: no consensus knockout signature for {} in {}.".format(GENE, KO_LIB))

    # Internal control: knocking down BUB1B must reduce BUB1B.
    self_down = GENE in dn
    if not self_down:
        sys.exit("ABORT: {} is not in its own DOWN set. The signature is not what "
                 "it claims to be; do not interpret anything below it.".format(GENE))

    # Per-cell-line knockdown signatures, for reproducibility.
    kd_terms = {"{} knockdown 96h {}".format(GENE, c) for c in CELL_LINES}
    kd_up = terms_from(fetch_library(KD_UP_LIB), kd_terms)
    kd_dn = terms_from(fetch_library(KD_DN_LIB), kd_terms)

    print("querying L1000CDS2 for signature reversers", file=sys.stderr, flush=True)
    runs = [("consensus", up, dn)]
    for c in CELL_LINES:
        t = "{} knockdown 96h {}".format(GENE, c)
        if t in kd_up and t in kd_dn:
            runs.append((c, kd_up[t], kd_dn[t]))

    hits, per_run = [], {}
    for label, u, d in runs:
        rows = [r for r in query_reversers(u, d, label)
                if str(r.get("pert_desc")) not in ("None", LINCS_NULL, "")]
        per_run[label] = rows
        for r in rows:
            hits.append({"run": label, "compound": r.get("pert_desc"),
                         "pert_id": r.get("pert_id"), "score": r.get("score"),
                         "cell_id": r.get("cell_id"), "dose": r.get("pert_dose"),
                         "dose_unit": r.get("pert_dose_unit"),
                         "pubchem_id": r.get("pubchem_id"),
                         "compound_class": classify(r.get("pert_desc"))})
        time.sleep(0.5)

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    cols = ["run", "compound", "pert_id", "score", "cell_id", "dose", "dose_unit",
            "pubchem_id", "compound_class"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for h in sorted(hits, key=lambda x: (x["run"], -float(x["score"] or 0))):
            w.writerow(h)

    # ---- everything below is computed ----
    cell_runs = [l for l in per_run if l != "consensus"]
    recur = Counter()
    for l in cell_runs:
        for c in {str(r.get("pert_desc")) for r in per_run[l]}:
            recur[c] += 1

    up_ifn = sorted(set(up) & INTERFERON)
    up_p53 = sorted(set(up) & P53_TARGETS)
    dn_prol = sorted(set(dn) & PROLIFERATION)
    cls_counts = Counter(h["compound_class"] for h in hits)
    n_classified = sum(v for k, v in cls_counts.items() if k != "unclassified")
    liab = Counter(CONTRAINDICATED_BY[h["compound_class"]] for h in hits
                   if h["compound_class"] in CONTRAINDICATED_BY)
    n_contra = sum(liab.values())
    frac_contra = n_contra / len(hits) if hits else 0.0
    frac_unc = cls_counts.get("unclassified", 0) / len(hits) if hits else 0.0
    unclassified_note = len({h["compound"] for h in hits
                             if h["compound_class"] == "unclassified"})

    L = ["=" * 88, " T2-06 - SIGNATURE-BASED REPURPOSING (Connectivity Map)", "=" * 88, "",
         "t2_03 searched drug-TARGET databases and found nothing, and named this",
         "omission as the largest gap in the analysis. This stage closes it.",
         "",
         "Method: take the transcriptional consequence of losing {}, and look for".format(GENE),
         "compounds whose own signature REVERSES it. Requires no druggable target.",
         "",
         "## The query signature", "",
         "  source                 LINCS L1000 CRISPR KO consensus ({})".format(KO_LIB),
         "  UP genes               {}".format(len(up)),
         "  DOWN genes             {}".format(len(dn)),
         "  internal control       {} present in its own DOWN set: {}".format(GENE, self_down),
         "  cell-line replicates   {} of {} shRNA knockdown signatures retrieved".format(
             len(cell_runs), len(CELL_LINES)),
         "",
         "## What the signature actually contains", "",
         "  UP, type-I interferon genes    {:>3}   {}".format(
             len(up_ifn), ", ".join(up_ifn[:10]) + (" ..." if len(up_ifn) > 10 else "")),
         "  UP, p53 target genes           {:>3}   {}".format(
             len(up_p53), ", ".join(up_p53[:10]) + (" ..." if len(up_p53) > 10 else "")),
         "  DOWN, proliferation genes      {:>3}   {}".format(
             len(dn_prol), ", ".join(dn_prol[:10]) + (" ..." if len(dn_prol) > 10 else "")),
         ""]

    if up_ifn and up_p53:
        L += ["  This is mechanistically coherent. Missegregation produces micronuclei,",
              "  micronuclei activate cGAS-STING, and mitotic stress activates p53. The",
              "  signature is capturing real consequences of losing BubR1.",
              "",
              "  IT IS ALSO THE PROBLEM. Those {} interferon genes and {} p53 targets are".format(
                  len(up_ifn), len(up_p53)),
              "  the cell's DEFENCE against aneuploidy, not the lesion. Connectivity Map",
              "  assumes reversing a disease signature is therapeutic. Reversing THIS",
              "  signature means suppressing interferon and p21 and pushing proliferation",
              "  genes back up - in a child with a cancer predisposition syndrome.",
              ""]

    L += ["## What the search returned", "",
          "  total hits across {} runs      {}".format(len(runs), len(hits)),
          "  distinct compounds             {}".format(len({h["compound"] for h in hits})),
          "", "  by compound class:", ""]
    L.append("    {:<36}{:>5}{:>7}   {}".format("CLASS", "N", "%", "LIABILITY HERE"))
    for cls, n in cls_counts.most_common():
        L.append("    {:<36}{:>5}{:>7.0%}   {}".format(
            cls, n, n / len(hits) if hits else 0, CONTRAINDICATED_BY.get(cls, "-")))

    L += ["",
          "  classified                      {:>4}  ({:.0%})".format(n_classified, 1 - frac_unc),
          "  CONTRAINDICATED for this child  {:>4}  ({:.0%})".format(n_contra, frac_contra),
          "     antiproliferative            {:>4}".format(liab.get("antiproliferative", 0)),
          "     immunosuppressive            {:>4}".format(liab.get("immunosuppressive", 0)),
          "  unclassified                    {:>4}  ({:.0%})".format(
              cls_counts.get("unclassified", 0), frac_unc),
          ""]

    if recur:
        L += ["## Reproducibility across independent cell lines", "",
              "  A compound recovered in one cell line is noise. These are the",
              "  compounds recovered in the most of the {} shRNA knockdown runs:".format(len(cell_runs)),
              "",
              "    {:<32}{:>8}  {}".format("COMPOUND", "RUNS", "CLASS")]
        for c, n in recur.most_common(15):
            L.append("    {:<32}{:>4}/{:<3}  {}".format(c[:31], n, len(cell_runs), classify(c)))
        top_n = recur.most_common(1)[0][1] if recur else 0
        L += ["",
              "  Best reproducibility: {} of {} cell lines.".format(top_n, len(cell_runs))]

    L += ["", "=" * 88, " VERDICT", "=" * 88, ""]

    if frac_contra >= 0.5:
        L += ["  {:.0%} of the hits ({} of {}) belong to a class that is contraindicated".format(
                  frac_contra, n_contra, len(hits)),
              "  for this patient, on one of two grounds computed above:",
              "",
              "    {:>4} antiproliferative  - they act on the axis the disease already".format(
                  liab.get("antiproliferative", 0)),
              "         compromises, and several are checkpoint inhibitors outright",
              "    {:>4} immunosuppressive  - the same liability that demoted amlexanox".format(
                  liab.get("immunosuppressive", 0)),
              "         against HP:0002859 in t2_04, arrived at independently here",
              "",
              "  Connectivity Map returns exactly what it is designed to return, and",
              "  what it returns is unusable. The reason is mechanistic rather than",
              "  incidental: {} of the UP genes are type-I interferon and {} are p53".format(
                  len(up_ifn), len(up_p53)),
              "  targets, so a large part of this 'disease signature' is the cell's",
              "  DEFENCE against aneuploidy. Reversing it means suppressing that",
              "  defence, and the most efficient way to reverse a mitotic-stress",
              "  response is to arrest or kill the cell.",
              "",
              "  CMap's founding assumption - that reversing a disease signature is",
              "  therapeutic - does not hold for a chromosomal instability syndrome."]
    else:
        L += ["  Only {:.0%} of hits ({} of {}) fall into a class this analysis can call".format(
                  frac_contra, n_contra, len(hits)),
              "  contraindicated, and {:.0%} remain unclassified. That is too little".format(frac_unc),
              "  coverage to characterise the result set as a whole, and no claim about",
              "  what CMap returned here should be made beyond the per-class counts",
              "  printed above."]

    if unclassified_note:
        L += ["",
              "  {} distinct compounds remain unclassified. They are listed in the TSV.".format(
                  unclassified_note),
              "  The taxonomy is keyword-based and therefore incomplete by construction;",
              "  the percentages above are lower bounds, not exact partitions."]

    vinc = sorted({h["compound"] for h in hits if "vincrist" in (h["compound"] or "").lower()})
    if vinc:
        L += ["",
              "  Note that {} appears among the hits. It is already in this".format(vinc[0]),
              "  child's rhabdomyosarcoma protocol - a spindle poison given to a patient",
              "  whose spindle checkpoint is deficient. The method surfacing it as a",
              "  'hit' shows why signature reversal cannot be read as a therapeutic",
              "  recommendation without a direction-of-effect argument."]

    L += ["",
          "  THIS IS THE THIRD INDEPENDENT SEARCH DIRECTION TO CONVERGE:",
          "",
          "    t2_03  drug-target databases   -> 0 compensatory, all inhibitors",
          "    t2_04  variant-class branch    -> 1 proposable, small expected effect",
          "    t2_06  signature reversal      -> {:.0%} of hits contraindicated by class".format(frac_contra),
          "",
          "  None of the three finds an approved drug that helps. They fail for",
          "  different reasons, which is what makes the convergence informative",
          "  rather than repetitive.",
          "",
          "## Limitations of this stage specifically", "",
          "  - LINCS knockdown in immortalised cancer cell lines is not a child's",
          "    biallelic hypomorphic state. Most L1000 lines are p53-mutant and",
          "    already aneuploid, which is precisely the axis under study.",
          "  - Signature reversal identifies compounds that oppose a transcriptional",
          "    CONSEQUENCE. It cannot identify a compound that restores BubR1, which",
          "    is what the disease actually requires.",
          "  - CMap scores are similarity measures, not efficacy estimates.",
          "  - A 96 h knockdown models chronic partial loss poorly.",
          "  - The class taxonomy declares no BENEFICIAL category, because no approved",
          "    compound class restores a spindle checkpoint. The absence of a helpful",
          "    hit here is therefore partly by construction; that conclusion rests on",
          "    t2_03, which searched drug-target space directly, not on this stage.",
          "",
          "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
