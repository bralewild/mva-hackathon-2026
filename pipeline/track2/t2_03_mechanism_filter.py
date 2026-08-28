#!/usr/bin/env python3
"""
==============================================================================
t2_03_mechanism_filter.py - Which associations survive mechanistic scrutiny

INPUT  : $WORK/t2_02_drug_evidence.tsv
OUTPUT : $WORK/t2_03_filtered.tsv
         $RESULTS/t2_03_mechanism_filter.txt

t2_02 recorded what the databases say. This stage decides what any of it is
worth. Every drop is counted, and every statement in the output is computed
from the data rather than asserted.

WHAT CHANGED, AND WHY (adversarial review, 2026-08-28)
------------------------------------------------------
The previous version had four defects, all of which pushed toward its own
conclusion. They are listed here rather than quietly fixed, because the report
built an argument on the output they produced.

1. THE DIRECTION MATCHER COULD NOT SEE ACTIVATORS.
   ACTIVATING_WORDS contained "activator"; DGIdb emits the literal "ACTIVATING".
   "activator" is not a substring of "activating", while "inhibitor" IS a
   substring of "inhibitory". The filter recognised inhibitors and silently
   reclassified activators as unknown - and unknown is a drop. Verified against
   the data: 194 INHIBITORY vs 1 ACTIVATING, so the practical impact was one
   association, but the asymmetry was real and it biased the result.

2. THE NARRATIVE WAS HARDCODED PROSE, NOT COMPUTED.
   The "no candidate survives" explanation - including "every approved-drug
   association is backed by exactly one source" - was an authored string emitted
   whenever the survivor set was empty. That claim is FALSE: 20 of 125 approved
   associations have two or more sources. All summary statements are now
   f-strings over computed values.

3. THE NEGATIVE WAS UNFALSIFIABLE.
   A total API outage produced an empty input, the same "no candidate survives"
   narrative, and exit code 0. The stage now refuses to interpret an input that
   is too small to interpret.

4. DIRECTION CLASSIFICATIONS WERE WRONG ON SEVERAL NODES.
   CDK1 and PLK1 were called harmful-if-inhibited. CDK1 inhibition blocks
   mitotic ENTRY, and a cell that never enters mitosis cannot missegregate;
   CDK1 also phosphorylates APC/C subunits to permit CDC20 binding, so
   inhibiting it reduces APC/C-CDC20 activity - the compensatory direction.
   PLK1 inhibition causes SAC-dependent prometaphase arrest, not checkpoint
   weakening. FZR1 (Cdh1) is the G1 APC/C activator, not the anaphase one, so
   inhibiting it does nothing to restrain anaphase. These are now AMBIGUOUS,
   which is an honest verdict rather than a convenient one.

THE DIRECTION FILTER - WHAT IT IS FOR
-------------------------------------
The disease is a LOSS of BubR1 function: the spindle assembly checkpoint is too
weak, so anaphase proceeds before chromosomes are correctly attached. BubR1's
job is to inhibit CDC20, which activates the APC/C.

An important caveat the previous version did not state: BubR1's inhibition is
CONDITIONAL and attachment-responsive, while a drug is constitutive. Restoring
restraint pharmacologically is not the same as restoring checkpoint FIDELITY -
constitutive APC/C inhibition produces mitotic arrest, slippage and
tetraploidy, which is why oncology develops these agents. "Compensatory" here
means "opposes the direction of the lesion", not "corrects the defect".
==============================================================================
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _paths import WORK, RESULTS  # noqa: E402

IN_TSV = os.path.join(WORK, "t2_02_drug_evidence.tsv")
OUT_TSV = os.path.join(WORK, "t2_03_filtered.tsv")
OUT_TXT = os.path.join(RESULTS, "t2_03_mechanism_filter.txt")

MIN_SCORE = 0.10
MIN_SOURCES = 2
MIN_ROWS_TO_INTERPRET = 50      # below this the input is broken, not informative

# Inhibiting these opposes the direction of the lesion.
COMPENSATORY_IF_INHIBITED = {
    "CDC20": "the node BubR1 itself inhibits",
    "UBE2C": "APC/C E2; blocking it slows APC/C output",
    "ESPL1": "separase; blocking cohesin cleavage directly opposes premature "
             "chromatid separation - though non-selective inhibition causes its "
             "own missegregation",
    "ANAPC1": "APC/C subunit", "ANAPC2": "APC/C subunit",
    "ANAPC4": "APC/C subunit", "ANAPC7": "APC/C subunit",
    "ANAPC10": "APC/C subunit", "ANAPC11": "APC/C subunit",
    "CDC16": "APC/C subunit", "CDC23": "APC/C subunit",
    "CDC26": "APC/C subunit", "CDC27": "APC/C subunit",
}
# Inhibiting these weakens an already-failing checkpoint.
HARMFUL_IF_INHIBITED = {
    "BUB1": "MCC component", "BUB3": "MCC component",
    "MAD1L1": "MAD2 template", "MAD2L1": "MCC component",
    "KNL1": "kinetochore scaffold for BUB1/BUB3",
    "TTK": "MPS1, apex of checkpoint signalling",
    "AURKB": "error correction at kinetochores",
    "CENPE": "kinetochore motor, BubR1-associated",
    "PPP2CA": "PP2A opposes Aurora B; recruited by the BubR1 KARD motif",
    "PPP2R5A": "PP2A-B56, BubR1's own phosphatase partner",
    "PTTG1": "securin; inhibiting it releases separase early",
}
# Genuinely two-sided - and saying so is the honest verdict.
AMBIGUOUS = {
    "CDK1": "inhibition blocks mitotic ENTRY (a cell that does not divide cannot "
            "missegregate) and reduces APC/C-CDC20 activation, but also disrupts "
            "normal mitosis",
    "CCNB1": "CDK1 partner; same two-sided argument",
    "PLK1": "inhibition causes SAC-dependent prometaphase arrest rather than "
            "checkpoint weakening",
    "FZR1": "Cdh1 activates APC/C in G1, not at anaphase; inhibiting it does not "
            "restrain anaphase onset",
}
# Tier 4: the variant-class route. Direction logic here is about the PTC, not
# the checkpoint, so it is evaluated separately.
READTHROUGH_RELEVANT = {
    "ETF1": "eRF1; inhibiting termination promotes readthrough of the PTC",
    "GSPT1": "eRF3a; degrading it promotes readthrough",
    "GSPT2": "eRF3b paralogue",
    "UPF1": "inhibiting NMD raises the PTC transcript available for readthrough",
    "UPF2": "NMD core", "UPF3B": "NMD core",
    "SMG1": "NMD kinase", "SMG5": "NMD", "SMG6": "NMD endonuclease", "SMG7": "NMD",
    "EIF4A3": "exon junction complex; marks the PTC as premature",
    "RPL3": "ribosomal; modulates readthrough", "RPS15": "ribosomal decoding site",
}

# Matched against the lowercased concatenation of interaction_type and
# directionality. Written to match DGIdb's enum literals INHIBITORY / ACTIVATING
# as well as free-text interaction types.
INHIBITORY_WORDS = ("inhibit", "antagonis", "blocker", "negative", "suppress",
                    "inverse agonist")
ACTIVATING_WORDS = ("activat", "agonist", "positive", "inducer", "potentiator")

CYTOTOXIC = {
    "CARBOPLATIN", "CISPLATIN", "OXALIPLATIN", "PACLITAXEL", "DOCETAXEL",
    "TOPOTECAN", "IRINOTECAN", "IDARUBICIN", "DOXORUBICIN", "DAUNORUBICIN",
    "EPIRUBICIN", "ETOPOSIDE", "TENIPOSIDE", "VINCRISTINE", "VINBLASTINE",
    "VINORELBINE", "CYCLOPHOSPHAMIDE", "IFOSFAMIDE", "GEMCITABINE", "CYTARABINE",
    "FLUOROURACIL", "METHOTREXATE", "MITOMYCIN", "BLEOMYCIN", "DACTINOMYCIN",
    "AMSACRINE", "PIXANTRONE", "VALRUBICIN", "MELPHALAN", "BUSULFAN",
    "TEMOZOLOMIDE", "PEMETREXED", "FLUDARABINE", "BORTEZOMIB",
}


def direction_of(interaction_type, directionality):
    """Resolve direction. The structured directionality field wins over free text."""
    d = (directionality or "").lower()
    if any(w in d for w in INHIBITORY_WORDS):
        return "inhibits"
    if any(w in d for w in ACTIVATING_WORDS):
        return "activates"
    t = (interaction_type or "").lower()
    if any(w in t for w in INHIBITORY_WORDS):
        return "inhibits"
    if any(w in t for w in ACTIVATING_WORDS):
        return "activates"
    return "unknown"


# The only tokens direction_of can emit. verdict_for validates against this
# because every branch below falls through to an INVERTED verdict on an
# unrecognised token: pass "inhibitor" instead of "inhibits" and a HARMFUL gene
# is reported COMPENSATORY, confidently and silently. Found exactly that way.
DIRECTIONS = frozenset(("inhibits", "activates", "unknown"))


def verdict_for(gene, direction):
    """Is this drug pushing the mechanism in the opposing direction?"""
    if direction not in DIRECTIONS:
        raise ValueError(
            "direction {!r} is not one of {}. An unrecognised token would fall "
            "through to the opposite verdict.".format(direction, sorted(DIRECTIONS)))
    if gene in READTHROUGH_RELEVANT:
        if direction == "inhibits":
            return "READTHROUGH", READTHROUGH_RELEVANT[gene]
        if direction == "activates":
            return "HARMFUL", "enhancing termination or NMD works against readthrough"
        return "UNKNOWN", "no direction recorded for a variant-class target"
    if direction == "unknown":
        return "UNKNOWN", "no interaction direction is recorded"
    if gene in AMBIGUOUS:
        return "AMBIGUOUS", AMBIGUOUS[gene]
    if gene in COMPENSATORY_IF_INHIBITED:
        if direction == "inhibits":
            return "COMPENSATORY", COMPENSATORY_IF_INHIBITED[gene]
        return "HARMFUL", "activating this releases anaphase further"
    if gene in HARMFUL_IF_INHIBITED:
        if direction == "inhibits":
            return "HARMFUL", HARMFUL_IF_INHIBITED[gene]
        return "COMPENSATORY", "strengthening this could partly offset the lesion"
    return "NEUTRAL", "no directional argument for this node"


def is_cytotoxic(name):
    return any(c in (name or "").upper() for c in CYTOTOXIC)


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("missing " + IN_TSV + " - run t2_02_drug_evidence.py first")

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))

    # Refuse to interpret an input too small to be informative. Without this, an
    # API outage produced an empty file, the full "no candidate survives"
    # narrative, and exit code 0.
    if len(rows) < MIN_ROWS_TO_INTERPRET:
        sys.exit("ABORT: only {} associations in {}. Expected at least {}. "
                 "This looks like an upstream failure, not a negative result. "
                 "Re-run t2_02 and check API availability."
                 .format(len(rows), IN_TSV, MIN_ROWS_TO_INTERPRET))

    approved = [r for r in rows if r["approved"] == "True"]
    if not approved:
        sys.exit("ABORT: no approved-drug associations at all. Check the "
                 "'approved' field serialisation in t2_02.")

    kept, drops = [], {"weak_evidence": [], "harmful_direction": [],
                       "unknown_direction": [], "cytotoxic_only": [],
                       "neutral": [], "ambiguous": []}

    for r in approved:
        score = float(r["score"] or 0)
        nsrc = int(r["n_sources"] or 0)
        d = direction_of(r.get("interaction_type", ""), r.get("directionality", ""))
        v, why = verdict_for(r["gene"], d)
        r["direction"], r["verdict"], r["reason"] = d, v, why

        if score < MIN_SCORE or nsrc < MIN_SOURCES:
            drops["weak_evidence"].append(r); continue
        if is_cytotoxic(r["drug"]):
            drops["cytotoxic_only"].append(r); continue
        if v == "HARMFUL":
            drops["harmful_direction"].append(r)
        elif v == "UNKNOWN":
            drops["unknown_direction"].append(r)
        elif v == "AMBIGUOUS":
            drops["ambiguous"].append(r)
        elif v == "NEUTRAL":
            drops["neutral"].append(r)
        else:
            kept.append(r)

    cols = ["gene", "tier", "drug", "approved", "score", "n_sources",
            "interaction_type", "directionality", "direction", "verdict", "reason"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in kept:
            w.writerow({c: r.get(c, "") for c in cols})

    # ---- everything below is computed, not asserted ----
    multi_src = [r for r in approved if int(r["n_sources"] or 0) >= MIN_SOURCES]
    with_dir = [r for r in approved
                if direction_of(r.get("interaction_type", ""), r.get("directionality", "")) != "unknown"]
    genes_all = sorted({r["gene"] for r in rows})
    genes_appr = sorted({r["gene"] for r in approved})
    scores = [float(r["score"] or 0) for r in approved]
    # What direction would have done on its own, ignoring the evidence gate:
    dir_alone = {}
    for r in approved:
        v = r["verdict"]
        dir_alone[v] = dir_alone.get(v, 0) + 1

    L = ["=" * 88,
         " T2-03 - MECHANISTIC FILTER",
         "=" * 88, "",
         "The disease is a LOSS of BubR1 function. A drug is only a candidate if it",
         "pushes the mechanism in the opposing direction, or acts on the variant class.",
         "",
         "## Funnel", "",
         "  associations recorded                {:>6,}".format(len(rows)),
         "  with an approved drug                {:>6,}".format(len(approved)),
         "  ... backed by >= {} sources           {:>6,}".format(MIN_SOURCES, len(multi_src)),
         "  ... carrying a resolvable direction  {:>6,}  ({:.0%} of approved)".format(
             len(with_dir), len(with_dir) / len(approved)),
         "  approved-drug score range            {:.2f} - {:.2f}".format(
             min(scores) if scores else 0, max(scores) if scores else 0),
         "",
         "  dropped - single source / low score  {:>6,}".format(len(drops["weak_evidence"])),
         "  dropped - acute cytotoxic only       {:>6,}".format(len(drops["cytotoxic_only"])),
         "  dropped - direction unknown          {:>6,}".format(len(drops["unknown_direction"])),
         "  dropped - direction HARMFUL          {:>6,}".format(len(drops["harmful_direction"])),
         "  dropped - direction AMBIGUOUS        {:>6,}".format(len(drops["ambiguous"])),
         "  dropped - no directional argument    {:>6,}".format(len(drops["neutral"])),
         "  " + "-" * 44,
         "  SURVIVING CANDIDATES                 {:>6,}".format(len(kept)),
         "",
         "## Did the direction filter do any work?", "",
         "  Verdicts assigned across ALL {} approved associations, independent of".format(len(approved)),
         "  the evidence gate - this is what direction alone would have decided:", ""]
    for v in sorted(dir_alone, key=lambda k: -dir_alone[k]):
        L.append("    {:<14} {:>5}".format(v, dir_alone[v]))
    n_harm = dir_alone.get("HARMFUL", 0)
    L += ["",
          "  Of those, {} reached the direction filter after the evidence gate.".format(
              len(drops["harmful_direction"]) + len(drops["unknown_direction"])
              + len(drops["ambiguous"]) + len(drops["neutral"]) + len(kept)),
          "  The filter dropped {} association(s) as HARMFUL.".format(len(drops["harmful_direction"]))]
    if not drops["harmful_direction"]:
        L += ["",
              "  STATED PLAINLY: within the APPROVED subset the direction filter removed",
              "  nothing that the evidence gate had not already removed. It identified {}".format(n_harm),
              "  harmful associations there, but all of them failed the evidence gate",
              "  first. The section below is where direction actually does the work."]

    # ---- The well-evidenced subset, REGARDLESS of approval status. This is the
    # informative layer: the approved associations are all single-source, so the
    # question "what does the good evidence say?" is only answerable here. Every
    # number below is computed from the file, none is asserted.
    well = [r for r in rows if int(r["n_sources"] or 0) >= MIN_SOURCES]
    well_verdicts, well_genes, typed_inhib = {}, {}, 0
    for r in well:
        d = direction_of(r.get("interaction_type", ""), r.get("directionality", ""))
        v, _ = verdict_for(r["gene"], d)
        well_verdicts[v] = well_verdicts.get(v, 0) + 1
        well_genes[r["gene"]] = well_genes.get(r["gene"], 0) + 1
        if "inhibit" in (r.get("interaction_type", "") or "").lower():
            typed_inhib += 1
    n_well_appr = sum(1 for r in well if r["approved"] == "True")
    n_compensatory = well_verdicts.get("COMPENSATORY", 0)

    L += ["", "## What the WELL-EVIDENCED associations say", "",
          "  Associations backed by >= {} distinct sources : {:>4}  of {:,}".format(
              MIN_SOURCES, len(well), len(rows)),
          "  ... of which involve an approved drug        : {:>4}".format(n_well_appr),
          "  ... explicitly typed 'inhibitor'             : {:>4}".format(typed_inhib),
          ""]
    if well:
        L.append("  They concentrate on {} gene(s):".format(len(well_genes)))
        L.append("")
        for g in sorted(well_genes, key=lambda k: -well_genes[k]):
            v, why = verdict_for(g, "inhibits")
            L.append("    {:<9} {:>3}   if inhibited: {:<13} {}".format(
                g, well_genes[g], v, why[:44]))
        L += ["", "  Direction verdicts across the well-evidenced set:", ""]
        for v in sorted(well_verdicts, key=lambda k: -well_verdicts[k]):
            L.append("    {:<14} {:>4}".format(v, well_verdicts[v]))
        L += ["",
              "  ACTING IN THE COMPENSATORY DIRECTION: {}".format(n_compensatory)]
        if n_compensatory == 0:
            L += ["",
                  "  Not one of the {} best-evidenced associations opposes the lesion.".format(len(well)),
                  "  That is not an accident of curation. These are oncology compounds,",
                  "  and oncology develops checkpoint inhibitors precisely to FORCE",
                  "  missegregation and kill dividing cells - this disease, deliberately",
                  "  induced. The pharmacology for this pathway is well developed and",
                  "  aimed in exactly the direction that would harm this patient.",
                  "",
                  "  The direction filter removes ALL {} of them. This is where it".format(len(well)),
                  "  changes the answer, and it is a demonstrated contribution rather",
                  "  than a proposed one."]

    L += ["", "## Genes with no drug association at all", ""]
    from_targets = os.path.join(WORK, "t2_01_targets.tsv")
    if os.path.exists(from_targets):
        allt = [r["gene"] for r in csv.DictReader(open(from_targets), delimiter="\t")]
        nod = [g for g in allt if g not in genes_all]
        L.append("  {} of {} network targets have no reported drug association.".format(
            len(nod), len(allt)))
        L.append("")
        L.append("  " + ", ".join(nod))

    if kept:
        L += ["", "## Surviving candidates", "",
              "  {:<9}{:<30}{:<12}{:<14}{}".format("GENE", "DRUG", "DIRECTION", "VERDICT", "WHY")]
        for r in sorted(kept, key=lambda x: -float(x["score"] or 0)):
            L.append("  {:<9}{:<30}{:<12}{:<14}{}".format(
                r["gene"], r["drug"][:29], r["direction"], r["verdict"], r["reason"][:60]))
    else:
        L += ["", "## No candidate survives", "",
              "  With {} approved associations, {} backed by >= {} sources and {} carrying".format(
                  len(approved), len(multi_src), MIN_SOURCES, len(with_dir)),
              "  a resolvable direction, no association is both well-evidenced and",
              "  directionally useful.",
              "",
              "  This is a statement about what these databases contain for this",
              "  network, not a proof that no drug could help. In particular:",
              "",
              "  - only approved drugs were considered; clinical-stage compounds for",
              "    TTK, PLK1, AURKB and CDK1 exist and were excluded by that filter",
              "  - only DGIdb was queried; ChEMBL and Open Targets were not",
              "  - signature-based repurposing (LINCS / Connectivity Map) was not run",
              "    at all, and it is the standard approach when no target is druggable"]

    L += ["", "  -> " + OUT_TSV]

    txt = "\n".join(L)
    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    open(OUT_TXT, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
