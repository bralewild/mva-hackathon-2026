# Track 2 — Drug Repositioning for biallelic *BUB1B*

**Participant:** `bralewild` (individual)
**Builds on Track 1:** biallelic *BUB1B* — `NM_001211.6:c.2210T>G` (p.Leu737Ter)
and `c.3006T>G` (p.Asn1002Lys), scored 100/100, F-max 1.000
**Code:** `pipeline/track2/` · evidence tables in `results/track2_evidence/`
**Date:** 2026-08-28

---

## 1. The proposal

**Escin** (β-aescin), a marketed triterpene saponin, is the single candidate that
is simultaneously (a) evidenced for translational readthrough, (b) currently
available as a medicine, and (c) free of any liability contraindicated by this
child's own phenotype.

Its readthrough activity was not a literature guess: escin was identified in an
**unbiased high-throughput screen of ~1,600 clinically approved compounds**
against CFTR premature termination codons, where it proved *"efficacious and
potent in a variety of primary human airway cells"* carrying G542X and W1282X.
It is marketed under a German Commission E monograph (1984, renewed 1994) and an
EMA traditional-use registration, with established human dosing of 50–75 mg twice
daily.

That ranking is produced by the pipeline, not asserted. Two better-known
candidates were demoted by an automated screen against the proband's HPO terms:

| Candidate | Mechanistic rank | After screen | Why |
|---|---:|---:|---|
| **Escin** | 2 | **1** | marketed; no phenotype conflict |
| Gentamicin / amikacin | 3 | 2 | **nephrotoxic** vs `HP:0000121` nephrocalcinosis |
| Ataluren | 4 | 3 | not marketed — reference only |
| Amlexanox | **1** | 4 | not marketed; **TBK1/IKKε immunosuppression** vs `HP:0002859` rhabdomyosarcoma |

**Amlexanox is the best mechanistic fit and it is not proposable.** It is the
only agent reported to induce readthrough *and* inhibit NMD — the pair of
activities this lesion needs — but it has no marketed product anywhere
(Aphthasol discontinued in the US; Solfa discontinued by Takeda in Japan in
2019), and its principal systemic pharmacology is TBK1/IKKε inhibition, which is
a poor thing to give a child with an active cancer predisposition syndrome.

**Expected effect size is small, and we say so up front.** Published in-cell
readthrough efficiencies for this drug class are single-digit percent. This is
not a proposal that a threshold will be crossed. It is the observation that
readthrough is the only pharmacological route that addresses the causal allele at
all, together with the experiment that would falsify it (§6).

---

## 2. Mechanism

### 2.1 Direction of effect: loss of function

*BUB1B* encodes BubR1, a core component of the mitotic checkpoint complex (MCC).
Its essential role is **stoichiometric**: BubR1, BUB3 and MAD2L1 sequester CDC20
and prevent it from activating the anaphase-promoting complex until every
chromosome is correctly attached.

A correction worth making explicitly, because an earlier draft of this analysis
got it wrong: **human BubR1 is widely characterised as a pseudokinase**, with no
demonstrable catalytic activity, and its kinase domain is dispensable for
checkpoint function. Attributing p.Asn1002Lys's effect to "reduced enzymatic
activity" would be inventing a mechanism. What can be said is that it is a
full-length protein carrying a substitution in that domain, of **unknown**
functional consequence.

| Allele | Change | What is established | What is not |
|---|---|---|---|
| 1 | `c.2210T>G` p.Leu737Ter | nonsense, exon 17/23, NMD-predicted, ClinVar Pathogenic/LP | — |
| 2 | `c.3006T>G` p.Asn1002Lys | absent from gnomAD, CADD 24.5; the same amino-acid change via `c.3006T>A` is a **ClinVar VUS** | **residual function — entirely unmeasured** |

```
BubR1 ↓ → MCC cannot restrain CDC20 → APC/C activates prematurely
      → premature sister chromatid separation
      → chromosome missegregation → mosaic aneuploidy
      → chromosomal instability → cancer predisposition
      → rhabdomyosarcoma (HP:0002859), the presenting event
```

### 2.2 The dose–response, and what it does and does not permit

BubR1 reduction produces a **graded** phenotype. In HeLa cells carrying an shRNA
gradient of residual expression, *all* transduced cells showed premature
chromatid separation, its severity correlated inversely with residual
expression, and detectable aneuploidy appeared below roughly 50 % (*Gradual
reduction of BUBR1 protein levels results in premature sister-chromatid
separation then in aneuploidy*, **Hum Genet 2008**, PMID 18932004).

An earlier draft treated 50 % as a **threshold to be crossed**. That framing does
not survive arithmetic, and the correction matters:

> Normalise each allele to 50 points. Allele 1 contributes ≈ 0. Therefore
> **a patient with one null allele cannot exceed 50 % even if the other allele
> were perfect.** The ceiling *is* the figure. Obligate carriers sit at exactly
> 50 % and are healthy — so "become a healthy carrier" is the maximum a perfect
> therapy could achieve, and single-digit readthrough moves perhaps 2–3 points
> toward it.

What survives is the **gradient**, not the cliff. Severity of premature chromatid
separation tracks residual expression continuously, so a small increase predicts
a small improvement — measurable in cells, not obviously meaningful in a child.
That is a weaker claim than the one this analysis started with, and it is the one
the data supports.

Two caveats on transferring even the gradient: HeLa is p53-null, hypotriploid and
already chromosomally unstable, so it is a poor reporter for a threshold in a
p53-competent child; and an shRNA gradient reduces *wild-type* protein, whereas
this patient has a full-length mutant that still occupies its slot in a
stoichiometric complex. A defective subunit that fills its position can be worse
than an absent one. Neither confound is resolvable from published data.

### 2.3 We looked for the cellular phenotype, and the result is ambiguous

Track 1 stage 08 screened for mosaic aneuploidy from B-allele fractions at
2,241,007 heterozygous SNVs across 285 windows.

An initial pass flagged eight autosomes, chr21 at an apparent 17 % mosaic
fraction. A uniformity test rejected it: a true whole-chromosome event affects
every window, and the signal concentrated in a minority — on the acrocentric and
heterochromatic chromosomes where mapping is hardest. But the rejection is not
clean, and the report says so:

- the uniformity test for chr21 ran on **5 windows**, which is little power;
- the GC-content confounder check returned **r = +0.475**, which the pipeline
  itself reports as *"GC bias contributes and cannot be excluded"*;
- MVA aneuploidy is **variegated** — different chromosomes in different cells —
  so a uniform whole-chromosome model may be the wrong detector entirely.

**The honest statement is that this screen is inconclusive**, not that
mosaicism is absent. An earlier draft also argued that blood is the wrong tissue
because aneuploid lymphocytes are "enriched in culture"; that is
self-contradictory, and MVA is in fact routinely diagnosed by karyotype on
cultured peripheral blood lymphocytes. Ex vivo work should use fibroblasts for
reasons of assay control, not because blood is uninformative.

---

## 3. The systematic search

Before proposing anything, we asked computationally whether a drug already acts
on this pathway. `pipeline/track2/` mirrors the Track 1 discipline: numbered
stages, traceable evidence, every discard counted, evidence tables committed.

### 3.1 Network (t2-01) — and a correction

52 targets, tiered by mechanistic distance:

| Tier | Members | Rationale |
|---|---|---|
| 0 | BUB1B | the mutated gene |
| 1 | MCC and CDC20 | what BubR1 acts in and on |
| 2 | APC/C, separase, securin, checkpoint kinases, PP2A | what the MCC restrains |
| 3 | additional physical STRING interactors | the immediate complex |
| **4** | **ETF1, GSPT1/2, UPF1/2/3B, SMG1/5/6/7, EIF4A3, RPL3, RPS15** | **actionable by variant CLASS, not by gene** |

Two corrections from an adversarial review of the first version:

**The network was contaminated.** Using STRING's *combined* score at ≥ 0.9 pulled
in BRCA2, TOP2A and BIRC5 — because mitotic genes are co-expressed and co-cited,
so the score saturates on those channels. Those three genes then supplied the
only well-evidenced drugs in the entire result: olaparib, etoposide,
trastuzumab. PARP inhibitors and topoisomerase poisons have nothing to do with
BubR1 dosage. Switching to the **physical** subnetwork removed them, and with
them 244 spurious associations.

**The search was blind to its own answer.** The proposal acts on the ribosome,
and no node of the translation-termination or NMD machinery was in the network —
so the search could not possibly have found it. A negative from a search that
cannot see the class of answer being proposed measures its own blind spot. Tier 4
fixes that, and it functions as a **positive control**: the corrected search now
independently recovers **ataluren** (on RPL3 and RPS15) and **ELX-02** from the
databases. The machinery works. Neither drug is marketed.

### 3.2 Evidence (t2-02) and filter (t2-03)

| | |
|---|---:|
| Associations returned | **424** |
| With an approved drug | **41** |
| Approved **and** backed by ≥ 2 distinct sources | **0** |
| Approved and carrying a resolvable interaction direction | **2 (5 %)** |
| Network targets with no association at all | **38 of 52** |
| **Surviving candidates** | **0** |

Every one of the 41 approved associations rests on a **single** source, so all
41 fail an evidence gate of score ≥ 0.10 and ≥ 2 distinct sources. Within the
approved subset, therefore, the direction filter removed nothing the evidence
gate had not already removed, and the pipeline says so itself rather than
claiming a decisiveness it did not have. An earlier draft claimed the opposite in
three separate sections; it was wrong, and the pipeline's own output contradicted
it.

### 3.3 The well-evidenced pharmacology exists — and all of it points the wrong way

The interesting result is one layer down. **43 of the 424 associations *are*
backed by two or more distinct sources** — none of them approved drugs. They are
not scattered across the network. They collapse onto four genes:

| Gene | Multi-source associations | Direction if inhibited |
|---|---:|---|
| AURKB | 21 | **harmful** |
| CDK1 | 13 | ambiguous |
| PLK1 | 8 | ambiguous |
| TTK | 1 | **harmful** |

**41 of the 43 are explicitly typed `inhibitor`.** The remaining two — rigosertib
and alisertib — carry no interaction type in DGIdb, so the pipeline refuses to
classify them. The verdicts it does assign are **21 harmful, 20 ambiguous, 2
unclassifiable**, and the number that matters is the one that is absent:

> **Compensatory: 0.** Not one of the 43 best-evidenced associations in the
> entire result opposes the lesion.

There is nothing accidental about that distribution. These are oncology
clinical-stage compounds, and oncology develops checkpoint inhibitors **precisely
to force missegregation and kill dividing cells** — which is the disease
mechanism here, deliberately induced. The pharmacology for this pathway is not
missing. It is well developed, well evidenced, and aimed in exactly the direction
that would harm this patient.

That is a stronger and more useful negative than an empty database, and it is
where the direction filter earns its place. All 43 are removed — 21 because
inhibiting AURKB or TTK would weaken an already-failing checkpoint, 20 because
CDK1 and PLK1 inhibition is genuinely two-sided (§3.4), and 2 because DGIdb
records no direction at all. Three different reasons, none of them
"insufficient evidence".

### 3.4 The direction assignments, and a correction

| Node | Effect of inhibiting | Verdict |
|---|---|---|
| CDC20, APC/C subunits, separase | restrains anaphase — opposes the lesion | compensatory |
| BUB1, BUB3, MAD2L1, KNL1, TTK, AURKB, CENPE, PP2A | weakens an already-failing checkpoint | harmful |
| **CDK1, CCNB1, PLK1, FZR1** | **two-sided** | **ambiguous** |

The ambiguous row is a correction to our own first version, which called CDK1 and
PLK1 harmful-if-inhibited. CDK1 inhibition blocks mitotic *entry* — a cell that
does not divide cannot missegregate — and CDK1 also phosphorylates APC/C subunits
to permit CDC20 binding, so inhibiting it reduces APC/C–CDC20 activity, the
compensatory direction. PLK1 inhibition causes SAC-dependent prometaphase arrest,
not checkpoint weakening. FZR1 activates APC/C in G1, not at anaphase. **Twenty
of the 43 land in that row**, so the ambiguity is not a technicality — it is
nearly half the well-evidenced result, and calling it harmful would have been
easier and wrong.

And a caveat the first version omitted: **BubR1's inhibition of CDC20 is
conditional and attachment-responsive; a drug is constitutive.** Constitutive
APC/C inhibition does not restore checkpoint *fidelity* — it produces mitotic
arrest, slippage and tetraploidy, which is why oncology develops these agents.
"Compensatory" here means *opposes the direction of the lesion*, not *corrects
the defect*. The distinction is the difference between a sign and a mechanism.

### 3.5 Is the negative real? A positive control

A zero is worthless on its own. A pipeline strict enough to return zero for
BUB1B might return zero for **every** gene — in which case the number measures
the filters, not the biology. `t2_05_positive_control.py` tests exactly that:
the same evidence gate, applied to loss-of-function rare diseases where an
approved drug does act on the deficient product.

For the mutated gene itself no curated direction knowledge is needed. The rule
is mechanical and gene-agnostic: **a drug that activates a deficient product is
compensatory; one that inhibits it is harmful.** That is the portable half of
the direction filter, and this is where it is tested.

| Control | Assoc. | Approved | Pass gate | **Survive** | |
|---|---:|---:|---:|---:|---|
| **CFTR** — cystic fibrosis | 49 | 19 | 7 | **5** | tezacaftor, elexacaftor, deutivacaftor, vanzacaftor, ivacaftor |
| **PAH** — phenylketonuria | 7 | 5 | 2 | **2** | sapropterin (both salt forms) |
| **GBA1** — Gaucher type 1 | 35 | 8 | 1 | 0 | imiglucerase absent from DGIdb |
| *AURKB* — **negative control** | 91 | — | **0** | **0** | correctly rejected |

Ivacaftor is admitted on its own annotation — `activator, positive modulator` /
`ACTIVATING` — with five distinct sources, and the pipeline reaches
COMPENSATORY without being told what cystic fibrosis is. AURKB, whose 91
associations include the best-evidenced compounds in the entire Track 2 result,
returns zero.

**Two of three positive controls pass and the negative control rejects. The
BUB1B zero is therefore a statement about BUB1B.**

The third is informative rather than embarrassing. Imiglucerase is an enzyme
replacement therapy, and DGIdb — built from drug–*target* interaction sources —
largely does not carry biologics. That is a systematic gap in the database, not
a failure of the gate, and it is a real limitation of any repurposing search
built on these resources.

**What this control does not cover.** It validates the gate on the *mutated
gene*. BUB1B posed the harder problem: BubR1 has no drug at all, so the search
had to go outward to network neighbours, and nothing here proves the gate would
recognise a good drug two nodes away. The control bounds the claim; it does not
remove it.

### 3.6 What this negative does and does not establish

It establishes that **DGIdb contains no well-evidenced, directionally useful,
approved drug for this network**, and that the well-evidenced pharmacology it
does contain is contraindicated by direction of effect. It does not establish
that no drug could help. Specifically not run:

- **only approved drugs** were carried forward; the clinical-stage TTK, PLK1,
  AURKB and CDK1 inhibitors of §3.3 were excluded by that filter, and — on the
  direction argument — rightly so;
- **only DGIdb** was queried; ChEMBL and Open Targets `knownDrugs` were not;
- **signature-based repurposing was not run at all** — LINCS / Connectivity Map
  is the field's standard approach when no target is druggable, and its absence
  is the largest gap in this analysis;
- **aneuploidy-selective compounds** (Tang et al., *Cell* 2011 — AICAR, 17-AAG,
  chloroquine) are a live hypothesis for the cancer-risk arm and were not
  pursued.

---

## 4. The variant-class route

A premature termination codon is actionable at the ribosome, which sidesteps
both the undruggability of BubR1 and the direction problem entirely.

### 4.1 The PTC, computed from the MANE CDS

```
transcript ENST00000287598 · CDS 3,153 nt · 1,051 codons
codon 737   TTA  (Leu)
c.2210T>G   TTA → TGA  (Ter)
context     ...CCAGAG [TGA] AGTGCC...
```

| Feature | Value | Rank | Effect size |
|---|---|---|---|
| Stop codon | **UGA** | 1 of 3 | UGA vs UAA ≈ an order of magnitude in reporter assays |
| +4 nucleotide | **A** | 3 of 4 | C vs A ≈ 2–3-fold |

Mixed, and the two are not commensurable — which is why they are reported with
effect sizes rather than as competing ranks.

### 4.2 What readthrough actually produces

At UGA, near-cognate incorporation inserts **tryptophan, cysteine or arginine —
not leucine**. Successful readthrough therefore yields **BubR1
p.Leu737Trp/Cys/Arg**, a *novel missense variant*, in a patient whose other
allele is already an uncharacterised missense.

This has a direct consequence for how the hypothesis must be tested: **a western
blot for full-length BubR1 cannot distinguish restored function from a
full-length non-functional product.** The functional readout is not a
confirmation step. It is the experiment.

### 4.3 Why the mechanism wants a dual activity

The PTC transcript is NMD-degraded, so readthrough has little substrate.
Patients with higher transcript levels respond better to readthrough drugs, and
co-administering an NMD inhibitor with gentamicin restored full-length protein in
a Hurler syndrome model. Readthrough itself partially inhibits NMD, which helps.

This requirement is what makes amlexanox mechanistically attractive and its
unavailability genuinely costly.

### 4.4 The class-level evidence is not encouraging, and that belongs here

**Ataluren is the only purpose-built PTC readthrough drug taken through a full
phase-3 programme, and it failed on efficacy.** Its EMA conditional
authorisation was not renewed — CHMP negative January 2024, re-examined and
confirmed October 2024, adopted by the Commission. Aminoglycoside readthrough
trials in cystic fibrosis and Duchenne produced marginal and inconsistent
results.

That is the single most informative datapoint about whether readthrough restores
clinically meaningful protein, and it argues **against** this strategy. Any
honest reading of this proposal has to weigh it. Ours is that the approach
remains worth testing *ex vivo* — where it is cheap and fast — precisely because
the clinical record says it should not be assumed.

---

## 5. The patient screen

Candidate liabilities are matched automatically against the proband's own HPO
terms, read from the phenotype file the pipeline extracts:

| HPO | Feature | Contraindicated liability |
|---|---|---|
| `HP:0000121` | Nephrocalcinosis | **nephrotoxic** — reduced renal reserve |
| `HP:0002859` | Rhabdomyosarcoma | **immunosuppressive**, **proliferative risk** |
| `HP:0001508` | Failure to thrive | GI intolerance |

It changed the answer twice. Gentamicin has the best readthrough evidence and the
closest precedent to BubR1 — aminoglycosides functionally restored *BRCA1*
nonsense alleles, a large nuclear genome-maintenance protein — and it is
contraindicated in a child with calcium deposits in the kidney since birth.
Amlexanox is the best mechanistic fit and its TBK1/IKKε inhibition suppresses
innate immune and type-I interferon signalling in a child with an active cancer
predisposition syndrome.

**A ranking by mechanism alone would have proposed gentamicin, then amlexanox.
Both are unusable in this patient, and it takes a phenotype file to see it.**

---

## 6. What would falsify this

In order, each step gating the next:

1. **Quantify *BUB1B* transcript** in patient-derived fibroblasts. Readthrough
   response tracks transcript abundance and NMD efficiency varies by allele. If
   the PTC transcript is absent, readthrough has no substrate and the hypothesis
   dies here.
2. **Measure baseline BubR1 protein**, and — critically — **measure allele 2's
   residual function**. The entire dose–response argument is unquantified without
   it. If p.Asn1002Lys is functionally null, restoring allele 1 lands at best at
   the carrier level, and that is the ceiling.
3. **Test readthrough ex vivo**: escin, with gentamicin as the mechanistic
   benchmark. Primary readout full-length BubR1 — but see §4.2: full-length is
   not the same as functional.
4. **Functional readout**: premature chromatid separation rate and aneuploidy in
   treated versus untreated patient fibroblasts. **This is the decisive
   experiment.** Protein restoration without functional rescue is a negative
   result and should be reported as one.
5. **Resolve phase.** The compound-heterozygous configuration is presumed, not
   proven (Track 1 §8). Parental testing or long-read sequencing settles it, and
   it is cheap. If the variants are in *cis*, the mechanism and this entire
   proposal collapse — that makes it the first thing to check, not a caveat.

### What clinical benefit could look like — and what it cannot

The IUGR, the growth restriction and the developmental phenotype were
established *in utero* and are not reversible by any of this. Existing aneuploid
cells stay aneuploid. **The only plausible benefit is a reduction in the rate of
future missegregation, and therefore in future cancer risk** — in one child, over
decades, with no short-term measurable endpoint.

That is a modest and slow claim, and it is the honest one. What demonstrably
helps MVA families today is tumour surveillance; nothing here displaces it, and
any intervention would sit alongside it. A further clinical consideration
deserves stating: rhabdomyosarcoma protocols contain **vincristine**, a spindle
poison, in a patient whose spindle checkpoint is already deficient. That
interaction is outside the scope of this proposal and inside the scope of this
child's care.

---

## 7. Scalability

The generalisable pattern:

```
loss-of-function variant in an undruggable target
   → build the mechanistic neighbourhood (physical interactions + curated core)
   → ADD the variant-class machinery: for a nonsense allele, the ribosome
   → query drug–gene evidence across both
   → filter by DIRECTION OF EFFECT relative to the lesion
   → screen survivors against the patient's own phenotype
```

Two components are worth reusing. The **direction filter** drops inhibitors that
would worsen a loss-of-function disease — a silent failure mode of naive
repurposing, and here not a hypothetical one: it removes **all 43** of the
best-evidenced associations in the result (§3.3), every one of them an inhibitor
of a checkpoint gene in a checkpoint-deficiency disease. The **phenotype screen** changed
the answer twice here and needs nothing but the HPO file every rare disease case
already has.

**One half of the direction filter is already gene-agnostic.** For the mutated
gene itself the rule needs no curation — *activating a deficient product is
compensatory, inhibiting it is harmful* — and §3.5 runs it unchanged against
CFTR, PAH and GBA1. That part ports to any loss-of-function gene as written.

**An honest scope statement for the rest**: the code is not gene-agnostic.
`t2_01`'s curated core, `t2_03`'s network direction sets and `t2_04`'s candidate
list are specific to this pathway and this variant class. What generalises is the *pattern*; porting it to
another gene means rewriting those three knowledge bases. Claiming otherwise —
as an earlier draft did — is contradicted by the files a reviewer can open.

---

## 8. Limitations

1. **No readthrough agent has been tested on *BUB1B*** in any system, and
   escin has no readthrough data in any nuclear or cell-cycle gene.
2. **Expected effect size is small.** Single-digit readthrough against a deficit
   that is an entire gene copy.
3. **The readthrough product is not wild-type protein** (§4.2), and its function
   is unknown.
4. **Allele 2's residual function is unmeasured**, and the dose–response argument
   is unquantified without it.
5. **Escin's registration is herbal/traditional**, not a full marketing
   authorisation, and no exposure calculation was performed — the concentrations
   producing readthrough in the CFTR screen versus achievable human plasma levels
   are the gating pharmacological question and remain open.
6. **The class's clinical record is poor** (§4.4).
7. **Phase was never proven**; if *cis*, this collapses.
8. **The search was partial**: DGIdb only, approved drugs only, no
   signature-based repurposing, no aneuploidy-selective compounds (§3.6). The
   positive control (§3.5) shows the gate works on the mutated gene; it does not
   show it would work two nodes out, which is what BUB1B required.
9. **DGIdb does not carry biologics.** Enzyme replacement therapy for Gaucher
   is absent from it entirely — a systematic blind spot of any repurposing
   search built on drug–target interaction databases.
10. **The mosaicism screen is inconclusive**, not negative (§2.3).

---

## 9. Methods description form

**Team name** — `bralewild` (individual).

**Approach (variant/mechanism → candidate)** — §2–5.

**Automated or manual?** — Both, and the split is worth stating precisely.
*Automated*: the 52-target network; all 424 drug–gene associations; the evidence,
direction and safety-class filters; the PTC sequence context from the Ensembl
MANE CDS; the phenotype screen that reads the HPO file and reorders candidates;
the positive and negative controls of §3.5; the zero-candidate result. *Manual*: the literature identifying readthrough as a
variant-class strategy; the curated candidate set; the direction classifications
in `t2_03`. **No database query returns escin or amlexanox as readthrough agents
— their activity is a primary-literature property, not a drug–target
annotation.** The corrected pipeline searched the right biology, recovered
ataluren and ELX-02 as a positive control, and still did not surface them. That
is now a demonstrated limitation of database-driven repurposing rather than an
assertion about it.

**Describe the manual review** — targeted searches on BubR1 mechanism and the
MCC; MVA management and prior therapeutic attempts; readthrough agents and their
regulatory status; NMD and combined readthrough/NMD strategies; stop-codon and
sequence-context effects; screens of approved compounds for PTC suppression;
aneuploidy-selective compounds. Every claim carried into this report is cited.
Claims that failed verification were removed — including our own statements that
ataluren was EMA-approved, that amlexanox was marketed, and that BubR1's kinase
activity was the relevant function. An adversarial review of the first draft is
what surfaced them; its findings are recorded in the repository history.

**Public or proprietary?** — Public only. STRING v12, DGIdb v5, Open Targets,
Ensembl REST, OMIM, PubMed/Europe PMC, EMA and FDA public records, the Human
Phenotype Ontology. Only gene symbols and genomic coordinates were transmitted —
no subject identifier.

**Proprietary sources** — none.

**How was the mechanism characterised?** — §2. Loss of function, from the
consequence type of both alleles, BubR1's stoichiometric role in the MCC, the
published dose–response, and phenotype coherence. Explicitly *not* from kinase
activity, which BubR1 is not established to have.

**Time and effort** — the pipeline runs in **under two minutes at zero cost**
(one STRING call, three batched DGIdb calls, one Open Targets, one Ensembl). The
substantive effort was mechanistic reading, verification, and a four-reviewer
adversarial audit of the first draft that invalidated several of its central
claims.

---

## 10. Method abstract

Mosaic Variegated Aneuploidy has no established treatment. Starting from the
Track 1 result — biallelic *BUB1B*, p.Leu737Ter and p.Asn1002Lys — this work asks
whether any available medicine could plausibly act on it.

The lesion is unambiguous loss of function. BubR1 stoichiometrically inhibits
CDC20, holding the anaphase-promoting complex inactive until chromosomes are
correctly attached; losing it produces premature chromatid separation,
missegregation, mosaic aneuploidy, chromosomal instability and the cancer
predisposition that presented here as rhabdomyosarcoma. The relationship is
graded: premature chromatid separation severity tracks residual expression
continuously, with detectable aneuploidy below roughly 50 %.

We asked the obvious question computationally. A pipeline built the mechanistic
neighbourhood — 52 targets across the checkpoint complex, the APC/C, and, because
one allele is a premature termination codon, the translation-termination and NMD
machinery — and retrieved 424 drug–gene associations, 41 involving an approved
drug. None survived filters on evidence quality, direction of effect and safety
class.

The shape of that negative is the informative part. Thirty-eight of the 52
targets have no reported drug at all. Every one of the 41 approved associations
rests on a single source. And the 43 associations that *are* well evidenced
collapse onto four genes — AURKB, CDK1, PLK1, TTK — with 41 of the 43 explicitly
typed as inhibitors and **not one acting in the compensatory direction**. That is
not a coincidence: these are oncology clinical-stage compounds, and oncology
develops checkpoint inhibitors precisely to force missegregation and kill
dividing cells, which is this disease deliberately induced. The pharmacology for
this pathway is well developed and aimed in exactly the direction that would harm
this patient. Including the ribosome in the network matters too: the search
independently recovers ataluren and ELX-02, confirming it can see readthrough
biology, and neither is marketed.

The proposal is therefore **escin**, a marketed triterpene saponin identified as
a readthrough inducer in an unbiased screen of ~1,600 approved compounds against
CFTR premature termination codons. It is the only candidate that is
simultaneously evidenced, available, and free of a liability this child's
phenotype contraindicates. An automated screen against the proband's HPO terms
demoted the two better-known alternatives: gentamicin, whose nephrotoxicity
conflicts with nephrocalcinosis, and amlexanox, whose TBK1/IKKε inhibition
conflicts with an active cancer predisposition — and which has, in any case, no
marketed product since 2019.

**Strengths.** The negative is computed, auditable and **controlled**: the same
gate admits the CFTR modulators and sapropterin from their own annotations while
rejecting all 91 AURKB associations, so the BUB1B zero is a statement about
BUB1B rather than about the filters. It is a shaped negative rather than an
empty one. The direction
filter and the phenotype screen address failure modes naive repurposing does not,
and both changed the answer here rather than sitting unexercised. The pipeline runs in under two minutes at zero cost,
and adding the variant-class machinery gives it a working positive control.

**Limitations, stated plainly.** No readthrough agent has been tested on *BUB1B*.
The expected effect is small: single-digit readthrough against a deficit of an
entire gene copy, and with one null allele the ceiling is the healthy-carrier
level. Readthrough at UGA inserts tryptophan, cysteine or arginine — the product
is a novel missense protein, not restored BubR1, so full-length protein on a blot
is not evidence of rescue. Allele 2's residual function is unmeasured. No
exposure calculation was performed. And ataluren, the only purpose-built agent
taken through phase 3, failed on efficacy — class-level evidence that argues
against this strategy and that we weigh rather than omit.

These are hypotheses for laboratory follow-up. The decisive experiment — premature
chromatid separation rate in treated patient fibroblasts — is cheap, fast, and
would settle it.

---

## 11. Acknowledgement

> This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON
> (The Benchmarking, Evaluation, and Assessment Consortium for Science), with
> prize sponsorship from AWS and Anthropic. We are deeply grateful to the child
> and their family who generously contributed their data and their story to
> advance research into this rare disease. We acknowledge their trust in making
> this Hackathon possible.

Protocol approved by **WCG IRB #20252010**. Released under **CC BY 4.0**.

### References

- Hanks S. *et al.* Constitutional aneuploidy and cancer predisposition caused by
  biallelic mutations in *BUB1B*. *Nat Genet* 2004.
- **Gradual reduction of BUBR1 protein levels results in premature
  sister-chromatid separation then in aneuploidy.** *Hum Genet* 2008,
  PMID 18932004 — the source of the dose–response and the ~50 % figure.
- Mutyam V. *et al.* Discovery of clinically approved agents that promote
  suppression of CFTR nonsense mutations. *Am J Respir Crit Care Med* 2016 —
  escin.
- Atkinson J. *et al.* Amlexanox enhances premature termination codon
  read-through in *COL7A1*. *J Invest Dermatol* 2017.
- Amlexanox: readthrough induction and NMD inhibition in a Charcot–Marie–Tooth
  hiPSC model harbouring a *GDAP1* nonsense mutation. *PMC10385573*, 2023.
- Functional restoration of *BRCA1* nonsense mutations by aminoglycoside-induced
  readthrough. *PMC9273842*.
- Translational readthrough of *BBS2* and *ALMS1* restores protein, ciliogenesis
  and function in patient fibroblasts. *PMC8353411*.
- Howard M. *et al.* Sequence specificity of aminoglycoside-induced stop codon
  readthrough. *Ann Neurol* 2000.
- Loughran G. *et al.* Evidence of efficient stop codon readthrough in four
  mammalian genes. *Nucleic Acids Res* 2014.
- Tang Y.-C. *et al.* Identification of aneuploidy-selective antiproliferation
  compounds. *Cell* 2011 — cited, and not pursued (§3.6).
- EMA / CHMP public assessment record, Translarna (ataluren), 2024.
- German Commission E monograph, *Aesculus hippocastanum* seed extract, 1984 /
  1994.
