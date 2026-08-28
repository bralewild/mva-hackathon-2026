# Track 2 — Drug Repositioning: Methods Report

**Participant:** `bralewild` (individual)
**Track 1 result this builds on:** biallelic *BUB1B* — `NM_001211.6:c.2210T>G`
(p.Leu737Ter) and `c.3006T>G` (p.Asn1002Lys), scored 100/100 with F-max 1.000
**Code:** see repository README; Track 2 pipeline in `pipeline/track2/`
**Date:** 2026-08-28

---

## 1. The proposal in one paragraph

The lesion is a loss of BubR1 function in an undruggable protein complex. A
systematic search across the entire spindle assembly checkpoint and APC/C —
61 targets, 668 drug–gene associations — returns **zero viable candidates**, and
the reasons are structural rather than accidental. That negative is the argument
for the proposal: do not act on the pathway, act on the ribosome.

**Lead candidate: amlexanox** — an FDA-approved anti-inflammatory reported to
act simultaneously as a translational readthrough inducer *and* a
nonsense-mediated decay inhibitor. One allele carries a premature termination
codon; the transcript is NMD-degraded; readthrough alone would have no
substrate. The dual activity is precisely what this lesion requires.

The therapeutic goal is not cure but a **threshold crossing**: aneuploidy
appears below roughly 50 % residual BUB1B expression, and this patient sits near
that line.

---

## 2. Mechanism characterization

### 2.1 Direction of effect: loss of function

*BUB1B* encodes BubR1, a core component of the mitotic checkpoint complex (MCC).
BubR1 is bifunctional: an enzymatic role requiring CENPE-dependent kinase
activation at kinetochores, and a **stoichiometric role as a direct inhibitor of
CDC20**. The MCC — BubR1, BUB3, MAD2L1 — sequesters CDC20 and prevents it from
activating the anaphase-promoting complex/cyclosome until every chromosome is
correctly attached.

Both of this patient's alleles reduce function, by different routes:

| Allele | Change | Consequence |
|---|---|---|
| 1 | `c.2210T>G` p.Leu737Ter | nonsense in exon 17/23; NMD-predicted; **≈ 0 % protein** |
| 2 | `c.3006T>G` p.Asn1002Lys | missense in exon 23/23, **within the kinase domain**; full-length protein, reduced enzymatic function |

This is unambiguously **loss of function**. There is no gain-of-function or
dominant-negative component to argue for.

### 2.2 Pathway disrupted, and the downstream consequence

```
BubR1 ↓ → MCC cannot restrain CDC20 → APC/C activates prematurely
      → sister chromatids separate before correct attachment
      → chromosome missegregation → MOSAIC ANEUPLOIDY
      → chromosomal instability → cancer predisposition
      → rhabdomyosarcoma (HP:0002859), the presenting event
```

The phenotype is coherent with the mechanism end to end: severe IUGR, growth
restriction, nephrocalcinosis, parental recurrent miscarriage, and a childhood
sarcoma are what a constitutional chromosomal-instability disorder produces.

### 2.3 The quantitative anchor

Published work reports **premature chromatid separation in all cells with
reduced BUB1B expression, with aneuploidy detectable below ~50 % residual
expression**. That converts a qualitative mechanism into a measurable target: a
therapy does not need to normalise BubR1, only to push residual function back
above a threshold.

With one null allele and one hypomorphic full-length allele, this patient sits
close to that line. **Recovering even a fraction of allele 1 may be sufficient.**

### 2.4 We looked for the consequence, and did not find it

Stage 08 of the Track 1 pipeline screened for mosaic aneuploidy directly from
the VCF, using B-allele fractions at 2,241,007 heterozygous SNVs across 285
genomic windows. **No whole-chromosome mosaicism is detectable above ~3 % cell
fraction in this blood sample.**

An initial pass appeared to show 5–7 % mosaic fractions on several chromosomes;
a uniformity test rejected it, because the signal concentrated in a minority of
windows on the acrocentric and heterochromatic chromosomes rather than spreading
uniformly as a true whole-chromosome event must.

This negative matters for Track 2: **blood is the wrong tissue for a functional
readout.** Lymphocytes carrying severe aneuploidy proliferate poorly, and
clinical cytogenetics measures aneuploid fractions in *cultured* lymphocytes,
where they are enriched. Any ex vivo test of a therapeutic hypothesis should use
patient fibroblasts, not bulk blood.

---

## 3. The systematic search — and why it returns nothing

Before proposing anything, we asked the obvious question computationally: **is
there already a drug that acts on this pathway?**

`pipeline/track2/` mirrors the Track 1 discipline — numbered stages, traceable
evidence, every discard counted.

### 3.1 The target network (t2-01)

61 targets, tiered by *mechanistic distance from the lesion* rather than
importance:

| Tier | Members | Rationale |
|---|---|---|
| 0 | BUB1B | the mutated gene |
| 1 | BUB1, BUB3, MAD2L1, CDC20, MAD1L1, KNL1 | the MCC and its direct target |
| 2 | APC/C subunits, TTK, AURKB, PLK1, CDK1, CENPE, CCNB1, NEK2 | what the MCC restrains, and the kinases regulating it |
| 3 | 38 further high-confidence STRING partners | the wider mitotic apparatus |

Source: STRING v12 at confidence ≥ 0.9, plus a curated core from the SAC
literature so the network does not depend on one API being available.

### 3.2 The evidence (t2-02)

DGIdb v5 across all 61 targets:

| | |
|---|---|
| Associations returned | **668** |
| With an approved drug | **125** |
| **Carrying a declared interaction direction** | **15 (12 %)** |
| Backed by more than one source | **0** |
| **Targets with no association at all** | **43 of 61** |

### 3.3 The filter, and the result (t2-03)

Three filters, applied in order: evidence quality, **direction of effect**, and
safety class. **Zero candidates survive.**

The direction filter is the part that matters. BubR1's job is to *inhibit*
CDC20, so:

| Node | Effect of inhibiting it | Verdict |
|---|---|---|
| CDC20, APC/C subunits | restrains anaphase onset — **substitutes for the lost BubR1 function** | compensatory |
| BUB1, BUB3, MAD2L1, KNL1, TTK, AURKB, PLK1, CDK1, CENPE | weakens an already-failing checkpoint | **harmful** |

**The direction that would help is the direction nobody builds drugs for.**
CDC20 and every APC/C subunit have zero reported drug associations. The
inhibitors that *do* exist — AURKB, PLK1, TTK, CDK1 — were built by oncology to
push cancer cells *past* a checkpoint. In a child whose checkpoint is already
failing, that is exactly backwards.

What remained after removing harmful directions was single-source text mining:
`PLK1–erythromycin`, `PLK1–lansoprazole`, `CDK1–sertraline`, `PLK1–oleic acid`,
at interaction scores of 0.01–0.11. Those are not pharmacology.

**This is a general problem, not a quirk of this gene.** Drug–target databases
are populated overwhelmingly with inhibitors, because inhibitors are what the
industry builds. For any loss-of-function disease, that makes naive target-based
repurposing structurally unlikely to succeed — and a search that does not filter
by direction will confidently return drugs that would make the patient worse.

---

## 4. The proposal: bypass the pathway

If the pathway is undruggable and the useful direction unavailable, the
remaining option is to stop acting on the pathway. **Act on the ribosome.**

One allele carries a premature termination codon. Translational readthrough
drugs allow the ribosome to insert an amino acid at a PTC and produce
full-length protein. This sidesteps both problems at once: BubR1 does not need
to be druggable, and no pathway node is pushed in the wrong direction.

### 4.1 The PTC context, computed rather than assumed

From the Ensembl MANE Select CDS (`ENST00000287598`, 3,153 nt, 1,051 codons):

```
codon 737 = TTA  (leucine)
c.2210T>G →  TTA → TGA
context:  ...CCAGAG [TGA] AGTGCC...
```

| Feature | Value | Rank | Meaning |
|---|---|---|---|
| Stop codon | **UGA** | 1 of 3 | the most readthrough-permissive (UGA > UAG >> UAA) |
| +1 nucleotide | **A** | 3 of 4 | less favourable (C > U > A > G) |

**Mixed, and reported as mixed.** The stop codon is the best of the three; the
downstream context is not. Readthrough efficiency at this specific allele is an
empirical question and the first thing any follow-up should measure.

### 4.2 Why the mechanism demands a *dual* activity

The p.Leu737Ter transcript is degraded by nonsense-mediated decay. Readthrough
drugs need transcript to act on: patients with higher transcript levels respond
better, and co-administration of an NMD inhibitor with gentamicin restored
full-length protein in a Hurler syndrome model.

So the requirement is not "a readthrough drug". It is **readthrough plus NMD
inhibition**. That requirement is what selects the candidate.

### 4.3 Candidates

#### 1. Amlexanox — the proposal

| | |
|---|---|
| **Approval** | **FDA-approved** (5 % oral paste, aphthous ulcers); oral anti-allergic in Japan |
| **Mechanism** | **dual**: translational readthrough inducer *and* NMD inhibitor |
| **Why it fits** | the lesion requires both activities; amlexanox provides them in one molecule |
| **Precedent** | *COL7A1* in recessive dystrophic epidermolysis bullosa: **8 of 12 PTC alleles responded, some reaching > 50 % of normal full-length protein**. *GDAP1* in patient-derived hiPSC neurons (Charcot–Marie–Tooth). |
| **Patient fit** | no known nephrotoxicity — decisive here (§5) |
| **Limitations** | the approved formulation is topical/oral-local; systemic exposure for a genetic indication would require reformulation and dose-finding. Never tested on *BUB1B*. Readthrough efficiency is allele-specific. |

The 50 % figure in the COL7A1 work is a striking coincidence with the ~50 %
BUB1B threshold. It is a coincidence, in a different gene and a different assay,
and should be read as an order-of-magnitude plausibility argument — not as a
prediction.

#### 2. Aminoglycosides (gentamicin, amikacin) — mechanistic reference, not a proposal

The best-characterised readthrough agents, approved worldwide for decades. The
closest published precedent to BubR1 is compelling: **aminoglycoside-induced
readthrough functionally restored *BRCA1* nonsense alleles** — a large nuclear
protein in genome maintenance, rescued by an approved drug. Also *BBS2*/*ALMS1*
ciliopathies, where both protein and ciliary function were restored in patient
fibroblasts.

**They are nevertheless contraindicated in this patient** (§5). Retained here as
an ex vivo tool compound and as the mechanistic benchmark against which
amlexanox should be compared, not as a therapy.

#### 3. Ataluren — listed for honesty, not proposed

Purpose-built for PTC readthrough, and the reason the class exists clinically.
**It is not currently market-approved.** The EMA's conditional authorisation for
nonsense-mutation Duchenne was not renewed: CHMP negative in January 2024,
annulled on procedural grounds and re-examined, confirmed negative in October
2024, adopted by the European Commission. Efficacy — not safety — was the
ground. Individual member states may still permit named-patient use.

We name it because a report that omitted the obvious candidate would look like
it had not checked.

---

## 5. Patient-specific safety — where a generic answer and this one diverge

Candidate drugs were cross-checked against the proband's own HPO phenotype:

| HPO | Feature | Consequence for drug selection |
|---|---|---|
| **HP:0000121** | **Nephrocalcinosis** | **avoid nephrotoxic agents** — aminoglycosides, and anything requiring renal dose adjustment |
| HP:0002859 | Rhabdomyosarcoma | avoid agents that could promote proliferation or confound oncological surveillance |
| HP:0001508 | Failure to thrive | avoid significant GI intolerance or appetite suppression |

The nephrocalcinosis flag is decisive and it reverses the ranking. The
**best-evidenced** readthrough agents are aminoglycosides — approved worldwide,
with the closest mechanistic precedent — and chronic systemic exposure is
**contraindicated in this specific child**, who has had calcium deposits in the
kidney since birth.

A repurposing proposal that ranked by mechanism alone would have put gentamicin
first and been clinically unusable. Amlexanox leads because it is the candidate
that survives both the mechanism filter and the patient.

---

## 6. What follow-up would look like

Nothing here is a treatment recommendation. These are hypotheses, and they are
falsifiable in a defined order:

1. **Quantify *BUB1B* transcript** in patient-derived cells. Readthrough
   response correlates with transcript abundance, and NMD efficiency varies
   between alleles. If the transcript is absent, readthrough cannot work and the
   hypothesis dies here.
2. **Measure baseline BubR1 protein** against the ~50 % threshold. This
   establishes how far there is to go.
3. **Test readthrough ex vivo in patient fibroblasts** — not blood, for the
   reason established in §2.4. Primary readout: full-length BubR1 on western
   blot. Compare amlexanox against gentamicin as the mechanistic benchmark.
4. **Functional readout**: premature chromatid separation and aneuploidy rate in
   treated versus untreated patient cells. Protein restoration without
   functional rescue would be a negative result, and should be reported as one.

Only if protein and function move together does this become a clinical question
rather than a laboratory one.

---

## 7. Scalability

The pipeline is not *BUB1B*-specific. Nothing in stages t2-01 to t2-04 hardcodes
this gene beyond a single entry point.

The generalisable pattern is:

```
loss-of-function variant in an undruggable target
   → build the mechanistic neighbourhood (STRING + curated core)
   → query drug–gene evidence across it
   → FILTER BY DIRECTION OF EFFECT relative to the lesion
   → if nothing survives: does the variant class itself offer a route?
        nonsense  → readthrough ± NMD modulation
        missense  → chaperone/stabiliser strategies
        splice    → antisense modulation
   → screen surviving candidates against the patient's own phenotype
```

Two components of that are, as far as we can tell, not standard practice:

**The direction filter.** Repurposing searches routinely return inhibitors for
loss-of-function diseases. Ours drops them explicitly and names them, which
turns a silent failure mode into a reported one.

**The phenotype-aware safety screen.** Cross-checking candidate toxicity against
the patient's own HPO terms is what moved gentamicin from first place to
contraindicated. It requires nothing but the phenotype file that every rare
disease case already has.

Roughly 10 % of rare disease alleles are nonsense variants. For every one of
them in an undruggable gene, this pipeline runs unchanged.

---

## 8. Limitations

Stated plainly, because a proposal that hides them cannot be evaluated:

1. **No readthrough agent has ever been tested on *BUB1B*,** in any system. The
   novelty and the risk are the same fact.
2. **Amlexanox's approved formulation is local, not systemic.** Reformulation
   and dose-finding would be required, and systemic pharmacokinetics for this
   indication are unknown.
3. **The +1 sequence context is unfavourable.** UGA is permissive, but the
   adenine immediately downstream is third of four. Efficiency at this allele is
   unmeasured.
4. **The second allele is untouched.** Readthrough addresses p.Leu737Ter only.
   p.Asn1002Lys remains a kinase-domain substitution of unknown residual
   function; if it is functionally null, restoring allele 1 may not be enough.
5. **The ~50 % threshold comes from other patients and cell systems,** not this
   one. It anchors the reasoning but has not been measured here.
6. **Phase was never proven** (Track 1, §8). The compound-heterozygous
   configuration is presumed, not demonstrated, in a singleton.
7. **The COL7A1 > 50 % result is a different gene and a different assay.** It
   supports plausibility, not prediction.
8. **We measured no mosaic aneuploidy in blood** (§2.4), which means the most
   direct functional readout is unavailable in the most accessible tissue.

---

## 9. Methods description form

Answers to `methods_description_form.xlsx`, Track 2 sheet.

**Team name** — `bralewild` (individual).

**Describe your approach in detail (variant/mechanism → candidate medication)** —
see §2–4. In brief: characterise the lesion as loss of function with a
quantitative expression threshold; build the mechanistic neighbourhood; query
drug–gene evidence across it; filter by evidence quality, **direction of effect**
and safety class; observe that nothing survives and why; then ask whether the
*variant class* rather than the *target* offers a route. It does: a premature
termination codon is actionable at the ribosome. Finally, screen the surviving
candidates against the patient's own phenotype.

**Was candidate identification automated, or manual literature review?** —
**Both, in a deliberate division of labour, and it is worth being precise about
which did what.**

*Automated* (`pipeline/track2/`, reproducible, re-runnable): construction of the
61-target network from STRING; retrieval of all 668 drug–gene associations from
DGIdb; the evidence, direction and safety-class filters; the computation of the
PTC sequence context from the Ensembl MANE CDS. **The negative result — zero
surviving candidates — is entirely automated output.**

*Manual* (documented in §9 below): the literature search that identified
readthrough as a variant-class strategy, and amlexanox specifically. No database
query would have returned amlexanox, because its readthrough activity is not a
drug–target annotation — it is a property reported in the primary literature.
The curation of which network nodes are compensatory versus harmful when
inhibited is also a human judgement, encoded explicitly in
`t2_03_mechanism_filter.py` so it can be inspected and disagreed with.

**Describe the manual literature review** — targeted searches on: BubR1
mechanism and the MCC; MVA clinical management and prior therapeutic attempts;
translational readthrough agents and their approval status; nonsense-mediated
decay and combined readthrough/NMD strategies; stop-codon and sequence-context
effects on readthrough efficiency; aneuploidy-selective compounds. Sources are
listed in §11. Each claim carried into the report is cited; claims that could not
be verified were removed — including our own initial statement that ataluren was
EMA-approved, which the EMA record contradicts (§4.3).

**Public or proprietary data sources?** — **Publicly available sources only.**
No proprietary database, no commercial licence, no restricted resource.

**Describe the public data sources** — STRING v12 (protein interaction network);
DGIdb v5 (drug–gene interactions, aggregating DrugBank, ChEMBL, TTD, PharmGKB,
Guide to Pharmacology and text-mined sets); Open Targets Platform (target
metadata and tractability); Ensembl REST (MANE Select transcript and CDS
sequence); OMIM (BUB1B, MVA); PubMed / Europe PMC and journal-hosted full text
for the mechanistic and readthrough literature; EMA and CHMP public records for
regulatory status; the Human Phenotype Ontology for the patient's phenotype
terms. Only genomic coordinates and gene symbols were transmitted to any API —
no subject identifier of any kind.

**Proprietary data sources** — none.

**How did you characterize the variant's mechanism?** — See §2. Loss of
function, established from: the consequence type of both alleles (nonsense with
predicted NMD; missense within the annotated kinase domain), the known
bifunctional role of BubR1 in the MCC, the published dosage threshold below
which aneuploidy appears, and the coherence of the patient's full phenotype with
a constitutional chromosomal-instability disorder. Pathway: spindle assembly
checkpoint → CDC20 → APC/C. Downstream consequence: premature chromatid
separation → mosaic aneuploidy → chromosomal instability → cancer predisposition.

**Estimate of time or effort** — the computational pipeline runs in **under
two minutes** end to end (61 STRING and DGIdb queries plus one Ensembl call) at
**zero cost**, on a laptop. The mechanistic and literature work behind it was
the substantive effort: roughly a day of focused analysis, most of it spent
verifying claims rather than generating them.

**Method abstract** — §10.

---

## 10. Method abstract (500 words)

Mosaic Variegated Aneuploidy has no established treatment. This proposal starts
from the Track 1 result — biallelic *BUB1B*, `p.Leu737Ter` and `p.Asn1002Lys` —
and asks whether any approved medicine could plausibly act on it.

The mechanism is unambiguous loss of function. BubR1 inhibits CDC20, holding the
anaphase-promoting complex inactive until every chromosome is correctly
attached. Losing it lets anaphase proceed early, producing missegregation,
mosaic aneuploidy, chromosomal instability, and the cancer predisposition that
presented here as rhabdomyosarcoma. Critically the relationship is
dosage-dependent with a published threshold: aneuploidy appears below roughly
50 % residual BUB1B expression. That converts an abstract mechanism into a
measurable target.

We first asked the obvious question computationally. A pipeline built the
mechanistic neighbourhood of the lesion — 61 targets across the checkpoint
complex, the APC/C and the mitotic kinases — and retrieved every drug–gene
association reported for them: 668 associations, 125 involving an approved drug.
Three filters followed: evidence quality, direction of effect relative to the
lesion, and safety class. **Zero candidates survived.**

The failure is structural. BubR1 has no drug, and 43 of 61 network members have
no reported association at all. The direction that would help — inhibiting CDC20
or the APC/C, substituting for the lost restraint on anaphase — is precisely the
direction nobody has built drugs for. The inhibitors that exist target BUB1,
TTK, AURKB, PLK1 and CDK1, because oncology wants to push cancer cells *past* a
checkpoint; in a child whose checkpoint is already failing, that is backwards.
Only 12 % of the approved-drug associations record a direction at all, and every
one rests on a single source. A repurposing search that does not filter by
direction will confidently return drugs that would worsen the disease.

That negative motivates the proposal: stop acting on the pathway; act on the
ribosome. One allele carries a premature termination codon, and translational
readthrough restores full-length protein without the target needing to be
druggable. The PTC is UGA, the most permissive stop codon, though the +1 context
is unfavourable; we report both. Because the transcript is NMD-degraded,
readthrough alone would lack substrate: the requirement is readthrough *plus*
NMD inhibition.

That requirement selects **amlexanox**, an FDA-approved anti-inflammatory
reported to do both. It restored >50 % of full-length protein in some *COL7A1*
nonsense alleles and showed activity in patient-derived *GDAP1* neurons, and it
carries no known nephrotoxicity — decisive here, because this child's
nephrocalcinosis contraindicates the aminoglycosides that are otherwise the
best-evidenced readthrough agents. Screening against the patient's own phenotype
is what moved gentamicin from first place to contraindicated.

**Strengths:** the negative is automated and auditable; the direction filter is
explicit; safety is patient-specific; the pipeline runs in under two minutes at
zero cost and generalises to any nonsense variant in an undruggable gene.

**Limitations:** no readthrough agent has been tested on *BUB1B*; amlexanox's
approved formulation is local rather than systemic; the +1 context is
unfavourable; the second allele is untouched; and the 50 % threshold comes from
other patients, not this one. These are hypotheses for laboratory follow-up, not
a treatment.

---

## 11. Acknowledgement

> This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON
> (The Benchmarking, Evaluation, and Assessment Consortium for Science), with
> prize sponsorship from AWS and Anthropic. We are deeply grateful to the child
> and their family who generously contributed their data and their story to
> advance research into this rare disease. We acknowledge their trust in making
> this Hackathon possible.

Data shared under a protocol approved by **WCG IRB #20252010**. This submission
is released under **CC BY 4.0**.

### Key references

- Hanks S. *et al.* Constitutional aneuploidy and cancer predisposition caused by
  biallelic mutations in *BUB1B*. *Nat Genet* 2004.
- OMIM 602860 — *BUB1B*, mitotic checkpoint serine/threonine kinase B.
- Physiological relevance of post-translational regulation of the spindle
  assembly checkpoint protein BubR1. *PMC8066494*.
- Atkinson J. *et al.* Amlexanox enhances premature termination codon
  read-through in *COL7A1* and expression of full-length type VII collagen.
  *J Invest Dermatol* 2017.
- Amlexanox: readthrough induction and NMD inhibition in a Charcot–Marie–Tooth
  model of hiPSC-derived neuronal cells harbouring a nonsense mutation in
  *GDAP1*. *PMC10385573*, 2023.
- Functional restoration of *BRCA1* nonsense mutations by aminoglycoside-induced
  readthrough. *PMC9273842*.
- Translational readthrough of ciliopathy genes *BBS2* and *ALMS1* restores
  protein, ciliogenesis and function in patient fibroblasts. *PMC8353411*.
- Nonsense-mediated mRNA decay efficiency varies in choroideremia, providing a
  target to boost small-molecule therapeutics. *Hum Mol Genet* 2019.
- Howard M. *et al.* Sequence specificity of aminoglycoside-induced stop codon
  readthrough. *Ann Neurol* 2000.
- Loughran G. *et al.* Evidence of efficient stop codon readthrough in four
  mammalian genes. *Nucleic Acids Res* 2014.
- Tang Y.-C. *et al.* Identification of aneuploidy-selective antiproliferation
  compounds. *Cell* 2011.
- EMA / CHMP public assessment record, Translarna (ataluren), 2024.
- Miller D.T. *et al.* ACMG SF v3.2 secondary findings list. *Genet Med* 2023.
