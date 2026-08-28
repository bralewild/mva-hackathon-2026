# 3-Minute Pitch Video — script and slides

**Deliverable:** required for Track 2, recommended for Track 1 ([02_compliance.md](02_compliance.md) §7).
**Runtime:** 3:04 at 150 wpm — never let it pass 3:00.
**Script length:** 458 spoken words. Slide timings below are derived from the
per-slide counts, not estimated.

The panel scores **rigour 35 % · impact 25 % · innovation 25 % · scalability 15 %**.
Every beat below is mapped to one of those. Nothing is in this script that does
not earn a rubric point.

---

## Recording notes

**Use your own voice.** A synthetic presenter costs credibility on a rubric that
weights rigour highest, and a panel that reads a report about verifying claims
will notice. An accent is not a defect; unclear audio is. Record in a small room
with soft furnishings, phone or headset microphone 15–20 cm off-axis from your
mouth so plosives miss it.

**Record audio and slides separately.** Read the script in one pass, re-record
individual paragraphs if needed, then time the slides to the audio. Trying to do
both at once is what makes people rush the ending.

**Speak slower than feels natural.** The word budget already assumes 150 wpm. If
you finish at 2:30, that is better than finishing at 3:10 and being cut off.

**Show the terminal at least once.** Real output on screen while you talk is
worth more than any diagram.

---

## Slide 1 — Title · 0:00–0:09

> **Rare Disease, Real Kid — MVA Hackathon 2026**
> A blind genome-wide search, and what it found
> `bralewild` · github.com/bralewild/mva-hackathon-2026

**Say:**

> One child. Five million variants. I wanted to know whether a pipeline could
> reach the diagnosis without being told what to look for.

*(23 words · 9 s)*

---

## Slide 2 — The design decision · 0:09–0:38

> **The temptation:** we know it's MVA → look at BUB1B, CEP57, TRIP13 → find it
> **The problem:** that method helps no future patient
>
> **So: nothing downstream of stage 01 knew the disease name, the gene list, or
> the inheritance pattern.**

**Say:**

> Here is the temptation. We already know this is Mosaic Variegated Aneuploidy.
> Three genes cause it. Check those three, find the variant, done.
>
> But a method that only works when you already know the answer helps no future
> patient. And reusability is the stated goal of this hackathon.
>
> So I ran it blind. No gene list. No disease name. No inheritance hint. Every
> stage saw only the VCF and the patient's phenotype terms.

*(73 words · 29 s)*

---

## Slide 3 — The funnel · 0:38–1:01

> `5,012,204` → `9,145` → `179 rare` → `140 genes ranked`
> **BUB1B — rank 1 · margin 22.8 % · 5 of 5 convergence criteria**
> `NM_001211.6:c.2210T>G` p.Leu737Ter + `c.3006T>G` p.Asn1002Lys
>
> **Scored 100.0 / 100 · F-max 1.000 · full match at rank 1**

**Say:**

> Four filters, and a ranking by semantic similarity to the child's phenotype.
>
> BUB1B came out first, by a twenty-three percent margin, with a compound
> heterozygous pair — a stop codon and a missense. The disease-name check ran
> only *after* the ranking was closed, so it could not have steered it.
>
> It scored one hundred out of one hundred.

*(57 words · 23 s)*

---

## Slide 4 — The pipeline caught itself · 1:01–1:31

> **First pass:** eight autosomes "mosaic" at 5–7 % — a plausible, exciting finding
> **The discriminator:** a real whole-chromosome event is *uniform*; an artefact *concentrates*
> **Result:** no chromosome reached 80 % uniformity → **finding withdrawn,
> reported as inconclusive rather than as a discovery**
>
> Also caught: 201 SERPINA1 "variants" at 13–20 % VAF — pseudogene mismapping

**Say:**

> This is the part worth seeing.
>
> The first mosaicism pass produced a beautiful finding: eight chromosomes,
> five to seven percent mosaic. It was wrong. A real whole-chromosome event
> affects every window uniformly; a mapping artefact concentrates in a few. Not
> one chromosome passed that test, so it went into the report as inconclusive —
> not as a result.
>
> A pipeline that detects its own false positives is worth more than one that
> never fails in a demonstration.

*(76 words · 30 s)*

---

## Slide 5 — Track 2: the search · 1:31–1:53

> **52 targets** — checkpoint complex, APC/C, and the ribosome
> **424** drug–gene associations · **41 approved** — *all single-source*
> **43 well-evidenced** — *none approved* — AURKB, CDK1, PLK1, TTK
>
> ### Acting in the compensatory direction: **0**

**Say:**

> Track 2 asks: is there anything we could actually give this child?
>
> I pulled every reported drug for the pathway. Nothing survived.
>
> But look at the *shape* of that negative. The forty-three best-evidenced drugs
> collapse onto four checkpoint genes, and forty-one of them are inhibitors. Not
> one acts in the direction that would help.

*(54 words · 22 s)*

---

## Slide 6 — Three searches, all pointing the wrong way · 1:53–2:20

> Oncology develops checkpoint inhibitors **precisely to force missegregation**
> and kill dividing cells. **That is this disease, deliberately induced.**
>
> Connectivity Map — needs no drug target — returned **vincristine**.
> A spindle poison already in his chemotherapy. For his spindle defect.
>
> **Three independent searches. All aimed the wrong way.**

**Say:**

> And that is not an accident. Oncology develops checkpoint inhibitors precisely
> to force missegregation and kill dividing cells — this child's disease, induced
> on purpose.
>
> So I ran the standard alternative too. Connectivity Map needs no drug target at
> all. It returned vincristine — a spindle poison already in his chemotherapy —
> as a therapy for his spindle defect.
>
> Three independent searches. Every one of them aimed the wrong way.

*(67 words · 27 s)*

---

## Slide 7 — So the proposal comes from the variant class · 2:20–2:42

> One allele is a **premature stop codon** — addressable at the ribosome, not the gene
>
> Safety screen vs. the child's own HPO terms:
> | | |
> |---|---|
> | gentamicin | best readthrough evidence — **nephrotoxic** vs nephrocalcinosis |
> | amlexanox | best mechanistic fit — **immunosuppressive** vs rhabdomyosarcoma |
> | **escin** | marketed · no phenotype conflict · **the only proposable candidate** |

**Say:**

> So the proposal comes from the variant class instead. One allele is a
> premature stop codon, and that is addressable at the ribosome.
>
> Then I screened the candidates against the child's own phenotype terms — and it
> changed the answer twice. The two best-known agents are both contraindicated
> by his own record. One candidate survives.

*(54 words · 22 s)*

---

## Slide 8 — What generalises · 2:42–3:04

> **Direction filter** — drops drugs that would worsen a loss-of-function disease
> **Phenotype screen** — needs only the HPO file every rare-disease case already has
>
> Expected effect size is **small**. The decisive experiment is **cheap**.
>
> `github.com/bralewild/mva-hackathon-2026`

**Say:**

> Two pieces generalise: a filter for direction of effect, and a screen against
> the patient's own phenotype.
>
> I am not claiming a cure. The expected effect is small, and I say so in the
> report's first paragraph. But the experiment that would settle it is cheap —
> and this child deserves someone to run it.

*(54 words · 22 s)*

---

## Totals

| Slide | Words | Seconds | Rubric |
|---|---:|---:|---|
| 1 Title | 23 | 9 | — |
| 2 Blind design | 73 | 29 | innovation, scalability |
| 3 Result | 57 | 23 | rigour, impact |
| 4 Self-caught error | 76 | 30 | **rigour** |
| 5 Track 2 search | 54 | 22 | rigour, innovation |
| 6 Three searches converge | 67 | 27 | **innovation** |
| 7 Phenotype screen | 54 | 22 | innovation, impact |
| 8 Generalises + close | 54 | 22 | **scalability**, impact |
| **Total** | **458** | **184** | |

**458 words → 3:04 at 150 wpm.** These counts and every slide timing above are derived from the script blocks, not typed by hand. If you edit a line, re-derive them.

If you run long on the day, cut slide 3's *"Four filters"* line. **Never cut slides 4, 6 or 7** — they carry rigour and innovation, which are 60 % of the score.

---

## Compliance

- No patient identifiers, no raw data, no screenshots of the clinical document
  may appear on screen. The HPO terms named aloud (nephrocalcinosis,
  rhabdomyosarcoma) appear in the challenge's own public case description.
- Upload unlisted to YouTube or Vimeo; a private link is acceptable, a
  login-walled one is not.
- Include the acknowledgement in the video description, not on a slide — it
  costs 15 seconds you do not have:

  > This work was made possible through the Hackathon, organized by Sage
  > Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON,
  > with prize sponsorship from AWS and Anthropic. We are deeply grateful to the
  > child and their family who generously contributed their data and their story.
