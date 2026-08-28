# 3-Minute Pitch Video — narration and slides

**Deliverable:** required for Track 2, recommended for Track 1 ([02_compliance.md](02_compliance.md) §7).
**Slides:** published deck — open full screen, press `N` to hide the timing chrome before you record.
**Runtime:** 2:57 at 150 wpm. Never let it pass 3:00.

The panel scores **rigour 35 % · impact 25 % · innovation 25 % · scalability 15 %**.
Every beat is mapped to one of those. Nothing here does not earn a rubric point.

> **On sounding human.** This narration is written to be *spoken*, not read — contractions,
> short sentences, natural breaks. The earlier draft was report prose, which is why it would
> have sounded robotic no matter who read it. The emotional weight is carried by the honesty
> (slides 4 and 8), not by adjectives: with rigour at 35 %, sentiment costs credibility and
> loses both halves at once.

---

## Recording

**Use your own voice.** A synthetic presenter costs credibility on a rubric that weights
rigour highest. An accent is not a defect; unclear audio is. Small room, soft furnishings,
microphone 15–20 cm off-axis from your mouth so plosives miss it.

**Record audio and slides separately**, then align. Doing both at once is what makes people
rush the ending. In the deck: `←` `→` to move, **`N`** hides the timing chrome, **`S`** shows
these notes on screen, **`F`** full screen.

**Speak slower than feels natural.** Finishing at 2:40 beats finishing at 3:10 and being cut.

**Show the terminal once** if you can — real output while you talk is worth more than any diagram.

---

## Slide 1 — Title — a boy, and a question · 0:00–0:11

**Say:**

> This is about a boy whose cancer arrived before his diagnosis did.
> Five million variants in his genome — and nobody allowed to tell the pipeline where to look.

*(28 words · 11 s · the human stake)*

---

## Slide 2 — The design decision · 0:11–0:33

**Say:**

> Here's the temptation. We already know this is Mosaic Variegated Aneuploidy. Three genes cause it. Check those three, done.
> But a method that only works when you already know the answer helps nobody's child.
> So I ran it blind. No gene list. No disease name. No inheritance hint. Every stage saw the VCF and his symptoms.

*(56 words · 22 s · innovation, scalability)*

---

## Slide 3 — The result · 0:33–0:51

**Say:**

> BUB1B came out first. Twenty-three percent clear of the runner-up, with a compound heterozygous pair — a stop codon and a missense.
> And the disease-name check ran only after the ranking was closed, so it couldn't have steered anything.
> It scored a hundred out of a hundred.

*(46 words · 18 s · rigour, impact)*

---

## Slide 4 — The part worth seeing · 0:51–1:24

**Say:**

> Here's the part I'd want you to see.
> The first mosaicism pass gave me a beautiful finding. Eight chromosomes. It looked real. It wasn't.
> A true whole-chromosome event shows up in every window; an artefact bunches into a few. Not one passed — so it went into the report as inconclusive, not a discovery.
> A pipeline that catches its own false positives is worth more than one that never fails in a demo. And a family may one day read what it says.

*(82 words · 33 s · **rigour**)*

---

## Slide 5 — Track 2 — is there anything we could give him? · 1:24–1:45

**Say:**

> Track two asks the next question. Is there anything we could actually give him?
> I pulled every reported drug for the pathway. Nothing survived.
> But look at the shape of that negative. The forty-three best-evidenced drugs land on four checkpoint genes — and forty-one are inhibitors. Not one pushes the direction that would help.

*(53 words · 21 s · rigour, innovation)*

---

## Slide 6 — Three searches, all aimed the wrong way · 1:45–2:09

**Say:**

> That's no accident. Oncology builds checkpoint inhibitors precisely to force missegregation and kill dividing cells — his disease, induced on purpose.
> So I ran the standard alternative. Connectivity Map needs no drug target at all.
> It handed me vincristine, a spindle poison already in his chemotherapy, as a treatment for his spindle defect.
> Three searches. Every one aimed the wrong way.

*(60 words · 24 s · **innovation**)*

---

## Slide 7 — Screened against his own phenotype · 2:09–2:30

**Say:**

> So the proposal comes from the variant class instead. One allele is a premature stop codon, and that you can act on at the ribosome.
> Then I screened the candidates against his own symptoms. It changed the answer twice. The two best-known agents are ruled out by his own record. One survives.

*(52 words · 21 s · innovation, impact)*

---

## Slide 8 — What generalises · 2:30–2:57

**Say:**

> Two pieces generalise: a filter for direction of effect, and a screen against the patient's own phenotype.
> I'm not claiming a cure. The expected effect is small — and that sentence sits in the report's first paragraph, not buried in the last. A family deserves the ceiling before the hope.
> But the experiment that would settle it costs almost nothing. And this child deserves someone to run it.

*(67 words · 27 s · **scalability**, impact)*

---

## Totals

| Slide | Words | Seconds | Rubric |
|---|---:|---:|---|
| 1 Title | 28 | 11 | the human stake |
| 2 The design decision | 56 | 22 | innovation, scalability |
| 3 The result | 46 | 18 | rigour, impact |
| 4 The part worth seeing | 82 | 33 | **rigour** |
| 5 Track 2 | 53 | 21 | rigour, innovation |
| 6 Three searches, all aimed the wrong way | 60 | 24 | **innovation** |
| 7 Screened against his own phenotype | 52 | 21 | innovation, impact |
| 8 What generalises | 67 | 27 | **scalability**, impact |
| **Total** | **444** | **177** | |

**444 words → 2:57 at 150 wpm.** Every count and timing on this page is generated from the
deck's own narration data, not typed by hand — edit the deck, regenerate this.

**Never cut slides 4, 6 or 8.** Slide 4 is the rigour beat, 6 is the innovation beat, and 8
is where the honesty about effect size lands. Those are 60 % of the score and the whole
human argument.

---

## Compliance

- No patient identifiers, no raw data, no screenshots of the clinical document on screen.
  The HPO terms named aloud (nephrocalcinosis, rhabdomyosarcoma) are in the challenge's own
  public case description; the deck names neither the child nor any family member.
- Upload unlisted to YouTube or Vimeo. A private link is acceptable; a login-walled one is not.
- Put the acknowledgement in the video description, not on a slide — it costs 15 seconds you
  do not have:

  > This work was made possible through the Hackathon, organized by Sage Bionetworks in
  > partnership with the MVA Society, Hugging Face, and BEACON, with prize sponsorship from
  > AWS and Anthropic. We are deeply grateful to the child and their family who generously
  > contributed their data and their story.
