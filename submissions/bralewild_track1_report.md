# Track 1 — Variant Prediction: Methods Report

**Participant:** `bralewild` (individual)
**Model / approach:** `blind-wgs-triage`
**Submission file:** `bralewild_blind-wgs-triage.csv`
**Code:** see repository README for the full pipeline
**Date:** 2026-08-27

---

## 1. Summary

Starting from the proband's WGS VCF (**5,012,204 variants**) and eight HPO terms,
a **blind, genome-wide pipeline** — one that at no point encodes the disease
name, a candidate gene list, or the organisers' compound-heterozygous hint —
ranked **BUB1B** first among 140 candidate genes, with a 22.8 % score margin
over the runner-up.

**Predicted causal pair (GRCh38, biallelic *BUB1B*, MVA type 1, MIM 257300):**

| | Variant 1 | Variant 2 |
|---|---|---|
| Position | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
| HGVSc (MANE `NM_001211.6`) | `c.2210T>G` | `c.3006T>G` |
| HGVSp | `p.Leu737Ter` | `p.Asn1002Lys` |
| Consequence | stop_gained | missense |
| Exon | 17 / 23 | 23 / 23 |
| CADD (phred) | **36** | **24.5** |
| gnomAD AF | 7.87 × 10⁻⁵ | **absent** |
| dbSNP | rs759242053 | — |
| ClinVar | **Pathogenic / Likely pathogenic** | not reported |
| Genotype | `0/1` DP 46 GQ 99 AD 21,25 (VAF 0.54) | `0/1` DP 28 GQ 99 AD 15,13 (VAF 0.46) |
| ACMG | PVS1, PM2, PP4 → **Pathogenic** | PM1, PM2, PP3, PP4 → **VUS, leaning Likely Pathogenic** |

The pair reproduces the canonical molecular architecture of *BUB1B*-related
Mosaic Variegated Aneuploidy: **one truncating allele plus one missense allele**
in compound heterozygosity.

---

## 2. Why a blind pipeline

The hackathon is named after the disease, and the public evaluator source code
states the answer key is a *"compound-heterozygous"* pair. Either fact alone
narrows the search to three genes (*BUB1B*, *CEP57*, *TRIP13*) and a single
inheritance model — an afternoon of work.

We deliberately did **not** do that.

A method that only works when you already know the answer demonstrates nothing
and helps no future patient. The stated goal of this hackathon is that the
methods developed here be **reusable for other undiagnosed individuals**. That
goal is only met by a pipeline that starts from the data.

Accordingly:

- Stages 01–05 contain **no gene list, no disease name, no inheritance hint**.
- The compound-heterozygous expectation is used **only** as a post-hoc
  validation gate (stage 06), applied **after** the ranking is closed.
- A separate, hypothesis-driven targeted analysis of the three known MVA genes
  was run **independently**, and is reported in §7 as an orthogonal check.

---

## 3. Data and environment

| | |
|---|---|
| Input | `WGS_EX2312012_HGWCNDSX7.vcf.gz` (301 MB) + index |
| Phenotype | `Challenge_Clinical_Phenotype_1.docx` → 8 HPO terms |
| Reference | GRCh38, `no_alt_analysis_set_plus_hs38d1`, **Ensembl contig naming** (`1`…`22,X,Y`) |
| Caller | GATK 4.2.4.0 `VariantFiltration`, hard filters (QD<2, MQ<40, ReadPosRankSum<−8, FS>60, MQRankSum<−12.5). No VQSR, no population annotation. |
| Design | **Singleton** — no parental samples |
| Compute | Local only: WSL2 / Ubuntu 24.04 on an Intel i9-14900HX (24 c / 32 t), 32 GB RAM. **No cloud cost.** |
| Raw reads | The 8 FASTQ files (~84 GB) were **not downloaded** — see §8. |

Software: bcftools/samtools/htslib 1.24, snpEff 5.4c with `GRCh38.115`,
Ensembl VEP REST, HPO ontology (`hp.obo`, `genes_to_phenotype.txt`),
Python 3.12, OpenJDK 21. Environment managed with micromamba; all versions
pinned in the repository.

---

## 4. Pipeline

```
5,012,204 variants
     │  01  QC baseline — build, caller, contig naming, FORMAT fields
     │
     │  02  snpEff GRCh38.115, genome-wide, no target regions
     ▼
     │  03  quality + genotype coherence + functional impact + inheritance
     ▼  9,145
     │  04  gnomAD AF, ClinVar, CADD, MANE HGVS (Ensembl VEP REST)
     ▼  179 rare
     │  05  Resnik semantic similarity over the HPO ontology
     ▼  140 genes ranked
     │  06  convergence gate (5 criteria, applied post-ranking)
     ▼  BUB1B
```

### Stage 03 — quality, coherence, impact, inheritance

1. `FILTER = PASS`
2. `GQ ≥ 20`, `DP ≥ 10`
3. **Genotype/VAF coherence** (see §6 — added after a QC finding):
   heterozygous `0.25 ≤ VAF ≤ 0.75`; homozygous `VAF ≥ 0.85`
4. snpEff impact `HIGH` or `MODERATE`
5. Inheritance model, grouped by gene:
   - `AR_COMPOUND_HET` — gene with ≥ 2 heterozygous variants
   - `AR_HOMOZYGOUS` — gene with ≥ 1 homozygous-alternate variant

*De novo* models were **not** evaluated: the case is a singleton, so phase
cannot be established from pedigree.

### Stage 04 — population frequency

An ultra-rare disease (< 50 individuals worldwide) cannot be caused by a
common variant: a 1 % allele frequency implies ~80 million carriers. Variants
with gnomAD AF ≥ 0.01 were removed. **Absence from gnomAD is recorded as ACMG
PM2 evidence.**

Annotation was retrieved through Ensembl VEP REST in batches of 200 rather than
by downloading multi-gigabyte local caches: faster for this scale, always
current, and no additional data at rest.

### Stage 05 — phenotype ranking

Resnik semantic similarity over the HPO ontology:

```
IC(t)      = −log( genes annotated to t or its descendants / total genes )
sim(a, b)  = max IC over the common ancestors of a and b
score(gene)= mean over patient terms of  max_b sim(patient_term, b)
```

The **mean** (not the sum) is used so that genes with many annotations are not
rewarded for volume. Information content weighting means matching
*rhabdomyosarcoma* contributes far more than matching *short stature*.

This is the principle behind Phenomizer and the Exomiser prioritiser,
implemented explicitly and auditably rather than delegated to a black box.

---

## 5. Results

### Reduction funnel

| Step | Remaining | Removed |
|---|---:|---:|
| Total variants | 5,012,204 | — |
| `FILTER = PASS` | 4,740,790 | 271,414 |
| GQ / DP thresholds | 4,676,417 | 64,373 |
| **VAF coherence** | 4,569,022 | **107,395** |
| Impact HIGH / MODERATE | 14,697 | 4,554,325 |
| In a recessive-model gene | 9,145 | — |
| gnomAD AF < 0.01 | **179** | 8,966 |
| Genes with HPO annotation | **140 genes** | — |

Every discard is counted in `results/0X_*_summary.txt`. **No filter is silent.**

### Ranking (top 10)

| # | Gene | Score | Model | Rare vars | Best impact | CADD | ClinVar |
|---|---|---:|---|---:|---|---:|---|
| **1** | **BUB1B** | **1.8476** | AR_COMPOUND_HET | **2** | **HIGH** | **36.0** | **Pathogenic/LP** |
| 2 | ZFP57 | 1.4271 | AR_COMPOUND_HET | 4 | MODERATE | 1.3 | VUS |
| 3 | ITGB4 | 1.2678 | AR_COMPOUND_HET | 1 | MODERATE | 17.9 | — |
| 4 | ATM | 1.2492 | AR_COMPOUND_HET | 1 | MODERATE | 25.4 | VUS/benign |
| 5 | PRDM16 | 1.1283 | AR_COMPOUND_HET | 1 | HIGH | 1.6 | benign |
| 6 | ARID1B | 1.1283 | AR_COMPOUND_HET | 1 | HIGH | 1.3 | — |
| 7 | RFWD3 | 1.0706 | AR_COMPOUND_HET | 1 | HIGH | 9.0 | — |
| 8 | BTNL2 | 1.0408 | AR_COMPOUND_HET | 1 | MODERATE | 22.3 | benign |
| 9 | ZNF808 | 0.9892 | AR_COMPOUND_HET | 1 | MODERATE | 0.1 | — |
| 10 | LAMA3 | 0.8832 | AR_COMPOUND_HET | 1 | MODERATE | 20.6 | — |

Note the **Rare vars** column: after frequency filtering, most ranked genes
retain a **single** rare variant and are therefore no longer viable
compound-heterozygous candidates. *BUB1B* retains its pair.

*ZFP57* matched all 8 HPO terms while *BUB1B* matched 7, yet *BUB1B* scores
higher — information-content weighting rewards **specificity, not count**.

### Convergence gate (stage 06)

| Criterion | Result |
|---|---|
| Margin over runner-up ≥ 15 % | **22.8 %** ✔ |
| Recessive model | AR_COMPOUND_HET ✔ |
| ≥ 2 rare variants retained | 2 ✔ |
| ≥ 1 HIGH-impact variant | 1 ✔ |
| Independent evidence | ClinVar Pathogenic + 2 variants CADD ≥ 20 ✔ |

**5 / 5 — convergence confirmed.**

### Phenotype concordance

| Feature (HPO) | *BUB1B*-MVA |
|---|---|
| Rhabdomyosarcoma `HP:0002859` | Characteristic malignancy of BUB1B-MVA |
| Small for gestational age `HP:0001518` (~1 kg) | Severe IUGR is a hallmark |
| Short stature `HP:0004322` | Hallmark |
| Failure to thrive `HP:0001508`, muscle atrophy `HP:0003202` | Reported |
| Nephrocalcinosis `HP:0000121` | Reported |
| Premature birth `HP:0001622` | Consistent |
| Parental recurrent miscarriage `HP:0200067` | Consistent with a chromosomal-instability disorder |

*BUB1B* encodes BubR1, the mitotic checkpoint kinase. Biallelic loss impairs the
spindle assembly checkpoint → chromosome missegregation → **mosaic aneuploidy**
→ chromosomal instability → cancer predisposition. The phenotype is coherent
with the mechanism end to end, supporting **ACMG PP4**.

---

## 6. Quality-control finding: a paralogue mismapping artefact

The secondary-findings stage surfaced **four distinct frameshift variants in
*SERPINA1* within 61 bp** — biologically impossible in a diploid genome.

Inspection of the locus revealed **201 variants across ~14 kb**, all
heterozygous, with allele fractions of **0.13–0.20** instead of the expected
~0.50, spaced every 2–10 bp, mixing SNVs and indels — **all with GQ = 99 and
DP 50–60**.

*SERPINA1* lies in the SERPINA cluster at 14q32.13, adjacent to the highly
homologous pseudogene *SERPINA2*. Reads originating from the paralogue mismap
onto the real gene and are called as low-fraction heterozygotes.

**The lesson, and the fix:** a high `GQ` means the *caller* is confident, not
that the *variant is real*. Our filters checked caller confidence but never
biological coherence. We added a VAF-coherence filter (het 0.25–0.75,
hom ≥ 0.85), which removed **107,395 variants genome-wide**.

**Robustness check:** after this substantially stricter filter, *BUB1B* remained
rank 1 with an **identical score (1.8476)** and 5/5 convergence criteria. The
finding did not depend on tolerating noise.

---

## 7. Orthogonal confirmation

Independently of the blind pipeline, a targeted analysis of the three known MVA
genes was performed: Ensembl REST gene coordinates → `bcftools` region query →
Ensembl VEP annotation. It found 7 PASS heterozygous variants in *BUB1B*, 2 in
*CEP57*, 3 in *TRIP13*. Of these, only two were not common intronic
polymorphisms (gnomAD 1.6 %–65 %) — **the same two variants** reported here.

Two methodologically independent routes converge on the same pair.

Physical phasing (`PGT`/`PID`) also excluded the *TRIP13* pair `5:893128` /
`5:893132`, which GATK phased **in cis** (same `PID`), and which additionally
showed DP 7 and QUAL 36.8.

---

## 8. Limitations

Stated explicitly, because a method that does not declare its limits cannot be
evaluated:

1. **Phase is not proven.** The case is a singleton. The two *BUB1B* variants
   are **10,911 bp apart** — beyond the reach of short-insert Illumina
   read-backed phasing (~350–550 bp fragments) and outside a single GATK
   assembly region (both carry empty `PGT`/`PID`). Re-processing the 84 GB of
   FASTQ would **not** resolve this. The pair is reported as **presumed in
   trans**, consistent with the recessive disease model and with both alleles
   being ultra-rare, but not demonstrated. Parental testing or long-read
   sequencing would be required.
2. **Small variants only.** The VCF contains SNVs and indels. **No CNVs,
   structural variants, or repeat expansions** were assessed.
3. **Coding-biased.** The HIGH/MODERATE impact filter discards deep intronic and
   regulatory variants — a deliberate sensitivity/noise trade-off.
4. **Annotation coverage.** A gene with no HPO annotations scores 0 even if
   causal. This is intrinsic to any phenotype-driven prioritisation.
5. **The second allele is the weaker one.** `p.Asn1002Lys` is absent from gnomAD
   and lies in the kinase domain (CADD 24.5), but is not in ClinVar. ClinVar
   does contain `c.3006T>A`, a *different nucleotide change producing the same
   amino-acid substitution*, classified as VUS. This is supporting context, not
   PS1 (which requires a previously established **pathogenic** variant).

## Secondary findings

**None reportable.** No variant in the **ACMG SF v3.2** gene list (Miller et
al., *Genet Med* 2023) met our criteria. The single candidate from a
self-curated treatable-disease list — *ATM* `p.Ser978Pro` — is classified
benign / likely benign by multiple ClinVar submitters (gnomAD AF 8.4 × 10⁻⁴)
and is therefore **not** reported.

A well-founded negative is a result. Padding a submission with a benign variant
would be the opposite of the rigour being assessed.

---

## 9. Reproducibility

The Official Rules state that submissions may be **rerun by the organisers**, so
reproducibility is treated as a functional requirement:

- `pipeline/run_all.sh` executes stages 01→07 end to end
- every stage is **idempotent** — valid existing output is skipped
- stage 04 keeps an **incremental cache** and resumes after network failure
  (3 of 47 VEP batches returned HTTP 500 on the first run; `04b_seed_cache.py`
  recovered them without repeating 40 minutes of valid queries)
- output is **validated by content, not exit code** — snpEff exits 0 even when
  it fails, so stage 02 checks file size and variant count
- `.gitattributes` enforces `LF` so scripts run on Linux
- every discard is logged in `results/0X_*_summary.txt`

**Anyone re-running this pipeline must supply their own authorised copy of the
dataset.** No patient data is included in the repository by design: stage `00b`
extracts the HPO terms from the original confidential `.docx`, which each
participant must obtain through the gated Hugging Face dataset.

## Data governance

Patient data is held exclusively on the ext4 filesystem inside WSL, never in the
repository. Three independent barriers: physical separation, `.gitignore`
(`*.vcf*`, `*.bam`, `*.cram`, `*.fastq*`, `*.docx`, `patient_hpo.tsv`), and an
audit script (`99_data_inventory.sh`) that inventories every location holding
data and classifies it as patient-derived or public resource.

Only genomic coordinates (e.g. `15 40209701 T G`) were transmitted to public
annotation APIs — no subject identifier of any kind. All data will be deleted
within 30 days of hackathon close, including derived datasets, with notification
to the organisers.

---

## 10. Methods description form

Answers to `methods_description_form.xlsx`, Track 1 sheet, reproduced here so the
panel has them alongside the report.

**Team name** — `bralewild` (individual participation).

**Model number** — 1 of up to 6. `blind-wgs-triage`.

**Describe your model/approach in detail** — see §4 and §5. In brief: an
eight-stage genome-wide triage that reduces 5,012,204 VCF variants to a ranked
gene list, using no disease name, no gene panel and no inheritance hint. Quality
and genotype/VAF coherence → snpEff functional impact → recessive inheritance
grouping → gnomAD/ClinVar/CADD via Ensembl VEP REST → Resnik semantic similarity
against the proband's eight HPO terms → a post-hoc convergence gate that applies
the compound-heterozygous expectation only after the ranking is closed.

**Is the submission file automated output, or manually curated?** —
**Automated.** `pipeline/08_make_submission.py` derives the CSV deterministically
from the pipeline's own results: gene order comes from the stage 05 ranking, the
variants from stage 04, and EPCR from an explicit evidence rule (rank base, plus
increments for a HIGH-impact allele, a ClinVar pathogenic classification, and a
complete compound-heterozygous pair). Re-running the pipeline regenerates a
byte-identical file. No row was hand-edited.

**Describe any downstream manual review** — none was applied to the submission
file. Human judgement entered upstream, in designing the method, and is
documented rather than hidden: the quality thresholds (§4), the choice of
Ensembl-wide over MANE-only annotation (§4), and the VAF-coherence filter added
after inspecting a mismapping artefact (§6). Each is stated with its rationale
and each is a parameter in `pipeline/00_config.sh`, not a per-variant decision.

**Public or proprietary data?** — **Publicly available data only.** No
proprietary source, no commercial database, no restricted resource.

**Public data used** — Ensembl 115 gene models and coordinates (REST API);
snpEff `GRCh38.115`; Ensembl VEP REST for consequence, MANE Select transcript,
HGVS, CADD phred, dbSNP identifiers, gnomAD allele frequencies and ClinVar
classifications; the Human Phenotype Ontology (`hp.obo`) and its
`genes_to_phenotype.txt` annotations; the ACMG SF v3.2 secondary-findings gene
list (Miller et al., *Genet Med* 2023). Only genomic coordinates were sent to
these APIs — no subject identifier of any kind.

**Proprietary data used** — none.

**Compound-heterozygous pairs, or single variants only?** — **Pairs, natively.**
Stage 03 groups variants by gene and classifies an `AR_COMPOUND_HET` model when a
gene carries two or more qualifying heterozygous variants; stage 08 emits both
alleles in the `chrom_2`/`pos_2`/`ref_2`/`alt_2` columns. Single-variant rows are
emitted for homozygous models or when only one rare variant survives filtering.
The pipeline reports such pairs as **presumed in trans**: the case is a
singleton, so phase cannot be proven (§8).

**How were secondary or incidental findings handled?** — Stage 07 cross-checks
every rare candidate against the ACMG SF v3.2 list and against a separate,
explicitly declared list of genes with treatable disease. **No reportable
secondary finding was identified**, so no `secondary` row appears in the
submission. The single ACMG-adjacent candidate, *ATM* `p.Ser978Pro`, is
classified benign / likely benign by multiple ClinVar submitters (gnomAD AF
8.4 × 10⁻⁴) and was excluded. A well-founded negative is a result; padding the
submission with a benign variant would work against the rigour being assessed.

**Run time and cost** — **~1.5 hours wall clock, USD 0.** Entirely local, no
cloud resources. On an Intel i9-14900HX (24 c / 32 t, 32 GB RAM) under WSL2:
snpEff database download ~1 min (770 MB), genome-wide annotation **58 min**
(single long step), stage 03 ~1 min, stage 04 ~40 min on the first run and
seconds thereafter from the incremental cache, stages 05–08 under 30 seconds
combined. Disk footprint ~2 GB beyond the input. The dominant costs are
annotation CPU time and Ensembl VEP REST latency, both of which scale linearly
with cohort size and parallelise trivially across samples.

### Method abstract (500 words)

Four teams had already reached a perfect leaderboard score when this work began,
so identifying the causal variant was not, on its own, a contribution. The
harder and more useful question was whether a method could find it *without
being told what to look for* — because that is the method that helps the next
undiagnosed child.

The hackathon names the disease in its title, and the evaluator's public source
states the answer is a compound-heterozygous pair. Either fact reduces the
problem to three genes and an afternoon. This pipeline deliberately encodes
neither. Stages 01–05 contain no disease name, no gene list and no inheritance
hint; the compound-heterozygous expectation is applied only in stage 06, as a
validation gate, after the ranking is final.

From 5,012,204 variants: quality filtering and a genotype/VAF coherence check
(added after diagnosing a paralogue-mismapping artefact) leave 4.57 M; snpEff
HIGH/MODERATE impact leaves 14,697; grouping by gene under recessive inheritance
models leaves 9,145 in 4,821 genes; a gnomAD frequency ceiling of 0.01 leaves
179; and Resnik semantic similarity between the proband's eight HPO terms and
`genes_to_phenotype` annotations ranks the 140 genes carrying HPO annotations.

*BUB1B* ranked first, 22.8 % clear of the runner-up, and satisfied all five
convergence criteria. Its two rare variants — `NM_001211.6:c.2210T>G`
p.Leu737Ter (nonsense, exon 17/23, ClinVar Pathogenic/Likely pathogenic, CADD
36) and `c.3006T>G` p.Asn1002Lys (missense in the kinase domain, absent from
gnomAD, CADD 24.5) — reproduce the canonical architecture of *BUB1B*-related
Mosaic Variegated Aneuploidy: one truncating allele plus one missense allele. An
independent hypothesis-driven analysis of the three known MVA genes converged on
the same pair.

**Strengths.** The result is demonstrated rather than asserted: two
methodologically independent routes agree, and the ranking survived a
substantially stricter filter (107,395 variants removed) with an identical
score, showing it never depended on tolerating noise. Every discard is counted
in a summary file, so no filter is silent. The submission CSV is generated
automatically from the pipeline's own output. Annotation runs against live
public APIs rather than frozen local caches, so the method stays current and
carries no additional data at rest. Everything is CPU-only and cost-free.

**Limitations.** Phase is not proven: the case is a singleton and the two
variants lie 10,911 bp apart, beyond the reach of short-read phasing even if the
FASTQ were reprocessed — the pair is *presumed* in trans. Only SNVs and indels
were assessed; no CNVs, structural variants or repeat expansions. The
HIGH/MODERATE impact filter is coding-biased and discards deep intronic and
regulatory variants. Phenotype ranking cannot score a gene that carries no HPO
annotations, however causal it might be. The second allele is the weaker one:
absent from gnomAD and computationally deleterious, but not present in ClinVar.

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

### References

- Hanks S. *et al.* Constitutional aneuploidy and cancer predisposition caused
  by biallelic mutations in *BUB1B*. *Nat Genet* 2004.
- Miller D.T. *et al.* ACMG SF v3.2 list for reporting of secondary findings.
  *Genet Med* 2023.
- Richards S. *et al.* Standards and guidelines for the interpretation of
  sequence variants (ACMG/AMP). *Genet Med* 2015.
- Resnik P. Using information content to evaluate semantic similarity. *IJCAI*
  1995.
- Stenton S.L. *et al.* Benchmarking of variant prioritisation in rare disease,
  2024 — scoring framework referenced by the challenge organisers.
- Köhler S. *et al.* The Human Phenotype Ontology.
- McLaren W. *et al.* The Ensembl Variant Effect Predictor. *Genome Biol* 2016.
- Cingolani P. *et al.* SnpEff. *Fly* 2012.
