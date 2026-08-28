# Submissions

Prediction files (`*.csv`) and the methods report are **not committed to this
repository** while the competition is open (close: 24 Oct 2026, 23:59 UTC).

**Why:** the challenge requires a public repository. Publishing the causal
coordinates before close would hand the answer to competing teams. What is
judged here is the **method** — reproducibility, innovation and scalability —
not a pre-chewed result.

The files are uploaded directly to the challenge submission form, which is
private until evaluation. They will be added to this repository after close so
the record is complete.

> **Operative rule:** this repository stays `PRIVATE` until the moment of
> submission. Its README, documentation and mirrored result files name the
> causal gene, so visibility — not selective redaction — is the protection.
> See [../docs/02_compliance.md](../docs/02_compliance.md) §3.

## Naming convention (required by the challenge)

```
<hf-username>_<model-name>.csv     ->  bralewild_blind-wgs-triage.csv
<hf-username>_track1_report.md     ->  bralewild_track1_report.md
```

## Format (verified against `tabs/submit_track1.py` and `evaluation.py`)

| Field | Type | Note |
|---|---|---|
| `proband_id` | string | **must be `PROBAND01`** — hardcoded in the submission handler |
| `chrom_1` / `chrom_2` | string | **`chr` prefix required** (e.g. `chr15`) |
| `pos_1` / `pos_2` | int | GRCh38 |
| `ref_*` / `alt_*` | string | |
| `epcr` | float | range `(0, 1]`, rows sorted descending |
| `finding_type` | string | `primary` or `secondary` |
| `notes` | string | optional rationale |

Maximum 10 rows. Maximum **6 submissions** per participant; only the
highest-scoring one appears on the leaderboard.

## The normalisation trap

`evaluation.py` performs **no normalisation** of chromosome names. It compares
exact tuples:

```python
(chrom.strip(), int(pos), ref.strip().upper(), alt.strip().upper())
```

The source VCF uses Ensembl contig naming (`15`), while the evaluator's own
fallback and the submission documentation use `chr15`. **A correct answer in the
wrong naming convention scores zero.** The submission therefore prepends `chr`,
with a lower-EPCR hedge row in the VCF-native naming as insurance — extra rows
do not reduce the automated score, and F-max takes the maximum across
thresholds, so the hedge costs nothing.
