# Submissions

Deliverables uploaded to the challenge submission forms.

## Track 1 — variant prediction · **SUBMITTED**

> **Scored 100.0 / 100** — F-max 1.000, full match at rank 1.
> 5 of 6 attempts remain.

| File | Size | Purpose |
|---|---|---|
| `bralewild_blind-wgs-triage.csv` | 1.7 KB | ranked variant predictions (7 lines) |
| `bralewild_track1_report.md` | 26.6 KB | methods report (555 lines) |

## Track 2 — drug repositioning · **NOT YET SUBMITTED**

> **One attempt only.** No resubmission, no correction after the fact.

| File | Size | Purpose |
|---|---|---|
| `bralewild_track2_report.md` | 33.9 KB | mechanism → candidate report (659 lines) |

Supporting evidence for Track 2 is in [../results/track2_evidence/](../results/track2_evidence/)
so the negative result can be audited rather than taken on trust.

Both are committed here. The Official Rules release submissions under
**CC BY 4.0**, and a repository with its results redacted reads as incomplete to
a reviewer.

> **What is not here:** raw patient data, and variant-level tables listing the
> hundreds of other coordinates the pipeline examined. Those fall under the Data
> Use Agreement and are blocked permanently by `.gitignore`. See
> [../docs/02_compliance.md](../docs/02_compliance.md) §2–3.

## Naming convention (required by the challenge)

```
<hf-username>_<model-name>.csv     ->  bralewild_blind-wgs-triage.csv
<hf-username>_track1_report.md     ->  bralewild_track1_report.md
<hf-username>_track2_report.md     ->  bralewild_track2_report.md
```

## CSV format (verified against `tabs/submit_track1.py` and `evaluation.py`)

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
wrong naming convention scores zero.**

The submission therefore prepends `chr`. Three independent sources confirm that
convention: the field description in `tabs/submit_track1.py` (*"e.g. chr15"*),
the `_LOCAL_FALLBACK` in `groundtruth.py` (`("chr2", ...)`, `("chr15", ...)`),
and the official `track1_submission_template.csv` (`chr1`, `chr7`).

An earlier draft carried a second row repeating the pair in VCF-native naming as
a hedge. It was dropped: with six attempts and instant scoring, a resubmission is
better insurance than an in-file duplicate, and a single confident row reads
better than one that signals uncertainty about the format.
