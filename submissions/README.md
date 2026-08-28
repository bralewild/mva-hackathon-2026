# Submissions

Track 1 deliverables, uploaded to the challenge submission form.

| File | Size | Purpose |
|---|---|---|
| `bralewild_blind-wgs-triage.csv` | 0.8 KB | ranked variant predictions |
| `bralewild_track1_report.md` | 15.3 KB | methods report (357 lines) |

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

The submission therefore prepends `chr`, with a second row in the VCF-native
naming at a lower EPCR as insurance. The hedge is free: extra rows do not reduce
the automated score, and F-max takes the maximum across confidence thresholds,
so the top row alone determines it when correct.
