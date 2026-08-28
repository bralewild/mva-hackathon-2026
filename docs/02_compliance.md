# Compliance with the Official Rules

Operational checklist derived from the *Official Rules* of the MVA Hackathon
2026. This document exists so that compliance is **verifiable** rather than
promised: every obligation has a concrete mechanism attached.

---

## 1. Privacy and recontact restrictions

| Obligation | How it is met |
|---|---|
| No recontact of the subject, their family, or the MVA Society | No pipeline stage queries identity, contact details, or social platforms. Only public scientific APIs (Ensembl, NCBI, HPO). |
| No release of, or access to, the data | Data exists only under `~/mva/` inside WSL. Never in the repository, never on an external service. |
| Safeguards against unauthorised use | Three barriers — see §3. |
| Report any unauthorised disclosure | Sage Help Center. |

**No patient data was uploaded to any service.** Queries to Ensembl VEP transmit
genomic coordinates (`15 40209701 T G`) with no subject identifier of any kind —
the intended use of an annotation API.

---

## 2. Data deletion — 30 days from close

The Official Rules are broader than the dataset gating form:

> *"All data must be deleted within 30 days of Hackathon close from all
> environments (local machines, cloud instances, notebooks, **private repos**,
> and any **intermediate or derived datasets**)."*

### What must be deleted

Run `pipeline/99_data_inventory.sh` for the exact inventory. Summary:

| Category | Path | Delete? |
|---|---|---|
| VCF, index, clinical document | `~/mva/data/raw/` | **Yes** |
| Extracted HPO terms | `~/mva/data/raw/patient_hpo.tsv` | **Yes** (derived) |
| Annotated VCF, candidates, VEP cache | `~/mva/work/` | **Yes** (derived) |
| Variant-level tables, per-variant reports | `~/mva/results/` | **Yes** (derived) |
| LINCS / Enrichr library cache | `~/mva/work/cache/` | No — public reference data, but it lives under `work/` and goes with it |
| Reference genome, snpEff DB, HPO ontology | `~/micromamba`, `~/mva/data/annot` | No — public resources |
| Pipeline code | git repository | No — contains no patient data |
| Committed aggregates and the published answer | `results/` in the repo | **See below — not a simple yes** |

### The repository copy needs a precise answer, not a blanket one

An earlier version of this table said `results/` in the repository was derived
data to be deleted, with `rm -rf "$PROJECT/results"` as the command. Both halves
were wrong, in opposite directions.

**It over-deletes.** Audited file by file, the committed `results/` contains
aggregate counts, pure drug-database output (DGIdb, STRING, LINCS — no patient
content whatsoever), the **causal variant** (`p.Leu737Ter`, `p.Asn1002Lys`,
`c.2210T>G`), and **three HPO terms** that appear in the challenge's own public
case description. It contains **no other patient variant and no coordinate from
the VCF beyond the submitted answer.**

The variant is the *published finding*. It is in the submitted CSV and in both
submitted reports, which the Official Rules release under **CC BY 4.0 with named
attribution** — the rules require publishing it. If naming `p.Leu737Ter` obliged
deletion, it would equally oblige deleting
`submissions/bralewild_track1_report.md`, which is a mandatory deliverable. That
reading is incoherent. The Data Use Agreement covers **the dataset and
variant-level derived tables**, not the result the hackathon exists to produce.

**And it does not delete.** `rm -rf` on a path that is committed to a public
repository removes it from the working tree only. Every file stays in git
history, in every clone and fork, and in GitHub's own caches. A deletion commit
*looks* like compliance without being it. If anything in this repository were
genuinely covered by the DUA, the remedy would be history rewriting
(`git filter-repo` or BFG), a force-push, and a request to GitHub support to
purge cached views — or taking the repository down. **Nothing currently
committed requires that**, and the safeguards in §3 are what keep it true.

### Command

```bash
# The data and everything derived from it - all outside the repository by design
rm -rf ~/mva/data/raw ~/mva/work ~/mva/results
```

The repository needs no deletion step. What it contains is the submitted,
CC BY 4.0-licensed result plus aggregate and public-database output. Re-run
`pipeline/99_data_inventory.sh` before notifying, to confirm that is still true.

### Notification

The official sources give **two different addresses**, so both are notified:

| Source | Address |
|---|---|
| Official Rules | `RarediseaserealkidMVAhackathon2026@synapse.org` |
| Dataset gating form | `MVAHackathon2026@synapse.org` |

> *"If you do not contact us, we may contact you directly in 30 days."*

---

## 3. Safeguards in place

Three independent barriers keep the child's genome out of GitHub:

1. **Physical separation.** Data lives on ext4 inside WSL; the repository is on
   NTFS. The files are not merely ignored — **they are not in the repository
   folder at all.**
2. **`.gitignore`.** Blocks `*.vcf*`, `*.bam`, `*.cram`, `*.fastq*`, `*.docx`,
   `patient_hpo.tsv`, `results/*.csv`, `results/*.tsv`, `results/*.html`,
   `submissions/bralewild_*`, `data/`, `work/`.
3. **Audit.** `99_data_inventory.sh` lists every location holding data and
   classifies it as patient-derived or public resource.

Verification before every push:

```bash
git ls-files | grep -E '\.(vcf|bam|cram|fastq|docx)$|patient_hpo'
# must return nothing
```

### Repository visibility — the operative rule

The repository **contains the causal gene and coordinates** in its README,
documentation, submission files and convergence report. Selectively redacting
them would leave an incomplete artefact that reads poorly to a reviewer.

**The repository is public**, published at the point of submission.

Two considerations were weighed:

* *Against publishing early* — while the leaderboard is open, a public
  repository exposes the method to competing teams.
* *For publishing early* — the challenge form requires a public repository URL,
  four teams already hold a perfect score (so the answer itself is not a
  secret), and submitting early verifies the CSV format while five attempts
  remain. Discovering a format defect on 24 October, with hours left, would be
  far more costly.

The second argument won. Note that this is a **strategic** decision, not a legal
one: the embargo explicitly permits publishing code at any time (see §5).

What is never published is raw patient data and the variant-level tables listing
the hundreds of other coordinates examined. That line is permanent and is
enforced by `.gitignore`, not by repository visibility.

---

## 4. Reproducibility — a functional requirement

> *"By registering, participants acknowledge that submissions **may be rerun**
> by the Hackathon organizers."*

The organisers can **execute this code**. Looking tidy is not enough:

- `pipeline/run_all.sh` runs the pipeline end to end
- every stage is **idempotent**: valid existing output is skipped
- stage 04 has an **incremental cache and resume** for network failures
- stage 02 validates **output content, not exit code** — snpEff exits 0 even
  when it fails
- the environment is declared in the README with exact versions
- `.gitattributes` enforces `LF` so scripts run on Linux
- all code and documentation are in English

**Requirement for anyone re-running:** they must hold their own authorised
access to the dataset. No patient data is bundled here by design; stage `00b`
extracts the HPO terms from the original `.docx`, which each participant must
obtain through the gated Hugging Face dataset.

---

## 5. Embargo and publication

| Action | Permitted? |
|---|---|
| Publishing **code, models and derived outputs** | **Yes, at any time** |
| Peer-review manuscripts during the embargo | No |
| Conference abstracts or posters | Yes, with **prior written approval** |

The embargo begins at hackathon close and ends when the organisers publish their
summary report or preprint.

**Practical consequence:** keeping the repository private until close is a
**strategic** choice, not a legal obligation.

---

## 6. Mandatory acknowledgement

Must appear **verbatim** in any publication, preprint, abstract or public
communication derived from the hackathon:

> "This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON
> (The Benchmarking, Evaluation, and Assessment Consortium for Science), with
> prize sponsorship from AWS and Anthropic. We are deeply grateful to the child
> and their family who generously contributed their data and their story to
> advance research into this rare disease. We acknowledge their trust in making
> this Hackathon possible."

In addition:

- **Subject privacy.** No publication may include information that could
  re-identify the patient or their family beyond what is already public through
  the family's own blog posts and communications.
- **Dataset citation.** Use the reference given on the hackathon Synapse page at
  the time of publication.
- **Ethics approval.** Protocol approved by **WCG IRB #20252010**.
- **Licence.** Submissions are released under **CC BY 4.0** with named
  attribution of the participant.

---

## 7. Deliverables per track

| | Track 1 | Track 2 |
|---|---|---|
| Submissions allowed | 6 | **1** |
| Prediction CSV | Yes | — |
| Written report | Yes | Yes |
| Public GitHub repository | Yes | Yes |
| 3-minute pitch video | See note | Yes |
| Evaluation | Automated (rank points + F-max) against the **NHS-validated** variant | Expert panel: rigour 35 %, impact 25 %, innovation 25 %, scalability 15 % |

> **Note on the video.** The Official Rules state generally that *"each team's
> submission includes a written report, a GitHub repository, and a 3-minute
> recorded pitch video"*, while the Track 1 submission tab asks only for
> CSV + report + repository. Given the ambiguity, it is safer to have the video
> ready.

---

## 8. Other requirements

- Participants must be 18 or older.
- Every member of a team must register **individually** and accept the rules.
  *(Individual participation as `bralewild` — not applicable.)*
- Compute is the participant's responsibility. *(Entirely local: WSL2 on an
  i9-14900HX, no cloud cost.)*
- Data may not be reshared through any channel.
