# Manual curation protocol — BMPR2 holdout gold (20 papers)

## What this is & why
You are building the **ground-truth answer key** for 20 papers so we can measure
how accurately an automated pipeline extracted BMPR2 variants from them. For each
paper you record every BMPR2 variant and its patient counts. We then compare the
pipeline against your answers. **Work independently** — you are NOT given the
pipeline's answers, and you should not try to find them. Honest, careful reading
is the whole point.

## What you need
Nothing but this folder and a spreadsheet program (Excel or Google Sheets). The
full text of every paper is included in `papers/<pmid>.md`, so you do **not** need
journal access or the internet. (A PubMed link is in `manifest.csv` if you want to
see figures/tables in the original — optional.)

## The task (per paper)
1. Open `papers/<pmid>.md` and read it (skim intro/methods, read results + tables).
2. Find every **BMPR2 variant** reported in the paper's current patient/family
   cohort. Do not copy variants mentioned only in background text, prior studies,
   database annotations, references, or assay controls.
3. For each distinct variant, add **one row** to `curation_template.csv`:

| column | what to put |
|---|---|
| `pmid` | the paper's PubMed ID (already filled in for you) |
| `variant` | the variant **as written** in the paper. Prefer HGVS: cDNA `c.5946delT` or protein `p.Ser1982fs`. Copy exactly; don't normalize. |
| `curation_status` | `complete` only after this paper has been checked exhaustively; use `needs_review` while unresolved. |
| `germline_or_somatic` | `germline` if inherited/constitutional (blood, saliva, family, carrier, proband, germline testing); `somatic` if tumor-only (tumor tissue, ctDNA, cell line, LOH); `unknown` if unclear. |
| `carriers` | number reported for this variant. Never convert families to individuals. Blank if not reported. |
| `carrier_count_role` | `individual`, `family`, `unknown`, or `not_reported`. This prevents family counts from being scored as people. |
| `affected` | of the individual carriers, how many had **pulmonary arterial hypertension or the pulmonary vascular phenotype defined by the paper**. Blank if not explicitly reported. |
| `unaffected` | of the individual carriers, how many did **not** have **pulmonary arterial hypertension or the pulmonary vascular phenotype defined by the paper**. Blank if not explicitly reported. |
| `evidence_note` | where you found it (e.g. "Table 2") + a short quote. |

## Rules that matter
- **One row per distinct variant.** If a variant appears in two cohorts in the
  same paper, sum the carriers into one row (note it).
- **Germline is the target.** Still record `somatic` variants (we use them to
  calibrate), but mark them clearly — they are not heritable carriers.
- **Never guess a number.** If a count isn't stated, leave it blank. A blank
  means "not reported", which is different from 0.
- **No BMPR2 variant in the paper?** Keep the paper's row and put `NONE` in
  `variant`, `complete` in `curation_status`, and explain the source check in
  `evidence_note`. (This matters — it tells us the pipeline shouldn't have found any.)
- **Family counts are not carrier counts.** Record the number with
  `carrier_count_role=family`; the scorer excludes it from person-level MAE.
- **Don't normalize notation.** Copy the variant string verbatim; we handle
  matching. If the paper gives both cDNA and protein, put both (one in `variant`,
  the other in the note).

## Worked examples
| pmid | variant | curation_status | germline_or_somatic | carriers | carrier_count_role | affected | unaffected | evidence_note |
|---|---|---|---|---|---|---|---|---|
| 12345678 | c.5946delT | complete | germline | 14 | individual | 9 | 5 | Table 2; "14 carriers, 9 with phenotype" |
| 12345678 | p.Cys1365Tyr | complete | germline | 1 | individual | 1 |  | Case report, Results para 3 |
| 99999999 | NONE | complete |  |  | not_reported |  |  | Results/tables checked; no current-cohort variant |

## When you're done
Save the spreadsheet as `curation_template_FILLED.csv` (keep CSV format) and send
it back. Budget ~3–4 days for 20 papers. Thank you!
