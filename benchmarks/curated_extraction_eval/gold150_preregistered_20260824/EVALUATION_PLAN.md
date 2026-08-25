# Preregistered evaluation plan

## Frozen cohort and split

- Cohort: the exact `*_pmids.txt` rosters under
  `results/vb_curated150_20260824/<GENE>/20260824_181400`.
- Genes: BRCA1, BRCA2, BMPR2; 50 papers per gene.
- Assignment: sort by
  `sha256("20260824|<GENE>|<PMID>")`; first 30 are calibration and the
  remaining 20 are holdout.
- Every selected full-text representation is bound by SHA-256 in `manifest.csv`.
- Calibration and holdout packets are separate. A calibration curator is not
  given the holdout packet or roster.

## Human answer key

The curator works from the included source, not pipeline predictions. Each paper
must be marked `complete` and contain either an exhaustive set of current-cohort
target-gene variants or one explicit `NONE` row. LLM evidence cards may help
locate passages but may not supply, approve, or edit the answer key.

Semantics are strict:

- blank count: not reported;
- `0`: the source explicitly reports zero;
- `NONE`: no in-scope target-gene variant after checking the current cohort,
  tables, and supplements;
- somatic: recorded for audit but excluded from the default heritable-carrier
  target; its PMID remains in the precision scope;
- family count: recorded with `carrier_count_role=family` and excluded from
  person-level carrier MAE;
- affected/unaffected: only explicit individual-level counts for the phenotype
  named in the packet protocol. No subtraction or phenotype inference.

The scorer rejects blank variants, incomplete status, missing evidence notes,
mixed `NONE` and variant rows, invalid count roles, and any PMID roster mismatch.

## Calibration and holdout use

1. Curate and score the 90 calibration papers (30 per gene).
2. Any code, prompt, matcher, or trust change may use calibration residuals.
3. Freeze one candidate commit and its runtime configuration.
4. Release and curate the 60 holdout papers (20 per gene).
5. Score holdout exactly once. Do not tune on holdout PMIDs or select the best
   model/run after seeing their results.

If a source is insufficient for exhaustive adjudication, record the paper as
`needs_review`; the scorer will refuse to run. Replace neither the paper nor the
split after seeing predictions.

## Metrics

Report each gene separately before any aggregate:

- variant-row and unique-variant recall;
- paper-exhaustive raw precision over the explicit PMID roster, including
  `NONE` and somatic-only papers;
- F1 and recall-weighted F2;
- carrier, affected, and unaffected supplied-value coverage;
- matched-row MAE/RMSE and end-to-end count error;
- count-role exclusions and explicit-zero strata.

The counted-extra precision view is secondary on this exhaustive gold. Do not
substitute it for raw paper-exhaustive precision.

## Cost and execution

Packet construction and scoring are deterministic and use no model calls. Do not
pay for another exact-150 extraction merely to create gold. After calibration,
human adjudication, and code freeze, one sealed holdout extraction may be run if
the existing frozen DB is not the candidate being evaluated. Its spend must be
metered separately and remain inside the active `$100` improvement envelope.

Scoring must pass the packet's exact PMID scope and exhaustive flag through
`score_curation_packet.py`; review-gold sync is disabled and the `all` gene tier
is selected so BRCA1, BRCA2, and BMPR2 cannot be silently filtered out.
