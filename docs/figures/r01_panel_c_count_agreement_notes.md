# Panel C - conditional carrier-count agreement

## What the current repository can support

The locked cardiac gold-118 evaluation supports count comparison at the
**matched variant-paper-measure** grain. It does not contain explicit cohort or
phenotype-definition identifiers, so it cannot yet support the originally
requested variant-cohort-phenotype claim.

The figure therefore separates two questions:

1. **Coverage:** how many gold-asserted count rows have both a manual and an
   automated integer?
2. **Conditional agreement:** among those paired rows only, how close are the
   two integers?

This separation is essential. Agreement on supplied pairs is not end-to-end
count recovery.

## Suggested manuscript caption

**C | Carrier-count agreement among paired variant-paper rows.** The analysis
uses the locked, source-bound gold-118 run. Each plotted observation is
a matched variant-paper row for which both the frozen manual gold standard and
the automated extraction supply an integer for the indicated measure; explicit
cohort and phenotype-definition identifiers are unavailable. Paired values were
available for 56 of 632 gold-asserted affected rows (8.9%) and 28 of 631
gold-asserted unaffected rows (4.4%). Among unpaired rows, 498 affected and 525
unaffected values are AI nulls on matched variant identities; 78 rows for each
measure are variant-identity misses. The denominators differ by one because RYR2
c.3599-9delT (PMID 30403697) has an adjudicated affected count of 1 and an
intentionally null unaffected count. Bubble area represents the number of rows
at each manual/automated coordinate, with circles for affected counts and
diamonds for unaffected counts. Coral outlines denote non-exact coordinates.
Among paired rows, 47/56 affected counts (83.9%) and 22/28 unaffected counts
(78.6%) were exact; mean absolute difference was 0.32 carriers for each
measure. Six affected and four unaffected rows differ by exactly one carrier;
the three affected differences exceeding one carrier and the two unaffected
differences exceeding one carrier are labeled.
H562R (PMID 25819988) exceeds one carrier for both affected and unaffected
counts; the two labeled K897T differences arise from different papers. Twenty-one
variant-paper rows contribute to both facets, so the facets are not independent
replicates. The displayed differences were not re-adjudicated specifically for
this figure.
These results measure numeric agreement conditional on count supply and variant
matching; they do not measure cohort-level validity or overall count-recovery
recall.

## Metric definitions

- **Paired coverage:** rows with a manual integer and an automated integer,
  divided by all gold rows with an asserted value for that measure.
- **Exact fraction:** paired rows where automated count equals manual count,
  divided by all paired rows.
- **Mean absolute difference:** mean of
  `abs(automated count - manual count)` over paired rows.
- **Absolute difference = 1:** paired rows differing by exactly one carrier.
- **Absolute difference > 1:** paired rows differing by more than one carrier;
  this is descriptive, not a prespecified clinical tolerance.

No correlation, regression slope, or unclustered confidence interval is shown.
Rows from the same paper are not independent, 21 rows contribute to both
facets, and a paper-clustered inferential analysis has not yet been
preregistered or run.

## Experiments needed for the original exact-grain claim

### 1. Build a cohort-aware gold-standard extension

Extend every count-bearing gold row with:

- stable `cohort_id` and `family_id` fields;
- a controlled `phenotype_definition_id`;
- count semantics (`person`, `family`, `allele`, or cohort total);
- explicit-zero versus unreported status;
- source location and source type (main text, table, supplement, or figure);
- overlap links for cohorts reused across papers or publications.

Have two curators annotate independently, followed by blinded adjudication.
Report agreement before adjudication and retain both original annotations.

### 2. Run a frozen, paper-held-out extraction evaluation

Freeze the extraction code and prompts before exposing the new answer key.
Split by paper - not by variant row - while balancing gene, source type, pedigree
content, table/supplement dependence, and count magnitude. Primary endpoints
should be paired-count coverage, exact fraction, and MAE at the
variant-cohort-phenotype grain. Report zero/nonzero disagreements separately.

### 3. Run a source-representation ablation

On the same held-out papers, compare:

- main text only;
- main text plus structured tables and supplements;
- full multimodal input including pedigree/figure interpretation.

This experiment identifies whether missing counts arise from retrieval,
representation, or extraction. Analyze paired differences within paper and use
paper-clustered uncertainty.

### 4. Create a targeted leakage and duplication stress set

Prospectively select difficult papers with multi-variant cohorts, pedigrees,
aggregate patient tables, repeated families, and overlapping publications.
Score cohort-total leakage, duplicate-cohort counting, phenotype reassignment,
and person-versus-family errors as separate adjudicated outcomes.

### 5. Repeat the final analysis on a sealed external set

After all rules are frozen, evaluate a newly curated paper set not used for
pipeline development. Preregister acceptance thresholds from clinical and
curation requirements rather than choosing them after observing the results.
Release the row-level audit table and paper-clustered uncertainty with the
figure.

## Reproducibility

Regenerate the editable SVG and its data files with:

```bash
.venv/bin/python scripts/build_panel_c_count_agreement.py
```

The script verifies the locked selection and prediction hashes, reproduces the
benchmark scorer's matching order, checks row counts and MAE against
`report.json`, and writes:

- `docs/figures/r01_panel_c_count_agreement.svg`
- `docs/figures/data/r01_panel_c_count_agreement.csv`
- `docs/figures/data/r01_panel_c_count_agreement.json`

The PDF is a vector conversion of the SVG; the PNG is a 2160 × 1350 raster
export (300 dpi at 7.2 × 4.5 inches).
