# Affected/unaffected failure panel: 22 papers, 2026-09-05

The next improvements should target **specific missing source components,
unconsumed phenotype tables, and the meaning of the counted endpoint**. More
whole-paper extraction alone will not resolve the largest discrepancies.

I ranked **336 completed gene/run attempts** across four opened locks and chose
**22 unique papers / 30 gene/run attempts**. Fourteen papers have at least some
supplied affected/unaffected counts in an inspected run; eight have none. Two
of the fourteen are exact positive controls. The selection includes large
count omissions, wrong supplied values, and high identity FP/FN papers that
count-error ranking would miss because their phenotype gold is zero.

There is also a concrete acquisition result: full papers for **SCN5A 20031634
and 25163546 were recovered from author repositories and added to the corpus**
through `build_source_corpus.py` after its dry run. For 20031634, Table 1 was
visually verified and its 13 source rows reconstructed into an audit CSV. This
is recovered source evidence, **not a new extraction score or accepted protocol**.

## What the numbers mean

- **TP/FP/FN** refer to variant identities against the fixture. An unmatched
  prediction is not automatically a biologically false statement: linkage,
  engineered variants, historical observations, and incomplete gold need review.
- **A/U supplied** counts populated affected/unaffected *fields on matched gold
  rows*, not people. A supplied zero counts as supplied; an NA does not.
- **A/U error** is the sum of absolute differences, separately by field. A
  missing identity or count evaluates as zero for this diagnostic only. Stored
  NULLs remain NULL. A missing prediction against gold zero has zero magnitude
  but remains explicitly missing, never a successful count capture.
- Gold-v2 overrides are honored. Gold zero can mean unreported; different gold
  papers also use different meanings of affected. Discrepancy units are **not
  proved misclassified people or recoverable yield**.
- `G118` = `20260902_false_zero_recovery_gold118`, legacy trusted projection,
  including linkage. `M01`, `M02`, and `C120` = the opened mixed01, mixed02,
  and cont120_01 candidate locks, respectively, paper-derived primary lane.
  Never pool their denominators or interpret their differences as treatment
  effects. All selected occurrences are retained in [panel.csv](panel.csv).

## The selected papers

The table shows one representative gene/run per paper. Additional selected
gene/run occurrences are in the CSV, not added together here. Source locators
refer to `corpus/<GENE>/<PMID>/<PMID>_FULL_CONTEXT.md` unless stated otherwise.

| PMID / gene / run | TP/FP/FN | A/U supplied | A/U error | Source finding and next action |
|---|---:|---:|---:|---|
| **30059973 SCN5A M02** | 182/3/3 | 0/0 | 445/0 | Full manuscript **already includes** supplemental Tables 5, 11, 14. Fixed-width variant parsing returns before phenotype-table reading. Preserve the symptom, ECG, and follow-up endpoints separately. |
| **14678125 KCNQ1 G118** | 35/9/6 | 0/0 | 289/0 | Both inspected runs used only an abstract. Current longer Wiley object is abstract + references/navigation, without the mutation table. M01 is 0/0/41. Retrieve the actual body/table; legacy linkage matches do not establish source coverage. |
| **20031634 SCN5A C120** | 0/0/13 | 0/0 | 60/48 | Abstract at extraction. **Now recovered from HAL.** Table 1 directly separates carrier ECG+, carrier ECG−, unknown phenotype, and ECG+ noncarriers. Replay on the recovered table after endpoint adjudication. |
| **33673806 MYBPC3 M01** | 73/1/1 | 55/0 | 93/0 | Additional file 2 is already converted. c.3181C>T = 24 and c.2373dup = 14 reach extraction but affected values are cleared as copied carriers. Study participants have suspected/presumed HCM, not uniformly verified diagnoses. Apply one ascertainment rule to singletons and larger groups. |
| **26496715 KCNQ1 G118** | 41/6/1 | 0/0 | 82/0 | `Extanded tables.doc` is already converted; its `No. of patients` column is present. **Body missing.** Acquire methods/ascertainment before deciding whether those counts support affected status. |
| **25814417 RYR2 G118** | 1/0/0 | 1/1 | 24/44 | G357S supplied 97/62 with 26 uncertain; gold 73/106. Supplemental patient rows are present and the deterministic audit already runs. Review endpoint/timepoint and unknown handling; do not force a gold-matching partition. |
| **19926015 RYR2 G118** | 35/4/5 | 1/1 | 45/13 | Table 2 mixes current cases, controls, and previously reported variants; supplemental material includes conservation rather than the missing clinical split. Bind study ownership and distinguish persons from control alleles/frequencies. |
| **22677073 RYR2 G118** | 21/5/4 | 2/0 | 41/0 | Molecular-autopsy case tables and exclusion lists are on disk. The same article has differing capture across KCNQ1/SCN5A/RYR2. Preserve case identity and study phenotype; do not infer syndrome diagnosis from every sudden-death case or turn exclusion-list cases into family counts. |
| **18808722 KCNQ1 G118** | 1/1/0 | 0/0 | 13/18 | Text gives 31 L187P carriers and 18 with normal-to-borderline initial QTc. Borderline is not proven unaffected; source defines an uncertain category. KCNH2 in the same paper supplies two affected values. Read definitions, table and pedigree jointly. |
| **21302287 MYBPC3 C120** | 25/5/3 | 11/0 | 29/0 | Main Table 3b gives explicit patient IDs, including `278, 251 268, 194`. IDs survive data-zone selection and appear in model requests; two inspected intronic rows end as identity-only regex records. Join distinct IDs to clinical observations, including synonymous/intronic variants. |
| **19398665 RYR2 G118** | 27/1/0 | 0/0 | 27/0 | Current source has Table 1/2 captions with no table bodies. Figure inventory has graph/topology assets; table entries have no image URL. Obtain actual clinical tables. Cohort-wide baseline/follow-up totals cannot be assigned to each mutation. |
| **30403697 RYR2 G118** | 21/0/0 | 6/2 | 17/9 | Table 1 has two variant columns, patient rows, shared families, and silent parents in prose. Counts are read but verifier language switches between CPVT enrollment and symptoms; some are then absent in the trusted projection. Build explicit variant–person–endpoint observations and trace their survival. |
| **15671429 SCN5A C120** | 4/1/0 | 2/0 | 8/15 | Figures are on disk; D1275N phenotype counts are still absent. G118 supplied none (4/8/0 identities). Audit genotype annotations and the DCM/AF legend per person. The historical 23/2 split is **not** the current prediction. |
| **25163546 SCN5A C120** | 0/0/20 | 0/0 | 20/0 | Abstract at extraction; current FULL_CONTEXT was empty. **Now recovered UCL manuscript.** Two alleged supplement files are HTML for the EHJ Supplements journal and OUP advertising, not this article’s data. Find article-bound variant supplements before claiming count recovery. |
| **19632626 KCNQ1 G118** | 1/63/0 | 1/1 | 6/5 | Largest identity FP paper; legacy citation linkage explains the large extra set in the prior lane audit. S209P source is familial AF, supplied 6/1 versus gold 0/6. Separate endpoint adjudication from linkage precision; this is not evidence to delete all extra variants. |
| **20181576 KCNH2 G118** | 2/1/0 | 2/2 | 4/6 | K897T and P926fsX splits disagree with gold. Text/figures describe different QT and symptom subsets, including modifier co-carriage. Reconcile variant ownership and target endpoint before altering values. |
| **33315912 RYR2 G118** | 1/0/0 | 1/1 | 2/5 | P2328S supplied 17/42 with three uncertain, gold 15/47. Source-bound patient derivation already exists. Review baseline versus follow-up and uncertain cells, not another blanket guard. |
| **25819988 KCNH2 G118** | 1/0/0 | 1/1 | 3/3 | H562R: 13 carriers, six symptomatic, nine diagnosed LQTS. Supplied 6/7 versus gold 9/4. Source provenance supports six symptoms but the inspected extraction has no unaffected fact supporting seven. Preserve endpoints and investigate the unsupported complementary value. |
| **29650123 KCNH2 G118** | 1/0/21 | 0/0 | 0/0 | High FN hidden by zero phenotype gold. Prior publisher audit found no mutation roster in the available supplement. Obtain an author/provenance-valid roster; do not reconstruct missing identities from gold. |
| **32533946 SCN5A G118** | 83/44/0 | 0/0 | 0/0 | Functional reclassification paper with many reference/assay identities and no supplied clinical counts. Keep as a false-count negative control and adjudicate study-observed versus reference/engineered identities. |
| **19398417 RYR2 G118** | 1/0/0 | 1/1 | 0/0 | Positive control: W4645R **4/2/2** is explicitly recoverable even from the abstract. Do not impose a blanket full-text-only or figure-only refusal. |
| **25435091 RYR2 G118** | 1/0/0 | 1/1 | 0/0 | Positive control: C2277R **8/7/1** survives the existing patient-row derivation and trust path. Preserve this while extending more complex table handling. |

## Download the component that is missing

**First priority: component-level acquisition status.** The previous identity
presence sweep cannot answer whether phenotype evidence was acquired. Replace
the operational inference “file exists / variants found / supplements downloaded”
with a component record: article identity, body, named table/figure, actual
payload type, source URL, download hash, converted block, row/header/legend
availability, and whether extraction consumed it. A table can be embedded in
an author manuscript; one supplement file can contain many tables.

Two concrete weaknesses are now reproducible:

1. `harvesting/supplement_reference_parser.py::parse_supplement_references`
   sees only `Table S14` in the rich 30059973 context: it misses the repeated
   **Online Table** references. Its gap logic compares table counts to file
   counts, which are different units. Resolve component IDs to assets/blocks.
2. The two HTML objects under 25163546’s supplements directory are journal and
   advertising pages. Validate parent DOI/title and supplementary-container
   linkage, MIME/magic bytes, and meaningful payload before calling a download
   a supplement. Existing `supplement_processing_service.py` identity checks
   are a foundation, but PDF checks alone do not address these HTML objects.

**Second priority: repository fallback with publication/version validation.**
The failed 20031634 abstract contains its DOI despite a failure reason saying
“No DOI or URL.” Resolve bibliographic identity before concluding no route
exists. Title/DOI search found the matching HAL author deposit; the publisher
returned 403 and CiteSeer timed out, while HAL returned the actual paper. UCL
provided 25163546’s manuscript. A similarly titled conference abstract has a
different DOI and different numbers: it is not an interchangeable substitute.

The two recovered bodies are deliberately marked `body_only`/retry in corpus
status, using the existing incomplete-supplement marker. The article-specific
supplement surface remains unresolved. Original PDFs and conversion records are
under `validation_runs/phenotype_failure_panel_20260905/`; hashes and URLs are
in [recovered_sources.json](recovered_sources.json) and the acquisition logs.

## Capture counts from material already acquired

**1. Complete phenotype-table processing after identity-table success.**
`pipeline/extraction.py` returns an `ExtractionResult` after the large
fixed-width table shortcut. For 30059973 this yields 185 extracted variant
records and bypasses LLM reading; the metadata explicitly says so. Tables 11
and 14 are still in FULL_CONTEXT and CLEANED, but the markdown router sees
zero tables. Finish the identity phase, then inspect any unconsumed clinical
components. Do this only for identified components and existing identities,
with bounded cost and no invented counts.

Preserve a PDF page’s orientation, cells, hierarchical headers, caption and
footnotes. `harvesting/format_converters.py::pdf_to_markdown` currently returns
on a successful text converter before its table-aware fallback; nonempty text
is therefore not proof of a faithful table. For 30059973 I visually checked
PDF pages 50 and 54 (supplement printed pages 24 and 28). Table 11 has a
**mutation presence/absence** column hierarchy, while Table 14 is landscape
with seven phenotype categories. Treating every integer nearby as a carrier
or affected count would be incorrect.

**2. Support explicit variant-to-patient joins as a separately validated
extension.** The current `patient_row_phenotype.py` intentionally requires
unique contiguous IDs in one complete variant-specific table. It does not
license arbitrary joins, noncontiguous lists, or two variant columns. Add a
separate path with a source-defined ID namespace, unambiguous table joins,
per-person genotype and endpoint evidence, duplicate detection, and an audit
for every included/excluded row. Test on 21302287 and 30403697. A patient can
carry two variants; count that patient once per appropriate variant, never
once per appearance or family mention. Related prose must establish which
parent carries which variant; “parents” alone cannot be copied to both alleles.

**3. Make endpoint scope survive parsing, verification and trust.** Preserve
carrier count, symptoms at presentation, diagnostic phenotype, ECG class,
follow-up events, and unknown status as distinct observations. Only project to
affected/unaffected when the declared target and evidence agree. This does not
require replacing the whole pipeline at once: first retain these facts for
the inspected components and expose unsupported mappings for review.

For 33673806, [stage_probes.json](stage_probes.json) records affected facts of
24 and 14 followed by `copied_carriers_onto_affected` flags; the values were
not lost during download. A fix that simply stamps them as source-backed would
still fail to resolve suspected versus confirmed HCM. For 25819988, six
symptoms do not establish seven clinically unaffected carriers. For 25814417
and 33315912, do not erase the 26 and three uncertain observations to improve
agreement with legacy gold. Stable gold-discrepancy labels should remain
separate from any curator-approved corrected endpoint standard.

## Acquisition proof: PMID 20031634 Table 1

The recovered [HAL paper](https://univ-brest.hal.science/hal-00750425/document)
matches the title and DOI on its cover and article first page. PDF page 3 /
printed page 553 supplies 13 family rows. Inspection of the original page and
coordinate-based transcription agree on totals: **115 carriers**, comprising
**54 BrS-ECG+, 55 not BrS-ECG+, and six undetermined**. Eight ECG+ noncarriers
are a different column and excluded from the carrier partition. G1408R is explicitly
14 = 4/9/1. [recovered_20031634_table1.csv](recovered_20031634_table1.csv)
retains original labels and raw text-layer notation, PDF hash, page and row
coordinates. Glyph-damaged cDNA is not silently relabeled as valid HGVS.

This establishes that appropriate count-bearing material is downloadable.
It also demonstrates why retrieval and score improvement are separate: this
lock’s gold sums to 60 affected and 48 unaffected, not the same partition.
The next replay must preserve the six unknowns and review per-row endpoint
disagreements. No database counts were written and no scored lock was changed.

## Reviewer challenge and disposition

Both requested CLIs were used. Grok’s first file-review attempt exhausted its
turn limit; a bounded facts-only review completed. Agy’s first headless read
was permission-denied; its authorized read-only review then completed, followed
by a facts-only challenge. The final [Grok review](grok_review.md),
[Agy review](agy_review.md), and [shared brief](reviewer_brief.txt) are retained.
The final CLI calls used their configured models; Grok’s model listing named
`grok-4.6` as default. No paid production extraction arm was opened.

- Accepted: separate acquisition from component consumption; endpoint-labelled
  evidence; unknown preservation; source-grounded joins; no blanket N→affected
  promotion or “symptomatic”/“unaffected” synonym rule.
- Rejected Agy’s initial unverified paper-level diagnoses: 26496715 and
  33673806 already have converted supplements; 14678125 had an abstract, not a
  failed table router; 15671429’s old 23/2 prediction is not in these locks.
- Qualified Agy’s claim that endpoint mismatch is the primary driver: this
  intentionally enriched panel does not establish its population prevalence.
- Qualified Grok’s “do not score” missing-source papers: retain them in the
  **end-to-end** denominator and report a separate source-ready capture view.
  Likewise unconsumed tables are count-capture failures, even when their bytes
  were acquired successfully.
- Both reviews preceded successful HAL recovery. Their statement that
  20031634 remains unverified is superseded by the PDF check above.

## Next implementation and acceptance order

1. Make source component status reliable (shells, wrong supplement HTML,
   missing clinical table bodies, and table references inside PDFs). Use
   14678125, 19398665, 25163546 and the recovered 20031634 as fixtures.
2. Replay a general clinical-table continuation after deterministic identity
   extraction on the fixed 30059973 PDF. Preserve Table 11/14 endpoint labels;
   **do not** optimize toward affected=69 merely because gold uses it.
3. Prototype audited patient joins on 21302287/30403697, with the two exact
   controls and adversarial duplicated IDs, second variant columns, missing
   genotype, ambiguous parent links, and incomplete status cells.
4. Validate component retrieval separately from count extraction. For each
   opened cohort report identity TP/FP/FN, A/U supply, supplied exactness,
   wrong supplied values, and all-row error including omissions. Keep
   source-only gold conflicts visible. After a candidate is frozen, use a
   separately authorized unopened cohort for confirmation; this panel is
   permanently calibration evidence.

## Artifacts and verification

- [all_papers.csv](all_papers.csv): all 336 ranked attempts.
- [panel.csv](panel.csv), [panel_rows.csv](panel_rows.csv): 30 selected
  attempts and 2,321 asserted carrier/phenotype field comparisons.
- [extra_predictions.csv](extra_predictions.csv): 164 unmatched identities
  across those attempts, with counts/provenance; no invented gold values.
- [source_inventory.csv](source_inventory.csv): pre-recovery current-corpus
  and run-file inventory, kept separate; 655 observed files with hashes.
- [stage_probes.json](stage_probes.json): source, request and extraction
  probes, count flags, decision metadata and input hashes for five examples.
- [manifest.json](manifest.json): prediction/selection lock checks,
  report/gold hashes, and per-cohort count statistics. It reproduces the locked
  TP/FN, supplied counts and conditional MAE for all three fields in all four
  runs. Inventory hashes describe observed audit-time bytes; they are not
  claimed to be independently bound to historical source locks.
- [verification.json](verification.json): independent rebuild reproduced all
  five diagnostic outputs byte-for-byte and all cohort metrics. All 34 input
  hashes still match; both recovered corpus bodies match their inspected source
  conversions. The rebuilt inventory intentionally reflects the new sources.

Rebuild into a **new output directory** to preserve the pre-recovery inventory:

```bash
.venv/bin/python scripts/recall_audit/phenotype_failure_panel.py \
  --panel docs/evidence/phenotype_failure_panel_20260905/panel_selection.json \
  --out-dir validation_runs/phenotype_failure_panel_20260905/recheck
```

Validation: Ruff passed; the four new denominator/NULL/lane tests and eight
existing difference-figure tests passed (**12 tests**). Both PDF table layouts
were inspected for 30059973; 20031634’s table was checked visually and all 13
carrier partitions and five printed column totals reconciled. Production
code, gold, locked predictions, accepted metrics and recovery defaults are
unchanged.
