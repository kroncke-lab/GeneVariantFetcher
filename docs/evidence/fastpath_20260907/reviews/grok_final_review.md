I'll do a read-only adversarial review of the current diff and the two named files, staying within the six-call cap.Diff is in; next I’ll read the two new coverage files and the cohort helpers they depend on.I still need the classify/exclusion and caption-lookup implementations to check for lost counts and cache bugs.**Verdict:** Ship the cohort-guard changes. One caption-uniqueness bug can drop valid stamps; the audit has coordinate/source and bound issues. Do not relax clinical/mixed-cohort exclusions or fall through to paper tier. 30059973 is 185 carrier-bearing variant rows (not 185 carriers).

## Blocking

**1. Caption “conflict” is string inequality, not two tables** (`pipeline/table_cohort_phenotype.py` `_resolve_caption`, ~498–502)

Uniqueness is `len({normalized expanded caption}) > 1` → `ambiguous_source_caption` and `continue` for every row with that label. Identical repeats are tested; **TOC vs full caption, `### Table N` plus legend vs `Table N. …`, or a wrapped extra clause** are different strings and refuse the whole table.

Gold calibration can stay flat while this still bites papers with a contents line plus the real caption. **Fix:** if expansions share a normalized prefix, keep the longest; treat conflict only when neither is a prefix of the other (or when `classify_table_cohort` roles differ). Add that case next to `test_repeated_identical_caption_remains_eligible`.

Do **not** use the first hit, and do **not** classify from `Cases (n)` while the Table N identity is still ambiguous.

## Non-blocking but fix before relying on the audit

**2. Audit source ≠ derivation source** (`pipeline/extraction.py` ~8435–8478)

Derivation uses `_augment_pdf_linearized_tables(paper.full_text or prepared_full_text)`. The audit is passed the unaugmented string. Locators and candidates can miss linearized tables or point at different text. Pass the same `source_text`.

**3. Coordinates / bounds** (`pipeline/table_phenotype_coverage.py`)

- `split("\n")` + form-feed handling is correct (`\f` is not a line break; `strip()`/`split()` drop it). `\r\n` leaves `\r` but `strip()` still matches; untested.
- Caption is **only the heading line** (240 chars), not wrapped legend. `source_line_1based` is the heading, which is fine if consumers do not treat `caption` as the full table title.
- Clinical window is 12 lines or next heading; a long legend with headers below that is missed (under-flag; allowed).
- All heading hits are stored, then sliced to 20; `clinical_signals` slices to 12 with **no truncation flag**. For a huge supplement this is extra allocations, not a completeness claim. Cap the scan (e.g. stop appending after 20, still count) and flag signal truncation.
- `review_candidates` only when `candidates and missing`; empty candidate lists never mean complete. Do not later gate extraction on `clinical_table_candidate_count`.

**4. Caption cache**

Keyed on raw `label`, not `_squash`/`_normalize_label`. `"Table 2"` vs `"Table  2"` rescans. Cache `(normalized_label, bare)` and short-circuit `_resolve_caption` once `len(candidates) > 1`.

## Intentional (do not reverse)

Selected-count exclusions (`affected`, `symptom`, `ecg`, `diagnos`, `followup`, `without`, `negative`, `positive`, plus `_CAPTION_EXCLUDE_RE` on the header) correctly block 30059973 Table 14 (negative ECG / LQT3 / PCCD / overlap) and Table 11 (asymptomatic at dx vs MCE follow-up). Keep A/U null.

`patient_count_without_disease_context` correctly stops RYR2 40875405 T4 “Number of Patients” (review of 221 papers / 964 patients, 263 variant rows) and uncaptioned BRCA1 catalogs. **Do not** fall through to paper-tier ascertainment; that re-stamps mixed rosters. Preserve `Cases (n)` and disease-in-column (`Number of LQT2 patients`). `_disease_hit` breadth / target specificity is out of scope.

Literal partitions are unchanged (`test_selected_count_column_exclusions_leave_literal_partition_untouched`).

## Test gaps

- TOC/short vs long caption (finding 1).
- Audit gets linearized `source_text`; `\r\n`; 20+ candidates **and** 12+ signals.
- Deterministic 185-row paper: carrier column eligible, ECG/timepoint columns not; audit metadata only; no LLM.
- `model_used` prefix `router+` actually matches production.

Leave `docs/evidence/fastpath_20260907/` untracked. 567/1118 exact affected (88.78% new-value exactness) is below 90%; no headline/holdout/promotion.
