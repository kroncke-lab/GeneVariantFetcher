### Advisory Evaluation: GeneVariantFetcher Phenotype Audit

---

### 1. Challenge to Previous Assumptions & Conclusions

1. **Prior Hypotheses Contradicted by Direct Evidence**:
   * The assumption that discrepancies were caused by unparsed supplements or stale pedigree errors is disproven.
   * **SCN5A 30059973**: Full 70-page text and supplements were already on disk.
   * **KCNQ1 26496715**: "Extanded tables.doc" was already converted to 382 Markdown lines (`No. of patients` available); missing text is the body ascertainment, not the supplement.
   * **MYBPC3 33673806**: Additional file 2 was already converted (`c.3181C>T` [24], `c.2373dup` [14]).

2. **Endpoint Mismatch is the Primary Error Driver**:
   * Extraction errors frequently stem from treating different clinical endpoints as equivalent.
   * **SCN5A 30059973**: Gold has Affected=69 / Unaffected=0 (carrier N). Main tables report presentation symptoms (0 cardiac arrest, 9 syncope, 60 asymptomatic) and ECG types (29 negative, 13 LQT3, 17 PCCD).
   * **KCNH2 25819988**: 13 carriers, 6 symptomatic, 9 LQTS diagnosed. Predicted 6/7 vs. gold 9/4 reflects distinct clinical endpoints, not parsing failure.

3. **Forbidden Heuristics**:
   * **No Blanket Carrier Promotion**: Copying total carriers to "affected" invalidates distinct phenotype extractions.
   * **No Residual Subtraction**: Unaffected $\neq$ Total Carriers $-$ Symptomatic.
   * **No ID-to-Count Coercion**: Main text strings like `278, 251 268, 194` (**MYBPC3 21302287**) are patient IDs, not clinical totals.
   * **Gold Zero $\neq$ Healthy**: Gold zeros frequently indicate unmentioned or uncharacterized phenotypes.

---

### 2. Top 3 General Implementation Steps

```
[Component Manifest] ───► [Atomic Observation Extraction] ───► [Deterministic Replay Engine]
(Register Body/Tables)     (Gene, Variant, Subject, Endpoint)   (Source-Backed Extractions)
```

1. **Structured Component Manifest & Provenance Registry**:
   * Register discrete document blocks (main text, converted supplement tables, page coordinates, table IDs).
   * Measure **Acquisition Before/After** (e.g., publisher shell vs. verified PDF) strictly separate from **Capture/Extraction Before/After**.
2. **Atomic Observation Schema with Explicit Endpoints**:
   * Extract low-level records: `(gene, variant, subject_id, genotype, clinical_status, endpoint_type, timepoint, table_id)` before aggregation.
   * Enforce patient-ID ownership joins (**MYBPC3 21302287**, **RYR2 30403697**) to prevent conflating family members, figures, or multi-variant columns.
3. **Deterministic Source Replay Engine**:
   * Replay extraction on structured, inspected source tables rather than executing non-deterministic full-paper reruns.
   * Maintain cohort boundaries and preserve explicit `unknown/unreported` states without force-closing partitions (**RYR2 25814417**, **33315912**).

---

### 3. Meaningful Regression Tests

* **Exact Match Controls**:
  * `W4645R`: Verify abstract-derived extraction yields exact 4/2/2 counts.
  * `C2277R`: Verify patient-row extraction yields exact 8/7/1 counts.
* **Patient ID vs. Count Guard**:
  * **MYBPC3 21302287**: Table 3b must parse `278, 251...` as individual patient ID tokens, rejecting them as counts.
* **Multi-Variant & Family Disambiguation**:
  * **RYR2 30403697**: Table 1 (2 variant columns, 15 patient rows) must correctly join variants to patient rows and isolate prose-described asymptomatic parents.
* **Endpoint Isolation Test**:
  * **SCN5A 30059973** & **KCNH2 25819988**: Ensure symptom presentation (e.g., syncope, arrest) and diagnostic labels (LQTS, ECG signs) are stored under separate endpoint tags rather than merged into a single affected count.
* **Unknown Partition Preservation**:
  * **RYR2 25814417**: Ensure 97/62 (+26 unknown) preserves the 26 unknown cases without force-closing into unaffected to match gold 73/106.

---

### 4. What Must Remain Uncertain

1. **Gold Standard Inconsistencies**: The exact endpoint criteria used by gold curators (e.g., carriers vs. verified diagnosis vs. symptoms) cannot be harmonized automatically across all historical papers.
2. **Orphaned Clinical Ascertainment**: For **KCNQ1 26496715**, while the 382 Markdown lines contain patient counts, the diagnostic criteria in the missing full body remain indeterminate.
3. **Unretrieved/Unverified Sources**:
   * **SCN5A 20031634**: Timed-out PDF host leaves web search snippets unvalidated.
   * **SCN5A 25163546**: Full manuscript retrieved, but missing supplements prevent complete count verification.
4. **Engineered / Reference Variants**: Identity and ascertainment status for high-FP sets (**19632626**, **32533946**) remain unadjudicated pending manual domain review.
