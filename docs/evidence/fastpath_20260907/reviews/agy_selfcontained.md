### 1. Diagnosis Verdict & Biological Falsification (PMID 30059973)

**Verdict:** Your diagnosis is **scientifically and structurally correct**; the prior recommendation to force `affected_count = 185` was an overfit to a flawed gold standard.

* **Biological Falsification of Gold:** *SCN5A* exhibits pleiotropic, overlapping, and incomplete penetrance (gain-of-function LQT3 vs. loss-of-function Brugada syndrome [BrS] and progressive cardiac conduction defect [PCCD]). In PMID 30059973:
  * **Table 5** enumerates carrier ascertainment ($N=185$).
  * **Table 14** breaks down 12 variants by ECG phenotype: for E1784K ($n=69$), 29 are ECG-negative, 13 LQT3, 0 BrS, 17 PCCD, and 10 overlap phenotypes.
  * **Table 11** stratifies presentation vs. longitudinal endpoints (67.9% asymptomatic at presentation; syncope vs. major cardiac events [MCE] during follow-up).
* **The Failure of Count Copying:** Forcing `affected_count = total_carriers_observed` conflates **genotype ascertainment** with **clinical penetrance**, ignores that 44.3% are ECG-negative, and arbitrarily collapses discordant diagnostic endpoints and timepoints. The deterministic extractor correctly extracted carrier counts while abstaining from guessing an undefined clinical endpoint.

---

### 2. Concrete Code Flaws in Cohort Projection & Audit (Ranked by Severity)

```
[CRITICAL] Multi-Column Split Failure in classify_table_cohort
    │
[HIGH]     False-Positive Header Exclusion on Covariates
    │
[MEDIUM]   Unaffected Field Bias & Fixed Window in Audit
    │
[LOW]      Silent Float Tolerance in Count Validation
```

#### Finding 1: Multi-Column Split Failure on Case/Control Captions
* **Location:** [`classify_table_cohort`](file:///Users/kronckbm/.gemini/antigravity-cli/scratch/pipeline/table_cohort_phenotype.py#L40-L60)
* **Defect:** If a table caption mentions the disease (e.g., *"Table 1: SCN5A variants in Brugada syndrome and healthy controls"*), `_disease_hit(caption, run_disease)` or `case_in_caption` evaluates to true. When evaluating a clean control column header (e.g., `"Controls (n=500)"`), the check:
  ```python
  if control_column or control_caption:
      if case_in_caption or _disease_hit(caption, run_disease):
          return TableCohort(None, None, "caption_mixes_cases_and_controls", **base)
  ```
  aborts immediately. It rejects legitimate control columns before ever checking column-level disambiguation or the downstream `case_and_control_columns_present` guard.
* **Fix:** If `control_column` is explicitly matched and contains no case nouns, bypass the caption-level mixed check.

#### Finding 2: Aggressive `clinical_column_present` Rejection
* **Location:** [`classify_table_cohort`](file:///Users/kronckbm/.gemini/antigravity-cli/scratch/pipeline/table_cohort_phenotype.py#L75-L85)
* **Defect:** Iterating `header_list` and rejecting on `_HEADER_EXCLUDE_RE.search(header)` causes any standard carrier/case table containing descriptive covariate columns (e.g., `Age at onset`, `Sex`, `QTc (ms)`) to return `clinical_column_present` and abort cohort assignment completely.
* **Fix:** Header exclusions must distinguish between *ineligible count columns* vs. *benign clinical covariate columns*.

#### Finding 3: Inherent Unaffected Bias & Sliding Window in Coverage Audit
* **Location:** [`audit_table_phenotype_coverage`](file:///Users/kronckbm/.gemini/antigravity-cli/scratch/pipeline/table_phenotype_coverage.py#L45-L75)
* **Defect:** `missing = [field for field in _FIELDS[1:] if supplied[field] < len(rows)]` always marks `unaffected_count` as missing for disease-only cohort papers, perpetually driving `status` to `"review_candidates"`. Additionally, `start + 12` lines for table candidate scanning is arbitrary and easily truncates wrapped multi-line markdown headers.

#### Finding 4: Permissive Float Count Validation
* **Location:** [`audit_table_phenotype_coverage`](file:///Users/kronckbm/.gemini/antigravity-cli/scratch/pipeline/table_phenotype_coverage.py#L35-L42)
* **Defect:** `math.isfinite(value) and value.is_integer()` treats `1.0` or `0.0` as valid discrete integer counts without raising a lint/type warning if upstream emits float data types.

---

### 3. Smallest Scientifically Defensible Repo Improvement

**Do NOT build an ad-hoc phenotype ontology parser.** The smallest defensible improvement is:

1. **Keep Refusal as the Ground Truth for Complex Clinical Phenotypes:** Deterministic extraction must extract explicit carrier/case/control columns and explicitly abstain from multi-phenotype / longitudinal tables.
2. **Expose Explicit Abstention Codes:** Replace silent zeros/nulls with machine-readable reasons (`multi_phenotype_table_unresolved`, `endpoint_ambiguity`, `mixed_timepoints`).
3. **Decouple Carrier Ingestion from Phenotype Assignment:** Maintain `total_carriers_observed` as an immutable anchor.

---

### 4. Field-Sensitive Fast-Path Audit vs. Bounded Enrichment

* **Audit Assessment:** The field-sensitive audit is lightweight and safe because it introduces zero side effects, makes no LLM calls, and produces diagnostic metadata without modifying the extraction payload.
* **Bounded Enrichment Architecture:**
  * If `status == "review_candidates"`, do **not** run an unconditional full-model fallback.
  * Instead, invoke a **field-bounded prompt** that takes *only* the extracted variant identities and candidate table markdown.
  * **Hard Constraint:** The enrichment step may only populate `affected_count` / `unaffected_count` under a strictly specified target endpoint; it is prohibited from mutating or deleting existing deterministic keys (`total_carriers_observed`, `source_notation`, `chromosome`, `pos`).

---

### 5. Adversarial Test Suite

```
Test Matrix:
├── Test 1: Fallback Non-Destructiveness (LLM API Error preserves 185 carriers)
├── Test 2: Multi-Column Disambiguation (Captions with both Case & Control)
├── Test 3: Pre-Supplied Key Immutability (Preserve existing count structures)
├── Test 4: Longitudinal / Timepoint Conflict (Presentation vs. Follow-up refusal)
├── Test 5: Pleiotropic Disease Collision (BrS vs. LQT3 without target filter)
├── Test 6: Large Table Token & Line Budgeting (>500 rows, candidate caps)
└── Test 7: Duplicate Supplement Disambiguation (Table 1 vs Table S1 caching)
```

1. **Fallback Non-Destructiveness:** Patch LLM enrichment to raise `RuntimeError("LLM unavailable")`. Assert pipeline returns `status="success"`, `variants` count = 185, `total_carriers_observed` populated on all rows, and `affected_count=None`.
2. **Multi-Column Disambiguation:** Input table: Caption `"SCN5A mutations in BrS patients and healthy controls"`, Columns: `[Mutation, Cases (n=100), Controls (n=200)]`. Assert `Cases` maps to `case` cohort and `Controls` maps to `control` cohort (verifying Fix for Finding 1).
3. **Pre-Supplied Value Invariance:** Feed row with `total_carriers_observed: 12`. Assert enrichment output retains `12` regardless of candidate table signals.
4. **Timepoint / Endpoint Conflict:** Table contains `Asymptomatic at presentation (n=60)` and `Syncope during follow-up (n=9)`. Assert deterministic classifier refuses to sum them into a single `affected_count`.
5. **Pleiotropic Disease Collision:** Variant table reporting counts under separate `LQT3`, `BrS`, and `PCCD` columns with no single specified `run_disease`. Assert pipeline rejects cohort derivation with `ambiguous_disease_endpoint`.
6. **Large Table & Candidate Throttling:** 1,000-row variant table across 50 markdown tables. Assert candidate collection cuts off cleanly at `MAX_CANDIDATES = 20` without memory explosion or multi-minute regex scanning.
7. **Duplicate / Contradictory Supplement Deduplication:** Feed identical Table 1 and Table S1 where captions conflict. Assert ambiguous resolution triggers refusal rather than nondeterministic order-dependent selection.

---

### 6. Strategic Verdict: What to Ship Now

| Component | Recommendation | Justification |
| :--- | :--- | :--- |
| **Carrier Extraction** | **Ship Now** | Deterministic Table 5 extraction (185 carriers) is exact, reproducible, and cost-free. |
| **Phenotype Recovery** | **Do Not Ship (Abstain)** | Extracting `affected_count` from Tables 11/14 without an endpoint ontology violates biological reality. |
| **Refusal Telemetry** | **Ship Now** | High-value observability; clearly communicates why phenotype fields are unpopulated. |
| **Cohort Classifier Fix** | **Ship Now (with Fixes)** | Patch Finding 1 & Finding 2 before deployment to avoid dropping valid multi-column tables. |

**Most Valuable Next Step:** Fix the multi-column split bug in `classify_table_cohort`, ship the read-only coverage audit to flag candidate tables in metadata, and preserve explicit carrier counts while refusing automated phenotype count synthesis.
