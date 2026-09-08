I'll treat this as a bounded adversarial review: read the named extraction and table-cohort paths first, then check the PMID 30059973 source tables against the claimed fast-path diagnosis.The short-circuit returns success before the model, but table-cohort still runs on that result. Next I’ll check whether 30059973 is refused as a carrier table and whether Tables 11/14 actually contradict a carrier-to-affected copy.I have Tables 11 and 14. I still need Table 5 plus `classify`/`derive` to test whether the lane would copy carriers onto affected or correctly refuse.**Verdict:** Your diagnosis of the *source* is right. The prior “identity-only” label is wrong. “Skips phenotype tables” is right, and that is the real protocol hole. Copying Table 5 `n` onto `affected` would match gold and falsify the paper.

## 1. Falsification

**Prior reviewers, half-wrong.** `_allow_deterministic_table_short_circuit` (extraction.py ~1256) only refuses when the census saw count columns and parsed rows have *no* `total_carriers_observed` / `affected` / `unaffected` / `patients.count`. PMID 26496715 is the type specimen. 30059973 is not that paper. Table 5 is a 185-row mutation roster (`n` = occurrence). Source has `c5350G>A pGlu1784Lys … 69`. The fast path is **carrier-complete, phenotype-blind**. TASKS.md’s “identity-only fixed-width fast path (30059973)” should be retired as a diagnosis.

**You, mostly right; one correction.** The fixed-width return (~8717–8773) does skip the router/LLM, so Tables 11/14 are never parsed. It does **not** skip the common success boundary (~8414–8458): patient-row then `derive_table_cohort_phenotype_counts` still run. They only see Table 5 rows. For this caption (`SCN5A mutations (442 patients, 445 mutations, 185 unique mutations)`) plus bare `n`, classify currently returns `no_cohort_evidence` (no builtin disease word; bare `n` is not a case-noun column). The copy is **not happening now**. Gold-fit is a trap, not a missing stamp.

**Source (checked):** Table 14 E1784K: 29+13+0+17+0+0+10 = **69**. Table 11 presence: 0 arrest + 9 syncope + 60 asymptomatic = **69**; 10 MCE in follow-up. Paper: 196/442 (44.3%) negative ECG; 67.9% asymptomatic at diagnosis (Figure 1). Those are different endpoints, times, and diseases. `affected = 69` is “ascertained mutation carriers,” not phenotype-positive.

## 2. Ranked findings

**P0 — Fast-path gate is the wrong predicate.** Any carrier `n` blesses the extract. Phenotype-bearing tables can sit unread. Unconditional full-model fallback on a 185-row PDF is the expensive wrong fix.

**P0 — Do not ship `affected = n` for this shape.** It would raise exact-match vs gold and violate EXTRACTION_CONTRACT / TASKS.md (“carrier N is not automatically affected N”).

**P1 — Table-cohort mutation-roster hole (latent).** TIER_COLUMN fires on a case-noun count column (`No. of patients`) with **no disease requirement**. Caption excludes `carriers?` but not “unique mutations” / functional-effect catalogs. Table 5 is saved by bare `n` + no disease token. If the parser’s `carriers_column_label` becomes a patient-noun, 185 mixed-phenotype frequencies become `affected`. `stampable_case_copy` would then **certify** a parser-emitted `affected == n`.

**P1 — Observability lies.** `table_cohort_phenotype_derivation` reports `attempted=True` on parsed rows only. Census/other-table coverage is invisible. Reviewers then call the paper “identity-only.”

**P2 — No affected integer is scientifically determined here without an endpoint policy.** ECG-positive (40), isolated LQT3 (13), isolated BrS (0), syncope-at-dx (9), MCE-in-FU (10) are all source-true.

**P3 — `_BUILTIN_DISEASE_RE` `syndrome`/`cancer`, fail-open settings `except Exception`, control-arm `affected=0` vs case-arm NULL unaffected.** Documented asymmetry; dangerous only if Table 5 is misclassified as `case`.

## 3. Smallest defensible repo change

Not a PMID patch. Not “disable short-circuit.” Not “project Table 5.”

**Field-sensitive coverage audit on the short-circuit, then bounded enrichment:**

1. Census/header/caption scan for phenotype-bearing tables (reuse `_HEADER_EXCLUDE_RE` / clinical-caption language already in `table_cohort_phenotype.py`).
2. If parsed rows have carriers but not those fields → **do not** return the large-table success object as complete.
3. Keep Table 5 identities and `n`. Run table-router / deterministic parse **only** on the unread clinical tables (here 11 and 14). Attach splits as structured facts, not as `affected`.
4. If enrichment fails or the model times out, **keep** the deterministic 185/69. Record `phenotype_tables_unparsed` / `enrichment_failed` in metadata.
5. Tighten classify: mutation-roster / “unique mutations” / functional-effect captions refuse `case` even when the count column says patients. Do not stamp `affected == n` without an explicit phenotype column or a true single-disease case-series caption.

That is the TASKS.md item already queued (~307–312), stated as a gate rather than a paper special case.

## 4. Audit vs full-model fallback

A field-sensitive audit is the right architecture. Unconditional full-text on this PDF is a cost and truncation risk and still invites the model to copy 69. Enrichment is justified only for tables the census marked phenotype-bearing and the parser did not consume. Target-disease, endpoint, and timepoint stay in the payload; they are not collapsed.

## 5. Adversarial tests (offline; no gold chasing)

| Case | Pass |
|---|---|
| Model timeout after enrichment | Table 5 `n=69` and 185 identities remain; A/U stay NULL |
| Mixed cohort (5+11+14) | `affected` not set to 69 or 185 |
| Already-supplied A/U from an explicit column | table-cohort does not overwrite |
| `--disease Brugada` | isolated LQT3 13 and PCCD 17 are not BrS-affected |
| Endpoint/timepoint | presentation (60 asymptomatic) ≠ FU MCE (10) ≠ baseline ECG (29 negative) |
| Huge table/cost | 185-row Table 5 never forces full-text LLM; only 11/14 (or equivalent) enter enrichment |
| Duplicate representations | 69 (T5) vs T14 partition vs T11 presence; one identity, no triple count |
| Gold-fit trap | exact gold `69/69/0` is **failure** if source shows 29 ECG-negative |
| Census | count column **and** phenotype tables → short-circuit must not fire on carriers alone |
| Bare `n` + “442 patients” caption | classify refuses `case` |

## 6. Can count recovery ship without a phenotype ontology?

**Not as `affected`/`unaffected` integers.** Table 14 is already a closed partition of the same 69 people; you can persist those columns as labelled counts (endpoint, timepoint, disease name, table id). Mapping them onto the two gold fields **is** an ontology (or a silent gold convention).

**Better justified now:** (a) document refusal for mixed-phenotype carrier catalogs; (b) coverage metadata so a fast path cannot look complete; (c) mutation-roster refusal in classify; (d) preserve deterministic carriers through enrichment failure. That is protocol, not a new disease map. Do not enable `COUNT_RECOVERY` or paper-ascertainment to “fix” 30059973.

**Gold note:** if gold stores affected = carriers for this PMID, that is curator policy. Extraction should not be optimized toward it (TASKS.md already says so).
