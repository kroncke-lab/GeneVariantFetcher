# GVF Extraction Contract (meta prompt)

**What this is.** The authoritative statement of *what GVF extracts and what it
must refuse to extract*, plus the reference prompt that encodes it. Use it three
ways:

1. As the source of truth when editing [`pipeline/prompts.py`](../pipeline/prompts.py)
   — the live prompts are instantiations of this contract, not independent text.
2. As the brief to hand a reviewer, a coworker, or another model that needs to
   judge an extraction without reading the pipeline.
3. As the map from each extraction rule to the code that catches it downstream —
   see [Rule → enforcement map](#rule--enforcement-map). Rules marked
   **prompt-only** have no downstream backup: if the prompt fails, nothing else
   catches them.

This doc carries rules and enforcement points only. Live metrics live in
[`RECALL_STATUS.md`](RECALL_STATUS.md); protocol changes are logged in
[`PROTOCOL_CHANGELOG.md`](PROTOCOL_CHANGELOG.md).

---

## Standing: clinical-grade output

Output is destined for clinician-facing variant-interpretation surfaces. That
sets the optimization target: **a wrong number is worse than a missing one.**
Five invariants follow, and every extraction prompt must encode all five.

| Invariant | Meaning |
|---|---|
| **Abstain over infer** | No explicit support ⇒ `null`. Never fill a field by inference, default, or symmetry. |
| **Provenance or nothing** | A count with no source column label and no verbatim quote is not a count. Emit `null`. |
| **This study only** | Count what this paper observed. Anything credited to another publication is out of scope. |
| **Identity is load-bearing** | No variant identifier ⇒ no record. A mutation *class* is not a variant. |
| **Human clinical carriers only** | Assay replicates, animals, cell lines, and population alleles are not carriers. |

### Open blockers to the clinical-grade claim

Recorded here because the contract is only as strong as its weakest consumer.
Tracked in [`TASKS.md`](../TASKS.md).

- **Publish and reports still read raw counts.** The scorer defaults to the
  trusted projection; the publish path — the clinic-facing one — does not. A
  quarantined count can still reach a display surface.
- **`validate_position` fails open** for genes with unknown length
  ([`utils/variant_normalizer.py`](../utils/variant_normalizer.py) returns `True`),
  so position validation silently does nothing on unregistered genes.
- **Headline metrics cover four cardiac genes only** (KCNH2, KCNQ1, SCN5A, RYR2).
  Anything served outside them is unmeasured output; per-gene-class calibration
  (`scripts/precision_sample.py`) is the gate.
- [`README.md`](../README.md) still declares research-use-only. That claim and
  this standing are in conflict; resolving it is a deliberate decision, not a
  side effect of an edit.

---

## Reference prompt

```text
ROLE
You extract clinical-grade genetic evidence from biomedical literature. Output feeds
clinician-facing variant-interpretation tools. A wrong number is worse than no number:
when support is not explicit, abstain (null). Never infer to fill a field.

TASK
For the target gene ONLY, extract every variant this study itself observed, with
per-variant carrier evidence and exact provenance.

PER VARIANT
- Identity (>=1 required): cdna_notation, protein_notation, genomic_position, or
  variant_class + structural_description. Plus source_notation = verbatim as printed.
- Counts: total_carriers_observed, affected, unaffected, uncertain.
- Individuals: id, age at onset/diagnosis/evaluation, sex, affected status, phenotype,
  ancestry, geographic origin, verbatim evidence sentence.
- Context: clinical significance, functional assays, segregation, population frequency.
- Provenance for EVERY count: verbatim column header + count_type (per_variant_carrier |
  family_count | proband_count | cohort_total | screened_N | case | control |
  unaffected_control | unknown) + table/row/column/section + short exact quote.
- Study level: study_type, study_design, ascertainment, cohort_source, population.

BEFORE READING ANY TABLE
Read caption, headers, footnotes, symbol legends, and the prose introducing it. Decide
(a) whose observations it reports — this study, a compilation of others, or both;
(b) what each column counts — families, individuals, alleles, probands, cases.
Count only what THIS study observed. Exclude anything credited to another publication.

NEVER COUNT AS CARRIERS
- Assay replicates ("n cells", current density, MPRA/saturation mutagenesis) ->
  functional study, counts null.
- Animal / model-organism / cell-line subjects.
- Denominators: "Total cases" when a "Carriers" column exists; screened N; MAF;
  allele counts (gnomAD/ExAC/TOPMed/1000G); "no. of occurrences".
- Row identifiers: patient ID, case no., adult/child number.
- Bare het/hom/genotype columns.
- Families reported as individuals.
- Study-wide totals pushed onto individual variants — likewise phenotype,
  demographics, penetrance %. Aggregates go in notes only.

ABSTAIN (do not guess)
- unaffected is null unless the paper explicitly describes carriers without disease.
  Never derive unaffected = 0.
- Case report with undescribed status: count the carrier, leave affected/unaffected null.
- "Novel mutation" with no described patient = a variant, not a carrier.
- count_type in {cohort_total, screened_N, unknown}: leave the count null, note the
  aggregate.
- Provenance ambiguous: say so, emit no count.
- Same variant in several tables: SUM if cohorts are independent, take the LARGER if
  they overlap.

REJECT AS NON-VARIANTS
Class-only cohorts ("34 patients with pore-region mutations"); single nucleotide or
amino-acid letters; gene symbols; nan/NA; significance strings in a notation field.

SELF-CHECK BEFORE RETURNING
- affected + unaffected + uncertain <= carriers; no negatives.
- Any count >10x the paper's median carrier count: denominator or allele count?
- 50+ variants with similar counts and ~100% penetrance: you read an assay table as
  patients. Re-read.
- Every populated count has a column label and a quote. If not, null it.
```

---

## Rule → enforcement map

The prompt is the first line of defense; the trust gate is the second. This table
says which is which. **Prompt-only** rows have no automated backup.

| Prompt rule | Downstream check | Where |
|---|---|---|
| Parts must fit the total; no negatives | `arith_inconsistent`, `negative_count` | [`trust_gate.py`](../pipeline/trust_gate.py) + final-check enforceable |
| Denominators are not carriers | `count_is_total` (`count_type` in cohort_total / screened_n / study_total / total) | `trust_gate.py` |
| Population alleles are not carriers | `population_count` (gnomAD/MAF/AC label regex; carriers > 100,000; > 50 in population/biobank/GWAS designs) | `trust_gate.py` |
| No study-wide totals per variant | `paper_outlier` (> 10× paper median, ≥ 50 absolute) — catches only the *large* cases | `trust_gate.py` |
| No clinical counts from reviews / pure functional / GWAS | `study_type_mismatch` | `trust_gate.py` |
| Never derive `unaffected = 0` | `implied_unaffected_zero` — masks **only** the `unaffected` field, in designs that enroll unaffected carriers (population, biobank, case-control, family segregation). Dormant on proband case reports and unknown design. | `trust_gate.py` |
| Right column, right gene, right phenotype | `wrong_column`, `wrong_gene`, `phenotype_misclassified` — high-severity, source-quoted, fact+field-bound findings only | [`paper_final_check_gate.py`](../pipeline/paper_final_check_gate.py) |
| Animal / model-organism subjects are not human carriers | **three layers.** Tier 2/Tier 3 relevance filters drop animal-only papers ([`filters.py:188`](../pipeline/filters.py:188), [`:475`](../pipeline/filters.py:475)); the prompt rule; then a deterministic post-extraction guard, `_apply_nonhuman_clinical_count_guard` ([`steps.py:322`](../pipeline/steps.py:322), wired at `:1946` and `:2012`), clears human clinical count fields when the source shows a strong species signal, preserving raw under `nonhuman_source_flags`. Conservative by design — background animal-model mentions do not trip it, and an explicit human section (`# HUMAN`, "cats and humans") is an escape hatch. | enforced |
| Families ≠ individuals; probands ≠ all carriers | **classified but not acted on.** `family_count` / `proband_count` *are* inferred deterministically from column headers ([`table_router.py:1532`](../pipeline/table_router.py:1532), [`extraction.py:2003`](../pipeline/extraction.py:2003)) and written into `count_provenance`. But no default consumer enforces them: neither gate references them, and `pipeline/count_classifier.py` — which would flag any non-`per_variant_carrier` type — has `COUNT_CLASSIFIER_POLICY` defaulting to `off` ([`settings.py:321`](../config/settings.py:321)). The signal is on disk and unused. | **prompt-only in effect** |
| Attribution: exclude counts credited to another publication | `attributed_to_other_study` — emitted when the source explicitly credits the row to a different publication (reference/citation column, footnote marker, "previously reported (ref 12)", compilation-table caption), and **enforceable**, so a high-severity source-quoted finding masks the field. Distinct from `unsupported_count` (thin or unlocated support), which stays advisory by design. A pinning test ([`test_codex_paper_eval.py:533`](../tests/unit/test_codex_paper_eval.py:533)) separately holds the prompt guidance in place. | enforced |
| **Row identifiers are not counts** | partial — only if a final-check finding lands as `wrong_column`. | **prompt-only in practice** |
| Identity required; reject cell artifacts | schema/notation validation at load; `validate_position` **fails open** on unregistered genes | [`variant_normalizer.py`](../utils/variant_normalizer.py) |

### How quarantine behaves

Soft, per-fact, reversible. The gate writes `trust_tier`, `trust_reasons`, and
`trust_rule_version` (a content hash of the rule set) onto the
`penetrance_data` row. It **never deletes the row and never NULLs the stored
value**, so a rule change re-tiers existing facts instead of destroying
evidence. Masking is per-field where a reason allows it.

Consumers read a projection, not a mutation:
`cli/compare_variants.py --trust-tier trusted` (the default) replaces quarantined
count fields with `NULL` at query time while keeping the identity row — so PMID
and unique-variant recall are unaffected and only the number is withheld.
`--trust-tier all` is the raw diagnostic mode.

**Exception worth knowing.** The non-human guard is a *pre-DB clear*, not a tier:
it zeroes the count in the extraction JSON before migration, keeping the raw value
in `nonhuman_source_flags` and a summary in
`extraction_metadata.nonhuman_count_guard`. So no `penetrance_data` row carries a
`trust_tier` explaining it — the audit trail lives in the extraction JSON, not the
database. This is the same shape as the legacy `carrier_guard` / `vf`-quarantine
guards that [`TASKS.md`](../TASKS.md) wants folded into the trust record.

---

## Changing the contract

1. Edit this doc first — it is the spec.
2. Update [`pipeline/prompts.py`](../pipeline/prompts.py). Both
   `EXTRACTION_PROMPT` and `COMPACT_EXTRACTION_PROMPT` must stay in sync; shared
   text belongs in `TABLE_ATTRIBUTION_GUIDANCE`. Compact mode triggers above
   `HIGH_VARIANT_THRESHOLD` (30) variants.
3. If a new rule deserves automated backup, add a reason code to
   `trust_gate.RULE_IDS` (and to `REASON_FIELDS` if it should mask one field
   rather than the whole fact). `rule_version()` changes automatically.
4. Regression-check on [`benchmarks/curated_extraction_eval/`](../benchmarks/curated_extraction_eval/README.md)
   — seconds, and free to score.
5. Confirm any headline claim with `scripts/run_recall_suite.py` before writing a
   number anywhere. Six recall axes are scored: PMIDs, variant rows, unique
   variants, patients (carriers), affected, unaffected — plus MAE/RMSE on
   carriers/affected/unaffected and the gold-PMID precision proxy.
6. Append a row to [`PROTOCOL_CHANGELOG.md`](PROTOCOL_CHANGELOG.md).
