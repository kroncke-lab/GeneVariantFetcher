# GVF Extraction Contract (meta prompt)

**What this is.** The authoritative statement of *what GVF extracts and what it
must refuse to extract*, plus the reference prompt that encodes it. Use it three
ways:

1. As the source of truth when editing [`pipeline/prompts.py`](../pipeline/prompts.py).
   The live prompts are instantiations of this contract, not independent text: as of
   2026-07-26 both share one `_CORE_RULES` block cut down to roughly the reference
   prompt below, so a rule is stated once and cannot drift between them.
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

Two derived clinical claims are explicitly out of scope for extraction:

- **Never calculate penetrance.** `penetrance_percentage` is populated only when
  the paper explicitly states a variant-specific percentage and the exact source
  quote is retained. `affected / carriers * 100` is a downstream analysis, not
  a paper extraction.
- **Never manufacture a phenotype partition.** A diagnosis label, enrollment in
  a disease cohort, or absence of a reported diagnosis does not by itself prove
  `affected` or `unaffected`. Do not fill either field by subtraction, by setting
  `affected = carriers`, or by treating an unselected cohort as unaffected.

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
- Identity (>=1 required): cdna_notation, protein_notation, genomic_position,
  strict source-only legacy_notation, or variant_class + structural_description.
  Plus source_notation = verbatim as printed.
- `legacy_notation` is reserved for BRCA1/BRCA2 strict prefixless BIC indels
  (3-5 digits, `del|dup|ins`, then 1-20 uppercase A/C/G/T bases or a 1-3 digit
  affected-base count). It is explicitly not HGVS; retain the same bare text in
  `source_notation` and never fabricate a `c.` prefix. In other genes, the same
  shape in an explicit cDNA/nucleotide-change column is omitted-prefix cDNA, not BIC.
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
- affected is null unless the paper explicitly assigns the target phenotype to
  those carriers. “N patients/probands carried the variant” is not enough when
  the source separately reports a smaller symptom/event subset.
- penetrance_percentage is null unless the paper explicitly states a
  variant-specific percentage and supplies an exact quote. Never calculate it
  from extracted counts.
- Case report with undescribed status: count the carrier, leave affected/unaffected null.
- "Novel mutation" with no described patient = a variant, not a carrier.
- count_type in {cohort_total, screened_N, unknown}: leave the count null, note the
  aggregate.
- Provenance ambiguous: say so, emit no count.
- Same variant in several tables: SUM if cohorts are independent, take the LARGER if
  they overlap.
- A frequency column that names its own within-study denominator ("Carrier frequency
  in 7,051 cases") IS a carrier count: emit round(freq x N), and read the role from
  the noun — cases/patients/probands to affected, controls to unaffected. This is
  the one frequency that converts. A population/reference allele frequency
  (gnomAD, ExAC, TopMed, 1000 Genomes, MAF) never converts, and neither does an
  *allele* frequency stated against a cohort, whose denominator is 2N chromosomes,
  not N people. Refuse a sub-cohort qualifier ("in 500 cases with a family
  history") and a cohort that never says which side of the phenotype split it is
  ("in 300 individuals"). If the printed precision does not pin a single integer,
  abstain — emit null, never 0.

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
| Do not copy `carriers` onto `affected` | Always-on `pipeline/phenotype_count_guard.py`: clears `affected == carriers` with `unaffected in {0, None}` when N≥2 or `source_layer=figure`, unless a distinct phenotype column sourced the split. One-proband 1/1/0 text/table rows are left alone. Wired from `steps._apply_phenotype_count_guard` and figure-reader parse. Eval still scores raw, so this changes the emitted integers, not only trust. | enforced |
| Do not calculate penetrance from extracted counts | Extraction keeps only an explicitly quoted, variant-specific percentage; `DataAggregator` never derives a percentage from raw pre-trust integers. | enforced |
| Do not complete phenotype partitions arithmetically or from cohort labels | Claim verification clears ambiguous symptom-vs-diagnosis splits and never fills `affected`, `unaffected`, or their complement by subtraction. | enforced |
| Right column, right gene, right phenotype | `wrong_column`, `wrong_gene`, `phenotype_misclassified` — high-severity, source-quoted, fact+field-bound findings only. **Dormant: Steps 3.8/3.9 are parked** (see below). | **parked** |
| Animal / model-organism subjects are not human carriers | **three layers.** Tier 2/Tier 3 relevance filters in [`pipeline/filters.py`](../pipeline/filters.py) drop animal-only papers; the prompt rule reinforces that decision; then `_apply_nonhuman_clinical_count_guard` in [`pipeline/steps.py`](../pipeline/steps.py) clears human clinical count fields when the source shows a strong species signal, preserving raw under `nonhuman_source_flags`. Conservative by design — background animal-model mentions do not trip it, and an explicit human section (`# HUMAN`, "cats and humans") is an escape hatch. | enforced |
| A non-human ortholog of the target gene is not a human clinical paper | The always-on title/full-text scope gate rejects the whole paper when a species adjective directly modifies the target gene (for example, "canine BRCA2"). Explicit PMID manifests cannot bypass it. The reason persists in extraction JSON and SQLite; source replay and ClinVar/PubTator/figure recovery honor it and purge legacy evidence links. This is narrower than the general animal-subject count guard and does not reject a human-variant study merely because it uses an animal model. | enforced |
| Families ≠ individuals; probands ≠ all carriers | **classified but not acted on.** `family_count` / `proband_count` are inferred deterministically by [`pipeline/table_router.py`](../pipeline/table_router.py) and [`pipeline/extraction.py`](../pipeline/extraction.py), then written into `count_provenance`. No default consumer enforces them: `pipeline/count_classifier.py` would flag any non-`per_variant_carrier` type, but `COUNT_CLASSIFIER_POLICY` in [`config/settings.py`](../config/settings.py) defaults to `off`. The signal is on disk and unused. | **prompt-only in effect** |
| **Attribution: exclude counts credited to another publication** | `attributed_to_other_study` is in `ENFORCEABLE_REASON_CODES`, so when the reviewer runs, a high-severity source-quoted finding masks the field. But **Steps 3.8/3.9 are parked**, so today this rests on the prompt alone. [`tests/unit/test_codex_paper_eval.py`](../tests/unit/test_codex_paper_eval.py) pins the prompt guidance. | **prompt-only while parked** |
| **Row identifiers are not counts** | partial — only if a final-check finding lands as `wrong_column`. | **prompt-only in practice** |
| Identity required; reject cell artifacts | schema/notation validation at load; `validate_position` **fails open** on unregistered genes | [`variant_normalizer.py`](../utils/variant_normalizer.py) |

### Parked: the per-paper final check (Steps 3.8 / 3.9)

**Default off since 2026-07-26.** One `@xhigh` reasoning call per paper cost more
time and money than its measured effect justified, for a step that only *records*
findings. The code, prompts, reason codes, and tests are all retained and stay
green in CI, so this is a dormant switch and not a deletion.

Revive **both together**:

```bash
PAPER_FINAL_CHECK_ENABLED=1 PAPER_FINAL_CHECK_GATE_ENABLED=1 \
  gvf gvf-run <GENE> --email <you> --output ./results
```

Enabling only the composer (3.9) is a trap: with no live reviewer it can only
refuse *stale* stored findings, and stale findings raise a stage failure
demanding a "reviewer replay" that can never happen
([`cli/gvf_run.py`](../cli/gvf_run.py), `run_gvf_pipeline`) — every run would fail
acceptance forever.

What parking costs: every reason code in the table above whose enforcement point
is `paper_final_check_gate.py` goes dormant, which is why those rows read
*parked* rather than *enforced*. What it does **not** cost: field trust already
composed into existing DBs survives, because the structural gate (Step 3.7)
deliberately preserves stored `paper_final_check:` reasons
([`trust_gate.py`](../pipeline/trust_gate.py)). The structural rules — the whole
gold-free set — are unaffected and still run.

The standing design verdict (three-model review, 2026-07-20) is that the
replacement is *not* a cheaper auditor but a **deterministic coverage-diff plus
closed-loop re-extraction**: an open-loop LLM reviewer that only writes findings
cannot fix what it finds.

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

**Exception worth knowing.** The general non-human subject/count guard is a
*pre-DB clear*, not a tier:
it zeroes the count in the extraction JSON before migration, keeping the raw value
in `nonhuman_source_flags` and a summary in
`extraction_metadata.nonhuman_count_guard`. So no `penetrance_data` row carries a
`trust_tier` explaining it — the audit trail lives in the extraction JSON, not the
database. This is the same shape as the legacy `carrier_guard` / `vf`-quarantine
guards that [`TASKS.md`](../TASKS.md) wants folded into the trust record.
The target-gene-ortholog exclusion is different: it is paper-level, persists as
`extraction_metadata.paper_scope_exclusion_reason`, and blocks all downstream
evidence recovery for that PMID.

---

## Changing the contract

1. Edit this doc first — it is the spec.
2. Update `_CORE_RULES` in [`pipeline/prompts.py`](../pipeline/prompts.py) — one
   edit reaches both `EXTRACTION_PROMPT` and `COMPACT_EXTRACTION_PROMPT`. Keep it
   brace-free (it is concatenated into `str.format` templates) and keep it terse:
   instruction overhead is paid on every paper and competes with the paper for the
   window. Only the JSON output schemas are per-prompt. Compact mode triggers above
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
