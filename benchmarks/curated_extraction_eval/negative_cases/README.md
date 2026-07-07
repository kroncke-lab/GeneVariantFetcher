# negative_cases/ — precision / guard cases (the other half of the benchmark)

The curated papers next door measure **recall**: did we *find* the gold variants?
They are gold-anchored and therefore **blind to garbage emission** — a run can
ace recall/MAE while minting thousands of bogus carrier rows on non-gold PMIDs,
wrong-gene variants, and annotation tables (that is exactly how the contaminated
data reached Variant_Browser staging; see `docs/RECALL_HISTORY.md` 2026-06-22).

These cases cover that gap. They assert the extraction protocol does **not** turn
annotation/population/in-silico-predictor data or wrong-gene rows into carrier
evidence. They are gold-free, self-contained (no network, LLM, or sibling DB),
and fast — run them whenever you change the extraction protocol, prompts, table
router, or count guards.

```bash
.venv/bin/python -m pytest benchmarks/curated_extraction_eval/negative_cases -q
```

## What's here

- **`annotation_table_pmid33013630.md`** — a real annotation table (`Gene |
  Variant | gnomAD allele count | SIFT | PolyPhen-2 | References`, gene-grouped
  with blank continuation cells, **no patient column**). Fed to a KCNH2 job the
  protocol must emit **0 carriers**: off-target blank-cell rows are wrong-gene,
  and the columns are annotation, not person counts.
- **`variantfeatures_residue_reference.tsv`** — KCNH2 reference residues
  (`aa_ref`) pulled from the **variantFeatures** warehouse
  (`data/variants.db`, `variant_consequences`). The local equivalent of the
  UniProt reference-AA check: `<ref><pos><alt>` assigned to a gene is a
  wrong-gene call when `<ref> != aa_ref`. Baked in so the test needs no DB.
- **`cases.tsv`** — every variant the bad extraction stamped KCNH2 from PMID
  33013630, labelled by failure mode and correct expected count:
  - `wrong_gene` (6) — residue mismatches KCNH2 (`V759I` claims V@759, KCNH2 is
    K → HCN4; etc.). Expected carriers: 0.
  - `real_gene_annot` (4) — genuine KCNH2 residue but a gnomAD/predictor row,
    not a person. Expected carriers: 0.
  - `real_gene_carrier` (2) — `R744*`, `G924A`, the only rows with explicit
    "identified in a SUDEP patient" text support. Expected carriers: 1 (SUDEP).
- **`test_negative_cases.py`** — three checks: the annotation table yields 0
  carriers (real guard, end-to-end); the variantFeatures residue slice flags
  exactly the wrong-gene rows; only the two text-supported rows expect a count.

## Adding cases

Append rows to `cases.tsv` (and the residue to
`variantfeatures_residue_reference.tsv`) for any new wrong-gene / annotation
false positive you find. To regenerate a residue from variantFeatures:

```sql
-- variantFeatures/data/variants.db
SELECT DISTINCT gene_symbol, aa_pos, aa_ref
FROM variant_consequences WHERE gene_symbol='KCNH2' AND aa_pos=759;
```

For a new annotation-table failure, drop the offending table as a `*.md` fixture
and add a `parse_routed_table(...) == []` assertion. Keep these **easy negatives
the fixed pipeline already passes**, so any future regression is a real signal.
