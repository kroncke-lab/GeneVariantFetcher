# gold_overrides/ — curator-supplied gold for papers the repo gold standard doesn't cover

When you add a paper to the benchmark (`registry.tsv`) whose **gold answer is not
already in** `gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv`
— e.g. a problem paper you just found, or any paper from a gene without a full
gold standard (BRCA1, APOE, MYBPC3, …) — put its expected variants here so the
benchmark can score it.

`build_fixture.py` reads gold for each paper from the **repo gold standard first**,
then falls back to `gold_overrides/<GENE>_recall_input.csv`. So you only need an
override row when the repo gold is missing the PMID. These files are committed
(they are curated variant + count data, not paper full text).

## Format

One CSV per gene, named `<GENE>_recall_input.csv`, **same schema as the repo
gold** (only the first five columns are required; leave `gold_v2_*` blank):

```csv
variant,pmid,carriers,affected,unaffected,gold_v2_carriers,gold_v2_affected,gold_v2_unaffected,gold_v2_status,gold_v2_note,gold_v2_source
R1234W,40123456,3,2,1,,,,,,
c.5678G>A,40123456,1,1,0,,,,,,
```

- `variant` — protein (`R1234W`, `p.Arg1234Trp`) or cDNA (`c.5678G>A`) notation,
  whatever the paper reports; the recall matcher normalizes both.
- `pmid` — must match the `registry.tsv` row you added.
- `carriers` / `affected` / `unaffected` — patient counts for that variant in that
  paper. Use `0` when the paper reports the variant but not that count (recall
  still credits the variant; MAE just won't have signal for it).

After editing, run `python build_fixture.py`; it will fold these rows into the
frozen `gold/` subset and stop warning about missing gold for those PMIDs.

A copy-paste starter for a new gene is in `_TEMPLATE_recall_input.csv`.
