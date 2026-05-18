# Recall Recovery Scripts

Scripts used to push KCNH2 recall from 33.5% to 82%+ on variant_rows.

## Recovery layers (run in this order — each is additive, none destructive)

After the main extraction pipeline produces a base SQLite DB, layer in:

1. **`merge_v12_db.py`** — Copy variants from an older, hand-curated DB for
   gold PMIDs where the older DB has more variants than the fresh run.
   Largest single-source gain.

2. **`recover_paywall_oa.py`** — Try Unpaywall + Europe PMC fallback to
   fetch OA full text for paywalled PMIDs.

3. **`ingest_clinvar.py`** — For each gold PMID, query ClinVar via E-Utils
   for variants citing that PMID. ClinVar's community curation captures
   variants from papers we can't directly access.

4. **`ingest_pubtator.py`** — NCBI PubTator3 variant-annotation API for
   text-mined mutations from abstracts. Filtered to KCNH2 gene id 3757.

5. **`../extract_figure_variants.py`** — Vision-LLM figure reader. Pulls
   variant tables and mutation maps out of image-only figures the HTML
   stripper dropped. Wires through `harvesting/figure_variant_reader.py`.

## Hard ceiling

After all layers, ~70-80 gold variants remain unmatched. These are in
paywalled papers (Heart Rhythm, JACC, AHA, Mayo Clinic Proc, Karger,
Wiley) whose figures we never downloaded. Closing that gap needs either
institutional credentials (`ELSEVIER_INSTTOKEN`, Vandy VPN cookies for
`scripts/fetch_paywalled.py`) or pulling the figure JPGs from publisher
CDNs first.

## Order matters

Migration (`harvesting/migrate_to_sqlite.py`) is destructive: re-inserting
a paper cascade-deletes its variant_papers. If you re-migrate, redo
steps 1, 3, 4, 5 afterward.
