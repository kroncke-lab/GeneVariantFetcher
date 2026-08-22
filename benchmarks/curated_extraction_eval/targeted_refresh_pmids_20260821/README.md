# BRCA1/BRCA2/BMPR2 targeted refresh

These manifests refresh the 25 missing or high-risk papers in the existing
50-paper-per-gene Variant Browser cohort. The other 125 completed extraction
artifacts are copied forward unchanged by `scripts/refresh_run_db.py
--stage-extractions`.

Selection records repair/audit risk only (missing extraction, source upgrade,
cross-gene table, legacy identity, or nameless-row cleanup). Extraction itself
remains blind and receives no expected variants or gold answers.

`source_overrides.csv` binds shared PMID 19949876 to its corrected Europe PMC
JATS full context. Normal source priority would otherwise choose an older
flattened/condensed file that destroys the grouped BRCA1/BRCA2 table structure.

`BMPR2_wrong_gene_retry.txt` is the source-audit retry for PMID 42001268. Its
panel table explicitly identifies non-BMPR2 rows in a named Gene column; the
retry is expected to decrease output after the open-vocabulary row-gene guard,
so it is the one replay that intentionally disables the count-regression gate.

`BRCA1_retry.txt` records the two still-unwritten targeted papers after the
first BRCA1 refresh was interrupted. PMID 30702160's full context contained a
duplicate archive expansion that made the scanner spend unbounded time over a
42 MB source. The retry follows the idempotent supplement fold and uses
`source_overrides_BRCA1_retry.csv` to bind the paper to the complete 2.7 MB
cleaned main-plus-supplement source. The 48 already-staged papers are copied
forward unchanged.

`BMPR2_staged_audit_retry.txt` is the post-migration JSON audit follow-up. Two
legacy artifacts retained rows whose explicit source Gene cells name ACVRL1,
ENG, KCNK3, or CBLN2; one router artifact repeated the same 96 cDNA identities
twice. These are structural source defects discovered without consulting an
answer key. Their corrected outputs are expected to contain fewer rows, so this
bounded retry intentionally disables the count-regression gate.

`BRCA2_staged_audit_retry.txt` is the equivalent complete-staging follow-up for
BRCA2. Five untouched legacy JSONs retained 221 rows whose named source Gene
cells explicitly identify other panel genes; one additional paper repeated an
exact identity. The retry is selected from source structure only, remains blind
to expected variants, and intentionally permits the resulting row-count drop.

`BRCA1_staged_audit_retry.txt` records the structurally analogous BRCA1 pass:
five legacy panel-table artifacts retained explicit non-BRCA1 Gene cells, and
nine papers retained repeated exact identities (one paper is in both groups).
Markdown-decorated target cells such as `*BRCA1*` are recognized by the same
alias matcher as production and are not false audit failures. This retry also
uses no answer key and permits the expected duplicate/off-target row decrease.

`BRCA1_migration_warning_retry.txt` is the single paper whose staged source
row reached SQLite migration without a usable normalized identity. It is a
bounded replay selected from the migration warning, so it does not disclose or
depend on a curated answer.

`BMPR2_router_duplicate_retry.txt` is the one-paper follow-up after the router
learned to collapse byte-equivalent table matrices exposed under different
archive IDs. The preceding staged-audit process had already loaded the old
router before that general fix landed, so only PMID 31727138 is replayed again.

The BRCA1/BRCA2 family-count preservation retries are selected solely from
staged extraction provenance where the old claim verifier cleared a typed
`family_count` observation. The current verifier preserves that raw source
observation while the trust gate still masks it from carrier-facing fields.
These lists are not selected from gold answers or benchmark residuals.
