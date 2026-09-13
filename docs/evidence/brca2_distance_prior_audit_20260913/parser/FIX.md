# Catalogue parser fix — validation receipt

The deterministic Markdown parser now refuses to treat subject-free variant
annotation catalogues as one-person clinical rows. It also prevents an unlabeled
trailing numeric annotation from becoming a count, preserves explicit conflicting
classifications, and uses patient phenotype cells rather than variant metadata to
infer a phenotype. Explicit clinical counts and genuine patient rows remain eligible.

The **exact archived source excerpt** previously produced one row with total=1,
affected=1 and unaffected unknown. The post-fix smoke returns `[]`: zero false
catalogue observations. The original source hash was rechecked and the pre-fix
`reproduction.json` remains unchanged. Its reproducer intentionally asserts the old
bug; use the new regression tests for post-fix validation.

| Validation | Actual result |
|---|---|
| New catalogue regressions | 15 cases included in both test runs |
| Focused extraction/parser/router run | 237 passed in 21.70 s |
| Complete offline unit suite | 3,064 passed in 123.99 s |
| Ruff check | Passed |
| Ruff format check | Two files already formatted |
| Diff whitespace check | Passed |

[fix_validation.json](fix_validation.json) records the exact commands, result
summaries, current code/test/document SHA-256 values, source/excerpt and baseline
hashes, and the exact post-fix smoke code/output. Test results are transcribed from
completed tool output; no separate full-stdout log exists. The complete suite was
not rerun merely to create this receipt.

The production change is in `pipeline/extraction.py`; the 15 regression cases are
in `tests/unit/test_catalog_table_safety.py`. The protocol changelog and architecture
notes are updated. No external LLM was called and no archived count, source database,
benchmark lock or headline was changed by this fix. A source-reviewed refresh and
registered scored extraction comparison remain necessary to measure recall/MAE;
that scored arm must produce its canonical and companion figures.
