# Evaluation provenance erratum

This erratum corrects only the provenance interpretation of the immutable
2026-08-21 lock. It does not alter `LOCK.json`, `predictions.json`, either
report, or any score.

The locked `report.json` says read-only recovery-layer gold scoring was
possible. That is a conservative schema-v1 fallback, not what occurred in this
run. The four generated production commands in `run_extraction.sh` each contain
`--gold-free-run`; each operator log records `Gold access disabled for this
run`; and there are no recovery-layer score directories under the production
runs. The run started from clean commit `1a5a997`, which introduced and tested
that run-wide gold boundary.

Later harness code records the gold-disabled mode explicitly in
`RUN_STATUS.json`, carries it into `predictions.json`, and lets an explicit
`read_only_layer_scoring_possible: false` override the old conservative
default. Those metadata improvements cannot be retroactively inserted into a
locked artifact.

The primary lock SHA-256 remains
`f74a86b47e8f1714b7c9e50c7e83d5fb920e705614cd9cb47b65dc2078ab8bae`.
