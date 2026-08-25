# Gold-150 preregistered human evaluation

This directory freezes the evaluation design for the 150 source-screened papers
in `results/vb_curated150_20260824`: 50 each for BRCA1, BRCA2, and BMPR2.
Each gene is split by SHA-256 rank with seed `20260824` into 30 calibration and
20 holdout papers.

The tracked files are the protocol, blank answer sheet, selected-paper manifest,
source hashes, and packet metadata. `papers/` and the handoff ZIPs are deliberately
gitignored: they are local source payloads, not answer keys. The calibration
packet metadata does not reveal the holdout PMID roster.

No precision, recall, F-score, carrier coverage, or MAE is defined for this
cohort until a human curator returns an exhaustive, source-grounded answer key.
The existing extraction output must not be used to construct that key.

See [EVALUATION_PLAN.md](EVALUATION_PLAN.md) for the custody and scoring rules.
