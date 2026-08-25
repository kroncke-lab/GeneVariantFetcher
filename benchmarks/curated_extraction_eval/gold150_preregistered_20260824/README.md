> **SUPERSEDED 2026-08-25 — do not curate from this tree.**
>
> This split has a cross-gene firewall breach: six PMIDs are calibration for one
> gene and holdout for another, so five of the twenty BRCA1 holdout papers are
> BRCA2 calibration papers. Seven of the eight cross-gene PMIDs are also bound to
> different source bytes per gene. Both defects, the cause, and the fix are
> documented in `../gold150_preregistered_20260825_amended/README.md`.
>
> This directory is preserved unchanged as the original preregistration. Use
> `../gold150_preregistered_20260825_amended/` for all curation and scoring.
> `EVALUATION_PLAN.md` and `SCORE_RUNBOOK.md` here remain in force for both.

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
