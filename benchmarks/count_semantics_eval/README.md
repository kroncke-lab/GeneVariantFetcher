# Count-semantics evaluation

This benchmark audits carrier-count scope separately from variant detection. It
uses locked predictions, source-grounded adjudication, and count recall beside MAE
so clearing predictions cannot masquerade as an improvement.

The active provenance-clean projection is
`runs/20260811_collaborator_gold_50/`: 48 cardiac manual-gold papers plus the two
BRCA2 papers with lead-approved Variant Browser adjudications. The original
56-paper run remains frozen at `runs/20260810_luna_xhigh_56/` and must be cited
as a historical mixed-provenance audit, not collaborator gold. Source-level decisions that
changed the answer key are recorded in `ADJUDICATIONS_20260810.md`; original
columns remain intact and authoritative decisions live in `gold_v2_*`. The
independent blind source review is recorded in `MULTIMODEL_REVIEW_20260810.md`.
The publication-oriented study design, metric definitions, intervention steps,
negative controls, production status, and reproducibility map are consolidated
in `METHODS_20260810.md`.
Status values use the exact vocabulary in `utils/gold_standard.py`; unknown
values fail closed in both scorers and the claim-verification pilot.
