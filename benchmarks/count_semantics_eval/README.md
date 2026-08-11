# Count-semantics evaluation

This benchmark audits carrier-count scope separately from variant detection. It
uses locked 56-paper predictions (48 cardiac manual-gold papers plus 8 BRCA2
curator-gold papers), source-grounded adjudication, and count recall beside MAE
so clearing predictions cannot masquerade as an improvement.

The first run is in `runs/20260810_luna_xhigh_56/`. Source-level decisions that
changed the answer key are recorded in `ADJUDICATIONS_20260810.md`; original
columns remain intact and authoritative decisions live in `gold_v2_*`. The
independent blind source review is recorded in `MULTIMODEL_REVIEW_20260810.md`.
The publication-oriented study design, metric definitions, intervention steps,
negative controls, production status, and reproducibility map are consolidated
in `METHODS_20260810.md`.
Status values use the exact vocabulary in `utils/gold_standard.py`; unknown
values fail closed in both scorers and the claim-verification pilot.
