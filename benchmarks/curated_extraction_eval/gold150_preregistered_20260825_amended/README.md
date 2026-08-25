# Gold-150 preregistration — amendment 1 (2026-08-25)

This packet tree supersedes `../gold150_preregistered_20260824/`, which is
preserved unchanged as the original preregistration. The amendment was made
**before any human answer key existed**, so no key, prediction, or score
influenced it. Nothing here was selected by looking at results.

Everything else — cohort, genes, 50 papers per gene, seed `20260824`, 30
calibration / 20 holdout per gene, curation protocol, count semantics, scoring
rules — is unchanged. `EVALUATION_PLAN.md` and `SCORE_RUNBOOK.md` in the
original directory remain in force.

## Why it was amended

The original split assigned papers by ranking
`sha256("<seed>|<GENE>|<PMID>")`. The gene is inside the hash, so the same PMID
gets an independent ranking under each gene it appears in. BRCA1 and BRCA2
share papers, and six of them landed on opposite sides of the firewall:

| PMID | Original BRCA1 | Original BRCA2 |
| --- | --- | --- |
| 18627636 | holdout | calibration |
| 19949876 | holdout | calibration |
| 21232165 | holdout | calibration |
| 22713736 | holdout | calibration |
| 22970155 | holdout | calibration |
| 23683081 | calibration | holdout |

Five of the twenty BRCA1 holdout papers were BRCA2 calibration papers. Tuning
on BRCA2 calibration would therefore have burned **25% of the BRCA1 holdout**
before it was ever scored, while every per-gene manifest still looked correct
on its own. This is not a curation error; it is a property of per-gene ranking
over a cohort that shares papers.

A second, independent defect: seven of the eight cross-gene PMIDs were bound to
**different source bytes** under each gene, because each gene harvested the
paper separately.

| PMID | BRCA1 bytes | BRCA2 bytes |
| ---: | ---: | ---: |
| 18627636 | 41,257 | 41,204 |
| 19949876 | 22,924 | 53,417 |
| 22713736 | 21,307 | 20,669 |
| 22970155 | 64,172 | 68,844 |
| 23683081 | 33,388 | 33,335 |
| 25948282 | 160,960 | 29,134 |
| 26848529 | 166,434 | 165,982 |

Curators would have adjudicated different documents under one PMID, and the
frozen SHA-256 would have been a per-packet guarantee rather than a per-paper
one.

## What changed

Two deterministic tables, both in this directory and both bound by SHA-256 in
every `packet_meta.json`:

- `split_lock.csv` — one global split per cross-gene PMID. Conflicts resolve to
  **holdout**, so a shared paper is never tuned on through another gene's
  calibration. The two PMIDs whose genes already agreed (25948282 holdout,
  26848529 calibration) keep their original assignment.
- `canonical_source.csv` — one document per cross-gene PMID. The **largest**
  capture wins, which is established practice here: per-run `FULL_CONTEXT`
  differs because supplement folding mutates it, so the largest capture is the
  most complete one.

Per-gene calibration size is preserved exactly at 30/20. Pinning shifts *which*
papers move, not how many, so 12 of 150 papers changed split:

| Gene | calibration | holdout | cal→hold | hold→cal |
| --- | ---: | ---: | ---: | ---: |
| BRCA1 | 30 | 20 | 1 | 1 |
| BRCA2 | 30 | 20 | 5 | 5 |
| BMPR2 | 30 | 20 | 0 | 0 |

BMPR2 is untouched: it shares no papers with the BRCA genes.

## Verification

`scripts/audit_split_firewall.py` reports both defect classes and exits
non-zero on either. Run it before curation and again before releasing holdout:

```bash
.venv/bin/python scripts/audit_split_firewall.py benchmarks/curated_extraction_eval/gold150_preregistered_20260825_amended
```

Against this tree it reports `split firewall: OK` and `source binding: OK`.
Against the original tree it exits 1 and names all six breaches.

## Known remaining limitation

The holdout roster is tracked in git, so this cohort is **answer-key blind, not
roster-sequestered**. Anyone with the repo can read holdout PMIDs and open
their source. That is acceptable for a single scored-once evaluation with a
disciplined operator, but it is not a sealed holdout, and any claim of a
"never-seen" test set should say so plainly. Sequestering the roster with an
external custodian is the stronger form if the claim needs to survive review.
