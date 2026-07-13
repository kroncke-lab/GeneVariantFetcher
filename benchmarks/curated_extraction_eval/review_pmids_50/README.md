# Fifty-paper-per-gene review cohort (current)

These PMID lists are the **initial review cohort: 50 papers per gene** across all
eight gene-disease pairs (**400 papers total**) — the operational set published
for human adjudication / re-review. They are **not** a gold standard.

## Selection (deterministic, per gene)

For each gene, the 50 papers with the most **penetrance-data rows carrying a real
carrier count** (`total_carriers_observed > 0`) in that gene's canonical run DB —
i.e. the papers with the most variant×carrier facts for a human to adjudicate.
Ranking on carrier-bearing penetrance rows (rather than raw variant count)
deliberately excludes ClinVar/PubTator aggregation "papers" (which carry many
variants but no primary carrier counts and make poor review targets).

Tie-break is ascending PMID, so the lists are reproducible.

| Gene | Papers | Gold standard for metrics |
|---|---:|---|
| KCNH2 | 50 | **human-curated** (counts toward recall/MAE/precision) |
| KCNQ1 | 50 | **human-curated** (counts toward recall/MAE/precision) |
| SCN5A | 50 | **human-curated** (counts toward recall/MAE/precision) |
| RYR2 | 50 | **human-curated** (counts toward recall/MAE/precision) |
| APOE | 50 | `gold_overrides` (curator/derived — **excluded** from headline metrics) |
| BRCA1 | 50 | `gold_overrides` (curator/derived — **excluded** from headline metrics) |
| BRCA2 | 50 | `gold_overrides` (curator/derived — **excluded** from headline metrics) |
| MYBPC3 | 50 | `gold_overrides` (curator/derived — **excluded** from headline metrics) |

Only the four cardiac genes have fully human-curated gold standards, so
recall / precision / MAE are reported **against those four only**. The other
four are review targets but not measured (see `../gold_overrides/README.md` and
`docs/RECALL_STATUS.md`).

## Relationship to earlier cohorts

Supersedes `../review_pmids_12/` (the frozen July 9, 2026 experiment cohort, 12
per gene). Per this benchmark's convention, a materially different cohort lives
in its own versioned directory rather than mutating the old one.

## Regenerate

```bash
python - <<'PY'
import sqlite3, pathlib
DBS = {  # canonical run DBs, kept in sync with run_benchmark.py CANONICAL_DBS
 "APOE":"results/APOE/20260621_072155_full_redo/APOE.db",
 "BRCA1":"results/BRCA1/20260616_132646/BRCA1.db",
 "BRCA2":"results/BRCA2/20260606_134517_hereditary_breast_cancer_500/BRCA2/20260606_134519/BRCA2.refresh_20260606_205358.db",
 "KCNH2":"results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db",
 "KCNQ1":"validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db",
 "MYBPC3":"results/MYBPC3/20260616_132646/MYBPC3.refresh_20260617_091043.db",
 "RYR2":"validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db",
 "SCN5A":"validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db",
}
out=pathlib.Path("benchmarks/curated_extraction_eval/review_pmids_50"); out.mkdir(exist_ok=True)
for g,p in DBS.items():
    c=sqlite3.connect(p)
    rows=c.execute("SELECT pmid, COUNT(*) n FROM penetrance_data WHERE pmid IS NOT NULL AND TRIM(pmid)!='' AND total_carriers_observed IS NOT NULL AND total_carriers_observed>0 GROUP BY pmid ORDER BY n DESC, pmid ASC LIMIT 50").fetchall()
    (out/f"{g}.txt").write_text("\n".join(str(r[0]) for r in rows)+"\n"); c.close()
PY
```

Keep these files stable for reproducibility; put any materially different future
cohort in a new versioned directory.
