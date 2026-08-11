# Fifty-paper-per-gene review cohort (historical initial cohort)

> **BRCA2 was superseded on 2026-08-11 by
> [`../review_pmids_20260811_brca2_provenance/`](../review_pmids_20260811_brca2_provenance/).**
> Its active queue has 46 papers after excluding six internally derived
> benchmark papers; four were removed from the live 50-paper snapshot and two
> were already absent. The other seven gene files remain the current 50-paper
> operational scopes.

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
| BRCA2 | 50 historical; 46 active | two active benchmark papers have Variant Browser collaborator gold; all BRCA2 remains **excluded** from headline metrics |
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
import pathlib
import sqlite3

from utils.pmid_utils import is_valid_pmid

DBS = {  # canonical run DBs, kept in sync with run_benchmark.py CANONICAL_DBS
 "APOE":"validation_runs/canonical_baseline/APOE.db",
 "BRCA1":"validation_runs/canonical_baseline/BRCA1.db",
 "BRCA2":"validation_runs/canonical_baseline/BRCA2.db",
 "KCNH2":"validation_runs/canonical_baseline/KCNH2.db",
 "KCNQ1":"validation_runs/canonical_baseline/KCNQ1.db",
 "MYBPC3":"validation_runs/canonical_baseline/MYBPC3.db",
 "RYR2":"validation_runs/canonical_baseline/RYR2.db",
 "SCN5A":"validation_runs/canonical_baseline/SCN5A.db",
}
out=pathlib.Path("benchmarks/curated_extraction_eval/review_pmids_50"); out.mkdir(exist_ok=True)
for g,p in DBS.items():
    c=sqlite3.connect(p)
    rows=c.execute("SELECT pmid, COUNT(*) n FROM penetrance_data WHERE pmid IS NOT NULL AND TRIM(pmid)!='' AND total_carriers_observed IS NOT NULL AND total_carriers_observed>0 GROUP BY pmid ORDER BY n DESC, pmid ASC LIMIT 50").fetchall()
    c.close()
    pmids=[str(pmid).strip() for pmid, _ in rows]
    invalid=sorted({pmid for pmid in pmids if not is_valid_pmid(pmid)})
    if invalid:
        raise ValueError(
            f"{g}: selected rows contain non-PMID identifiers: {invalid}; "
            "refresh it from canonical extraction filenames before regenerating"
        )
    if len(pmids) != 50 or len(set(pmids)) != 50:
        raise ValueError(f"{g}: expected 50 unique PMIDs, found {len(set(pmids))}")
    (out/f"{g}.txt").write_text("\n".join(pmids)+"\n")
PY
```

Keep these files stable for reproducibility; put any materially different future
cohort in a new versioned directory. The fail-fast validation is intentional:
repair a malformed identifier in the canonical extraction/DB rather than
filtering it out or aliasing it only in the review cohort. For example, an older
KCNH2 snapshot stored PMCID `PMC9522753` as its PMID; its canonical extraction
filename identifies the paper as PMID `34546463`, while migration preserves
`PMC9522753` in `papers.pmc_id`.
