# Coverage Improvement Plan

This plan replaces the early 50-PMID prototype plan. Use the 2026-05-15 KCNH2
manual-recovery score as the current baseline.

## Baseline

| Metric | Current | Target |
|---|---:|---:|
| Unique variants | 323/530 (60.9%) | 90% |
| Variant rows | 542/991 (54.7%) | 90% |
| PMIDs | 184/262 (70.2%) | 90% |
| Patients/carriers | 1758/2674 (65.7%) | 90% |

Current variant-row gap to 90%: 350 rows. Current missing-in-SQLite file has
449 rows across 107 PMIDs. The top 10 PMIDs account for 227 missing rows.

## Highest-Value PMIDs

| PMID | Missing rows | Action |
|---|---:|---|
| 15840476 | 86 | Unlock Elsevier/Heart Rhythm with `ELSEVIER_INSTTOKEN`, VPN/IP access, or valid manual PDF |
| 14661677 | 24 | Re-fetch/re-extract; current DB has only 2 variant links |
| 29650123 | 21 | Re-extract full-ish context; current DB has one matched variant |
| 24667783 | 20 | Recover supplements/table detail; current extraction has placeholder positions |
| 16922724 | 20 | Wiley access or valid manual PDF |
| 23098067 | 16 | Open access body exists; improve table/supplement extraction |
| 23631430 | 12 | Sage/Liebert access or valid manual PDF |
| 17905336 | 10 | Re-extract full-ish context |
| 12402336 | 9 | Wiley Human Mutation access or valid manual PDF |
| 27871843 | 9 | Re-fetch/re-extract |

## Plan

1. Run paywall recovery from an institutional network or Vanderbilt VPN.
2. Set `ELSEVIER_INSTTOKEN` and retry PMID 15840476 first.
3. Re-extract the highest-loss PMIDs that already have usable context:
   29650123, 24667783, 23098067, and 17905336.
4. Audit `count_mismatches=117`; separate cohort over-counting from true recall
   misses.
5. Re-score after every recovery batch with `scripts/run_recall_suite.py`.

## Commands

```bash
.venv/bin/python scripts/fetch_paywalled.py \
  --input results/KCNH2/20260506_102238/pmc_fulltext/paywalled_missing.csv \
  --output results/KCNH2/20260506_102238/pmc_fulltext \
  --no-headless
```

```bash
.venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
  --db KCNH2=results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db \
  --outdir recall_metrics/kcnh2_after_recovery
```

## Do Not Revive

The older plan targeted a 50-PMID prototype set and proposed new ad-hoc
normalizers and browser-fetch rewrites. Those are obsolete. Use the current
paths:

- `scripts/fetch_paywalled.py` for recovery.
- `harvesting/browser_html/content_quality.py` for stub rejection.
- `pipeline/table_router.py` for clinical table routing.
- `cli/compare_variants.py` for comparison and matching.
