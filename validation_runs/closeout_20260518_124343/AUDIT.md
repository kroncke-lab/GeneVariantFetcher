# Multi-Gene Closeout Audit

Run directory: `validation_runs/closeout_20260518_124343`. All scored outputs are under this directory. Publisher credential variables were unset; no Elsevier/Wiley/Springer institutional or SSO access was used.

## Progression

### KCNH2

| layer | pmids | variant_rows | unique_variants | patients | affected | unaffected |
|---|---:|---:|---:|---:|---:|---:|
| 0_current_final | 90.84% | 82.85% | 83.21% | 86.01% | 87.52% | 81.84% |

KCNH2 was scored as the pre-merged current high-recall DB per instruction; this closeout did not reconstruct layer-by-layer deltas for that database.

### KCNQ1

| layer | pmids | variant_rows | unique_variants | patients | affected | unaffected |
|---|---:|---:|---:|---:|---:|---:|
| 0_baseline | 80.00% | 27.68% | 37.71% | 45.21% | 30.24% | 68.40% |
| 1_clinvar | 94.75% | 61.20% | 61.82% | 69.84% | 64.49% | 82.21% |
| 2_pubtator | 95.41% | 61.32% | 61.82% | 69.90% | 64.58% | 82.21% |
| 3_figures | 95.41% | 62.06% | 61.98% | 70.72% | 66.05% | 82.28% |

### RYR2

| layer | pmids | variant_rows | unique_variants | patients | affected | unaffected |
|---|---:|---:|---:|---:|---:|---:|
| 0_baseline | 57.30% | 30.94% | 34.52% | 47.32% | 40.05% | 79.47% |
| 1_clinvar | 77.53% | 55.70% | 59.11% | 66.11% | 61.82% | 85.07% |
| 2_pubtator | 78.09% | 55.81% | 59.26% | 66.16% | 61.88% | 85.07% |
| 3_figures | 78.09% | 55.81% | 59.26% | 66.16% | 61.88% | 85.07% |

### SCN5A

| layer | pmids | variant_rows | unique_variants | patients | affected | unaffected |
|---|---:|---:|---:|---:|---:|---:|
| 0_baseline | 69.09% | 33.76% | 41.69% | 39.73% | 34.29% | 59.49% |
| 1_clinvar | 79.52% | 49.49% | 56.12% | 51.49% | 46.88% | 68.21% |
| 2_pubtator | 80.58% | 49.94% | 56.71% | 52.24% | 47.83% | 68.28% |
| 3_figures | 80.71% | 50.10% | 57.13% | 52.24% | 47.83% | 68.28% |

SCN5A figure layer was intentionally stopped for the wall-clock budget after 38 of 385 discovered figure-PMID reports; the score above includes partial figure DB inserts.

## 90% Target

| gene | final pmids | variant_rows | unique_variants | patients | affected | unaffected | metrics >=90% |
|---|---:|---:|---:|---:|---:|---:|---:|
| KCNH2 | 90.84% | 82.85% | 83.21% | 86.01% | 87.52% | 81.84% | 1/6 |
| KCNQ1 | 95.41% | 62.06% | 61.98% | 70.72% | 66.05% | 82.28% | 1/6 |
| RYR2 | 78.09% | 55.81% | 59.26% | 66.16% | 61.88% | 85.07% | 0/6 |
| SCN5A | 80.71% | 50.10% | 57.13% | 52.24% | 47.83% | 68.28% | 0/6 |

Only KCNH2 and KCNQ1 clear 90% on article coverage. No gene clears 90% on variant rows, unique variants, patient counts, affected counts, or unaffected counts in this closeout.

## Layer Contribution

| gene | largest variant_rows jump | delta |
|---|---|---:|
| KCNH2 | not measured in closeout; current DB includes KCNH2-specific manual/v12 recovery | n/a |
| KCNQ1 | 1_clinvar | 33.52 pp |
| RYR2 | 1_clinvar | 24.76 pp |
| SCN5A | 1_clinvar | 15.73 pp |

Cold-start layer averages across KCNQ1, RYR2, and SCN5A:

| layer | avg pmids delta | avg variant_rows delta |
|---|---:|---:|
| 1_clinvar | 15.14 pp | 24.67 pp |
| 2_pubtator | 0.76 pp | 0.23 pp |
| 3_figures | 0.04 pp | 0.30 pp |

## Top Missing Examples

### KCNH2

| PMID | variant | likely reason |
|---|---|---|
| 21216356 | G572S | abstract-only fallback; likely paywall/no full text |
| 21216356 | G628R | abstract-only fallback; likely paywall/no full text |
| 21216356 | V115M | abstract-only fallback; likely paywall/no full text |
| 21216356 | A193T | abstract-only fallback; likely paywall/no full text |
| 21216356 | G21D | abstract-only fallback; likely paywall/no full text |

### KCNQ1

| PMID | variant | likely reason |
|---|---|---|
| 34135346 | G148R | figures exist; likely image/table-only content or figure reader miss |
| 34135346 | R387Q | figures exist; likely image/table-only content or figure reader miss |
| 33632346 | S140G | abstract-only fallback; likely paywall/no full text |
| 33614747 | A344A | figures exist; likely image/table-only content or figure reader miss |
| 33382510 | R518X | figures exist; likely image/table-only content or figure reader miss |

### RYR2

| PMID | variant | likely reason |
|---|---|---|
| 19926015 | M81L | figures exist; likely image/table-only content or figure reader miss |
| 19926015 | EXON 3 DELETION | variant string is in cached text; extraction, normalization, or count matching missed it |
| 15466642 | P164S | figures exist; likely image/table-only content or figure reader miss |
| 15466642 | F4499C | figures exist; likely image/table-only content or figure reader miss |
| 15466642 | G4671R | figures exist; likely image/table-only content or figure reader miss |

### SCN5A

| PMID | variant | likely reason |
|---|---|---|
| 29514831 | Y1976A | variant string is in cached text; extraction, normalization, or count matching missed it |
| 19279983 | c.3480delT | figures exist; SCN5A figure pass was stopped after 38/385 reports |
| 17897138 | c.5040_5042delTTAinsC | figures exist; SCN5A figure pass was stopped after 38/385 reports |
| 24363796 | c.5445_5446insT | figures exist; SCN5A figure pass was stopped after 38/385 reports |
| 12106943 | D951X | figures exist; SCN5A figure pass was stopped after 38/385 reports |

## SCN5A Denylist Finding

The original `Media-Pack-2024.pdf` and `JournalCatalog*.pdf` failures did not recur as a crash condition. A new PMC administrative-form pattern did fire: `Electronic_Copyright_Form_for_Author.pdf` and `Electronic_Disclosure_Form_for_Author.pdf`. The minimal fix adds the generic filename pattern at `harvesting/supplement_scraper.py:1398-1403`, filters the PMC supplement lists in `harvesting/orchestrator.py:831-837` and `harvesting/orchestrator.py:920-925`, and covers both paths in `tests/unit/test_generic_scraper.py:267-328`.

SCN5A acquisition completed with 757 cached full-context files. Migration loaded 750 of 754 extraction JSONs; four malformed JSONs were skipped. The Azure-backed full CLI extraction hit quota/time limits, so the extraction was completed via the shared extraction step using Anthropic Haiku and then migrated normally. The incomplete piece is the full SCN5A figure-reader sweep, not text acquisition or migration.

## Assessment

The new gene-agnostic stack is useful but not enough for 90% variant recall. ClinVar is the dominant cold-start recovery layer, adding +33.52 pp variant-row recall for KCNQ1, +24.76 pp for RYR2, and +15.73 pp for SCN5A. PubTator adds small article/variant gains. Figure reading can add rows, but full sweeps are slow when there is no cache and did not materially change RYR2 or partial SCN5A recall.

Credential-gated gaps remain real for paywalled publisher full text and figures. Separate from credentials, there are pipeline gaps in slow/unbounded figure sweeps, stripped tables/supplements, malformed JSON migration, and variant/count normalization. KCNH2 remains a special case because the current DB includes KCNH2-specific manual/v12 recovery that a cold-start gene will not have.
