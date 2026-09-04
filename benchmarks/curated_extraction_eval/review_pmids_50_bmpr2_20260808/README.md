# BMPR2 50-paper Variant Browser review cohort

> **Frozen historical manifest (2026-08-08).** This is the seeded 50-paper cohort
> that went live for collaborator adjudication as Variant Browser pair 16
> (`bmpr2_review_50_20260808`). For collaborator review it has since been
> superseded by `../review_pmids_50_20260824_curated/`, whose
> `membership_changes.tsv` records the nine BMPR2 removals; this directory is the
> provenance record and is not an input to new runs. Rescued on 2026-09-04 from
> the unmerged PR #179 snapshot (f797260d). `scripts/prepare_review_cohort.py` is
> the selector that built it, against the `results/BMPR2/20260807_163246` run DB.


This pinned cohort was selected from `results/BMPR2/20260807_163246` for private human adjudication
in Variant Browser. It is not a gold standard.

## Contract

- All 50 papers have at least one extracted `BMPR2` variant.
- 44/50 papers (88%) have at least one
  trusted variant row with explicit carrier, affected, and unaffected counts.
- Publication years span 2001–2026.
- Source origins: abstract-only=8, browser-html-generic=2, elsevier-api=5, pmc_xml=32, publisher-free=1, wiley-api=2.
- The manifest distinguishes pipeline extraction selection, secondary LLM adjudication,
  and the human adjudication that will occur after publishing.

## Prepared staging commands

Run this only after BMPR2 variantFeatures readiness is verified and the EZproxy-backed
source gaps have been handled or explicitly accepted:

```bash
cd /Users/kronckbm/GitRepos/Variant_Browser
set -a; source .env; set +a
venv/bin/python manage.py import_features \
  --gene BMPR2 \
  --disease 'pulmonary arterial hypertension' \
  --source local \
  --db /Users/kronckbm/GitRepos/variantFeatures/data/variants.db \
  --database staging
```

That feature import creates the private gene–disease pair and its first snapshot.
Then attach the pinned GVF paper cohort:

```bash
GVF_PMID_FILE=/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/curated_extraction_eval/review_pmids_50_bmpr2_20260808/BMPR2.txt \
GVF_DATASET_LABEL=bmpr2_review_50_20260808 \
GVF_DATASET_NOTE='BMPR2/PAH initial 50-paper human-review cohort' \
bash /Users/kronckbm/GitRepos/Variant_Browser/scripts/gvf_publish.sh \
  BMPR2 /Users/kronckbm/GitRepos/GeneVariantFetcher/results/BMPR2/20260807_163246/BMPR2.db 'pulmonary arterial hypertension'
```

Both commands write only to the private Variant Browser staging/review database.
The gene–disease pair does not currently exist. If a feature import is deliberately
deferred, add `GVF_CREATE_PAIR=1` to the carrier command to create a carrier-only
snapshot; the preferred path is features first. Do not set `GVF_FULL_DB_REPLACE`
for this cohort.
