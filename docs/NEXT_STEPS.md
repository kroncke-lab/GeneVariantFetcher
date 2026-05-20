# What's Left to Hit 90% on Every Metric

## Current State

Current live metrics and the authoritative baseline artifact are in
`docs/CURRENT_RECALL_STATUS_2026-05-20.md`.

Do not maintain a second recall table here. This page is a short interpretation
of the current blocker classes.

## What The Long Run Proved

The pipeline is now turnkey enough to run all four requested genes from one
command, survive VPN churn, score genes with gold inputs, and produce per-gene
reports. It is not yet close to the recall target.

The recovery layers are now honest cold-start layers: ClinVar and PubTator use
PMIDs already observed in the DB by default. Gold-PMID enrichment remains only
as an explicit diagnostic mode.

KCNH2 degraded relative to the older closeout baseline because this run did not
auto-merge the KCNH2 v12 manual-recovery DB and did not use gold-PMID-conditioned
recovery. That is the right behavior for turnkey claims.

RYR2 and SCN5A are dominated by source/table coverage gaps. ClinVar helped both,
but PubTator and figure reads did not recover enough variant rows or carrier
counts to approach 90%.

## Highest ROI Blockers

1. `ELSEVIER_INSTTOKEN` is still unset. Elsevier API key alone is not enough for
   Heart Rhythm, JACC, Mayo Clinic Proceedings, Clinica Chimica Acta, and
   European Heart Journal full text.
2. Several supplemental files are present only as failed/challenge artifacts.
   The harvester now detects Cloudflare/reCAPTCHA/browser-challenge downloads,
   but the missing source data still has to be acquired.
3. KCNE1 needs a normalized per-PMID gold input before recall can be claimed.
4. SCN5A count extraction is inflated. The run reported 27,692,860 carriers in
   extraction statistics, which points to cohort/table parsing inflation that
   should be separated from missing-row recall work.
5. Oversized raw contexts can waste model time. Data Scout now prefers an
   existing deterministic `*_CLEANED.md` sibling when raw `*_FULL_CONTEXT.md`
   exceeds the configured size guard.

## Next Work

1. Get `ELSEVIER_INSTTOKEN` or manually acquire real PDFs/supplements for the
   largest KCNH2/RyR2/SCN5A missing-row PMIDs.
2. Audit top `missing_in_sqlite.csv` rows per gene from the current-status
   baseline path.
3. Fix cohort/table count semantics before trusting patient-level totals.
4. Build or import `gene_variant_fetcher_gold_standard/normalized/KCNE1_recall_input.csv`.
5. Add focused source-specific recovery for browser-challenge supplemental PDFs
   instead of retrying the same blocked LinkOut URLs.
