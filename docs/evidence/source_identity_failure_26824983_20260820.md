# PMID 26824983 Scholar source-identity failure and repair

Date: 2026-08-20. Publication status: **HOLD**. No recovered source from this
incident was replayed into SQLite or synchronized into the shared corpus.

## What happened

Source recovery searched for PMID 26824983, the 2016 Oncotarget paper “Multiple
gene sequencing for risk assessment in patients with early-onset or familial
breast cancer” (DOI `10.18632/oncotarget.7027`). Google Scholar returned a PDF
for the unrelated 2021 Frontiers in Genetics paper “Characteristics of Germline
Non-BRCA Mutation Status of High-Risk Breast Cancer Patients in China and
Correlation with High-Risk Factors and Multigene Testing Suggestions” (DOI
`10.3389/fgene.2021.674094`). The latter cites the requested paper, so a
text-only title/DOI search could mistake citation text for source identity.

The refresh process was stopped before it changed the run DB. The original 60
variant-paper links for PMID 26824983 were unchanged. Corpus synchronization was
also stopped before the foreign text could land. The correct Oncotarget context
remains in `corpus/BRCA2/26824983/26824983_FULL_CONTEXT.md`.

The foreign source was moved, not deleted, to
`results/quarantine_wrong_source_20260820/BRCA2_26824983/`, including its source
URL record, both context copies, per-PMID fetch sidecars, and the original fetch
summary. The run-local summary was regenerated without the false Scholar
success; its other three failed-source diagnostics were preserved.

## Fail-closed repair

- Scholar candidates must carry the requested DOI in the candidate PDF URL; if
  no DOI is known, the URL must carry the requested PMID.
- The article header before Abstract/Introduction/Methods must also match the
  compact exact target title. When a DOI is known, it must match there too.
- `fetch_paywalled.py` repeats the identity check before writing/stamping a
  Scholar success; only then may it set `source_identity_verified=true`.
- Each fetch starts with an authoritative empty `summary.json`, so interrupted
  reruns cannot inherit stale completion claims.
- Outcome summarization treats an existing root summary—even an empty one—as
  authoritative. Legacy sidecars cannot override it; a Scholar success is never
  promoted without explicit source-identity verification.

## Verification

The actual quarantined Frontiers URL/text is rejected because its URL carries a
foreign DOI. Focused source-identity and acquisition-outcome tests pass 26/26.
The final repository unit suite passes 2,041 / 2,042 tests with one intentional
skip.
Grok 4.6 `xhigh` returned three successive NO-GO verdicts for concrete bypasses
(citation text beyond a fixed window, stale result/FULL_CONTEXT promotion, and
a citing paper quoting both target title and DOI), then returned **GO** after the
stable-identifier URL requirement and authoritative-empty-summary behavior were
implemented and exercised against the live functions.
