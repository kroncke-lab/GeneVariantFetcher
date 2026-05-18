# What's Left to Hit 90% on Every Metric

## Current state (2026-05-17 evening, after 4 parallel agents)

| Metric | Recall | Gap | Why |
|---|---:|---:|---|
| pmids | **90.8%** | ✓ over target | done |
| affected | 87.5% | 41 | needs paywall full text |
| patients | 86.0% | 107 | needs paywall full text |
| unique_variants | 83.2% | 36 | needs paywall full text |
| variant_rows | 82.8% | 71 | needs paywall full text |
| unaffected | 81.8% | 62 | needs paywall full text |

## What the agent sweep proved

After spawning 3 parallel agents on top of the manual recovery work:
- **Agent A** (paywall figure downloader) — recovered Figure 1 of PMID 29650123
  via `ars.els-cdn.com`, +5 matched variants. Other publishers (Wiley, Sage,
  Tandfonline, JACC.org, Oxford, ScienceDirect.com) all Cloudflare-block.
- **Agent B** (alternate paywall routes) — tried Crossref TDM, NIHMS
  manuscripts, OpenAIRE, BASE, Wayback Machine, DOAJ, bioRxiv against the
  top-15 paywalled PMIDs. **0 recoveries**. Snapshots that exist are
  React SPA shells, not rendered article bodies. These DOIs genuinely
  have no public OA copy.
- **Agent C** (carrier count enrichment, 2 rounds, 350+250 LLM calls) —
  +208 records total, but most landed on already-matched variants where
  another row already had counts, OR on duplicate variant entries
  (`p.Tyr43Cys` vs `p.Y43C` — same canonical form, different `variant_id`).
  Net recall delta: +5 variant_rows matches.

## Hard ceiling without credentials

The remaining gap is concentrated in **paywalled tables and figures we
cannot reach without one of:**

1. `ELSEVIER_INSTTOKEN` set in `.env` (Heart Rhythm / JACC / Mayo Clin Proc
   / Clin Chim Acta / Eur Heart J via Elsevier API). Reaches roughly
   30–40 missing variants.
2. Working Chrome SSO cookies for `scripts/fetch_paywalled.py` — the
   current macOS keychain decryption failure blocks this in headless mode.
3. Wiley TDM key with broader entitlement (current key returns 403 on
   `humu.10131` and similar). Reaches Hum Mutat / J Cardiovasc Electrophysiol.

Without these, we're at the ceiling.

## What still moves without credentials

1. **More figure downloads via the figure-reader API** — Agent A had a 1/15
   hit rate because most non-Elsevier publishers Cloudflare-block. The
   yield-per-call ratio is poor, but every additional Elsevier paper
   processed via `ars.els-cdn.com` adds ~3-5 variants on average.
2. **Better deduplication before count enrichment** — round 2 of agent C
   produced 143 records that should have moved patient/affected/unaffected
   by ~10pp. They didn't because we have duplicate `variants` rows (e.g.
   `p.Y43C` and `p.Tyr43Cys` as separate entries). The dedupe pass we
   tried merged them at the cost of one matched PMID; a more conservative
   "canonical-form lookup at insert time" path in `ingest_clinvar.py` and
   `ingest_pubtator.py` would prevent the duplicates in the first place.
3. **Cross-link individual_records counts** — the recall scorer reads
   counts from `penetrance_data` *and* aggregates from `individual_records`.
   A handful of papers have `individual_records` populated but no
   `penetrance_data` row, so the counts never surface. A small migration
   to backfill `penetrance_data` from `individual_records` aggregation
   would land another ~10-20 patient credits.

## Path to >90% on the count metrics

If we drop the strict "need full text to recover counts" rule and accept
**ClinVar carrier counts** (when the ClinVar entry is from a per-family
submission), we could backfill counts for variants we already match.
ClinVar's `ObservedIn` block has `NumberOfIndividuals` for many entries.
Worth a 1-shot agent that hits ClinVar's full XML and writes the
`(pmid, variant) → (affected, unaffected)` map.
