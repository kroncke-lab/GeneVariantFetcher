# Vanderbilt Library API Access for Text/Data Mining

**Contact:** Kayce Gill, Health Sciences Librarian for Collections at Eskind
**Date:** 2026-02-03

## Overview

Vanderbilt has institutional agreements with several publishers for API access to full-text content for text/data mining purposes.

**Main resource page:** https://researchguides.library.vanderbilt.edu/textmining/resources

## Publisher Access

### Elsevier/ScienceDirect ✅ ACTIVE (insttoken issued 2026-05-21)
- **Status:** Institutional agreement active; API key + institutional token both configured.
- **Token issuance:** Elsevier Data Support (Jun Bautista) generated `X-ELS-Insttoken`
  against the `@vanderbilt.edu`-registered API key labeled `GeneVariantFetcher` on
  2026-05-21. Token is stored in `.env` as `ELSEVIER_INSTTOKEN` (file perms `600`,
  gitignored, sent only as a header by `harvesting/elsevier_api.py`).
- **Measured unlock:** 242/246 (98.4%) of previously-paywalled Elsevier candidates
  across KCNH2/SCN5A/RYR2/KCNE1/KCNQ1 now return full text. See
  `docs/RECALL_STATUS.md` ("2026-05-21 Elsevier Insttoken Activation") for details.
- **Out-of-scope hits:** `10.1016/j.cardiores.*` (Cardiovascular Research, migrated
  to Oxford in 2014) and `10.1016/j.jogc.*` (Journal of Obstetrics and Gynaecology
  Canada) are not served by the Elsevier API even with the token.
- **GVF handler:** `harvesting/elsevier_api.py`.
- **Config:** `ELSEVIER_API_KEY` + `ELSEVIER_INSTTOKEN` in `.env`.
- **Token hygiene (per Elsevier ToS):** server-side only, never in browser code,
  never in URLs/address bar, HTTPS only, may be revoked at any time without notice.

### Springer Nature ✅ AVAILABLE
- **Status:** Institutional agreement exists
- **Requirement:** Register for account with `@vanderbilt.edu` email
- **GVF handler:** Springer/BMC handler exists

### Wiley ✅ AVAILABLE
- **Status:** Institutional agreement exists
- **Requirement:** Register for account with `@vanderbilt.edu` email
- **GVF handler:** `harvesting/wiley_api.py` (ready)

### Oxford ❌ LIMITED
- **Status:** Limited to historical treaties only
- **Not useful** for biomedical literature

### Karger ⚠️ CONTACT REQUIRED
- **Status:** No institutional agreement
- **Action:** Submit TDM request via https://karger.com/pages/contact-us?mailbox=520255&mailtype=group
- **Note:** Currently blocked by Cloudflare in GVF

## Proxy Authentication

For URLs without direct API access, use Vanderbilt proxy prefix:
```
http://proxy.library.vanderbilt.edu/login?url=
```

Example: To access `https://example.com/article`, use:
```
http://proxy.library.vanderbilt.edu/login?url=https://example.com/article
```

## Action Items

### Status (2026-05-21):
- [x] Elsevier API key configured and WORKING
- [x] **Elsevier institutional token (`ELSEVIER_INSTTOKEN`) active as of 2026-05-21** — 242/246 paywalled Elsevier candidates unlocked across 5 genes
- [~] Wiley API key was working; current `WILEY_API_KEY` in `.env` is revoked — needs reissue
- [ ] Springer API key NOT configured
- [ ] Karger TDM request not submitted

### For Brett:
1. [ ] Reissue Wiley TDM API key
2. [ ] Register for Springer Nature API account
3. [ ] Submit TDM request to Karger
4. [ ] Kick off re-extraction for KCNH2/KCNQ1/RYR2/SCN5A so the unlocked Elsevier bodies actually become SQLite rows

### For Boswell (after credentials received):
1. [x] `ELSEVIER_API_KEY` + `ELSEVIER_INSTTOKEN` configured in GVF
2. [ ] Add Springer credentials to config
3. [ ] Refresh Wiley credentials when reissued
4. [ ] Implement proxy prefix fallback for remaining publishers

## Impact on GVF Recall

Pre-2026-05-21 baseline: ~70% of gold-standard papers were paywalled, with
Elsevier estimated at ~40% of that paywalled share.

**Measured impact of `ELSEVIER_INSTTOKEN` (2026-05-21):**
- 242/246 (98.4%) of previously-paywalled Elsevier candidates now return full
  text through the API tier — across KCNH2 (69/70), SCN5A (75/78), RYR2 (51/51),
  KCNE1 (22/22), and KCNQ1 (25/25).
- Bodies saved as `{PMID}_FULL_CONTEXT.md` into each run's `pmc_fulltext/`.
- PMID recall will not move on the existing DBs until a re-extraction pass is
  run; see `docs/RECALL_STATUS.md` "Next Run Plan".

Remaining paywall share (estimates, pre-re-extraction):
- Wiley: ~25% — blocked by revoked `WILEY_API_KEY`.
- Springer: ~20% — Springer key not yet configured.
- Karger/Sage/other: ~15% — Cloudflare and no-institutional-agreement issues.
