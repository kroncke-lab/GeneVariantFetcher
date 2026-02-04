# Vanderbilt Library API Access for Text/Data Mining

**Contact:** Kayce Gill, Health Sciences Librarian for Collections at Eskind
**Date:** 2026-02-03

## Overview

Vanderbilt has institutional agreements with several publishers for API access to full-text content for text/data mining purposes.

**Main resource page:** https://researchguides.library.vanderbilt.edu/textmining/resources

## Publisher Access

### Elsevier/ScienceDirect ✅ AVAILABLE
- **Status:** Institutional agreement exists
- **Requirement:** Register with `@vanderbilt.edu` email to get institutional token
- **Registration:** See text/data mining guide link above
- **GVF handler:** `harvesting/elsevier_api.py` (ready, needs API key)
- **Config:** Set `ELSEVIER_API_KEY` environment variable

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

### Status (2026-02-03):
- [x] Elsevier API key configured and WORKING
- [x] Wiley API key configured and WORKING
- [ ] Springer API key NOT configured
- [ ] Karger TDM request not submitted

### For Brett:
1. [ ] Register for Springer Nature API account
2. [ ] Submit TDM request to Karger

### For Boswell (after credentials received):
1. [ ] Configure `ELSEVIER_API_KEY` in GVF
2. [ ] Add Springer credentials to config
3. [ ] Add Wiley credentials to config
4. [ ] Implement proxy prefix fallback for remaining publishers

## Impact on GVF Recall

Current bottleneck: ~70% of Excel baseline papers are paywalled

Estimated breakdown:
- Elsevier (Circulation, Heart Rhythm, JACC): ~40%
- Wiley: ~25%
- Springer: ~20%
- Others: ~15%

**With institutional API access, we could potentially unlock 85%+ of paywalled content.**
