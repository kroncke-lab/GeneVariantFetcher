# Vanderbilt Institutional Access - Implementation Guide

**Date:** 2026-02-07  
**Created for:** GVF (GeneVariantFetcher) Paper Acquisition Pipeline

## Current Status: Configuring Institutional Access

### Publisher Access Matrix

| Publisher | API Status | Proxy Required | Email Registration | Notes |
|-----------|-------------|----------------|-------------------|--------|
| **Elsevier/ScienceDirect** | ✅ **VALID** | No | brett.kroncke@gmail.com | Key: 1870c1053224933f7b3237ca98abe9ef |
| **Wiley** | ✅ **VALID** | No | brett.kroncke@gmail.com | Key: a4f6e330-7027-4e71-a467-8df766f935c7 |
| **Springer/Nature** | ⚠️ **NEED KEY** | No | @vanderbilt.edu required | Requires institutional registration |
| **Oxford Journals** | ❌ **LIMITED** | Yes | N/A | Historical treaties only, not biomedical |
| **Karger** | ❌ **NO AGREEMENT** | Yes | N/A | Manual request required |

## Installation & Configuration Steps

### 1. Springer API Registration (Runtime Config)
**Important:** Springer requires institutional domain registration

```bash
# Register at: https://dev.springernature.com/
- Use institutional email: brett.kroncke@vanderbilt.edu
- Select "I am affiliated with an institution"
- Choose "Vanderbilt University" from dropdown
- API key will be emailed within 24 hours
```

### 2. Proxy Configuration for Non-API Publishers

For publishers without institutional API agreements (Oxford, Karger, etc.):

```python
# Wilson Library proxy prefix
VANDERBILT_PROXY_PREFIX = "http://proxy.library.vanderbilt.edu/login?url="

# Usage example:
original_url = "https://journals.elsevier.com/heart-rhythm/fulltext/S1547512706041299"
proxy_url = f"{VANDERBILT_PROXY_PREFIX}{original_url}"
```

### 3. Environment Configuration Update

Add to `.env` (already partially configured):
```bash
# Already existing
ELSEVIER_API_KEY=1870c1053224933f7b3237ca98abe9ef
WILEY_API_KEY=a4f6e330-7027-4e71-a467-8df766f935c7
SPRINGER_API_KEY=<ADD_AFTER_REGISTRATION>

# Proxy configuration
VANDERBILT_PROXY_ENABLED=true
VANDERBILT_PROXY_PREFIX=http://proxy.library.vanderbilt.edu/login?url=
GATEWAY_INTERFACE=library_access
```

## 224-Paper Queue Implementation

### Testing Protocol (Phase 1: 10 papers)
1. **Test Elsevier API** - 4 papers (Circulation, Heart Rhythm)
2. **Test Wiley API** - 3 papers (AJHG, Ann Hum Genet)
3. **Test Springer API** - 3 papers (J Med Genet, Nat Med)

### Scaling Protocol (Phase 2: 214 remaining)
- **Parallel processing** across 3 APIs
- **Rate limiting**: 5 req/sec (Elsevier), 2 req/sec (others)
- **Error handling**: automatic proxy failover
- **Progress tracking**: real-time status updates

## Download Tracker Schema

### Status Codes
```
code | meaning               | action
-----|----------------------|------------------------
API  | Direct API success   | Success, no further action
PROX | Proxy access success | Success via Vanderbilt
UNPA | Unpaywalled          | Free access
FAIL | All methods failed   | Manual intervention
```

### File Structure
```
gvf_output/
└── institutional_downloads/
    └── 2026-02-07/
        ├── target_pmids.txt
        ├── paper_acquisition_tracker.md
        ├── success_log.json
        ├── failure_log.json
        ├── papers/
        │   ├── elsevier/
        │   ├── wiley/
        │   ├── springer/
        │   └── proxy_fallback/
        └── reports/
            ├── final_summary.md
            └── access_analysis.json
```

## End-to-End Validation Script

### Test Steps
```bash
# 1. Validate API connectivity
python validate_institutional_access.py

# 2. Test 10-paper sample
python test_sample_download.py --count=10 --provider=all

# 3. Generate access report
python generate_access_report.py

# 4. Scale to full queue
python scale_downloads.py --total=224 --batch-size=25
```

## Critical Dependencies

1. **Network**: Vanderbilt VPN for proxy (if off-site)
2. **Time**: Springer key generation (24-48 hours)
3. **Volume**: 224 papers across 3 months of literature (biomedical)
4. **Rate limits**: Respect publisher APIs to avoid bans

## Next Steps (Requiring User Action)

### Immediate (Today)
- [ ] Register for Springer API key using brett.kroncke@vanderbilt.edu
- [ ] Wait for email confirmation (check spam folder)
- [ ] Run validation script to confirm institutional access

### When keys arrive
- [ ] Add Springer key to .env
- [ ] Run full 224-paper queue
- [ ] Monitor progress via tracker

### Proxy Alternative
For non-API publishers:
- Use Vanderbilt browser session
- Export authenticated cookies
- Apply proxy prefix to URLs
- Manual download + metadata extraction

## Contact Information

- **Librarian**: Kayce Gill (kayce.gill@vanderbilt.edu)
- **GVF Support**: brett.kroncke@vanderbilt.edu
- **Emergency**: Text mining assistance available Mon-Fri 8AM-5PM CT