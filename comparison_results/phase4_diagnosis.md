# Phase 4 — Paywall / access diagnosis (2026-05-04)

## 4.10 — Wiley TDM key

**Status: key authenticates but is consistently denied (HTTP 403).**

Test against the live Wiley TDM endpoint with the configured key
(fingerprint `a4f6...35c7`, len 36):

| Probe | Status | Notes |
|-------|--------|-------|
| `GET /tdm/v1/articles/10.1002/humu.10131` (with `Wiley-TDM-Client-Token`) | **403 Forbidden** | empty body, no `WWW-Authenticate` |
| `GET /tdm/v1/articles/10.1002/jcb.10412` | 403 | empty body |
| `GET /tdm/v1/articles/10.1002/cphy.c150042` | 403 | empty body |
| same URL with no token | 400 | endpoint validates token presence |
| same URL with `CR-Clickthrough-Client-Token` (legacy) | 403 | empty body |

The endpoint **does** recognize the token (omitting it gives 400, not 403) — it
just refuses to authorize it for any DOI. There is no `WWW-Authenticate` header
explaining why, which is the Wiley-typical signature for a revoked or
out-of-scope agreement.

**Action for Brett:** the key needs to be renewed through the Wiley TDM portal.
Two likely causes:

1. The institutional TDM agreement has expired or been revoked at the library
   level — confirm with VU library liaison.
2. The token's per-key allowance has rolled over and Wiley auto-suspended it —
   re-request a new token via <https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining>.

Once renewed, drop the new value into `.env` as `WILEY_API_KEY=...` and
the orchestrator will pick it up automatically — no code changes needed.

For reference, a working request looks like:

```
curl -L -H "Wiley-TDM-Client-Token: <new-key>" \
  https://api.wiley.com/onlinelibrary/tdm/v1/articles/10.1002/humu.10131
```

A 200 with PDF/XML body confirms authorization. Keep the test DOI handy
because Wiley's status page does not surface TDM-specific outages.

## 4.11 — Browser-fetch headless feasibility

**Verdict: feasible overnight, with one-time setup.** The wiring is all in
place; what's missing is a dedicated, pre-authenticated Chrome profile.

### What `cli/browser_fetch.py` already supports

- `--headless` flag → uses `chromium.launch(headless=True)` /
  `launch_persistent_context(headless=True)`.
- `--use-profile` + `--profile-path <dir>` → persistent Chrome profile.
  Cookies, saved logins, OpenAthens / institutional SSO state all survive
  across runs.
- `--use-claude` → Claude-assisted PDF link discovery on complex sites.
- Stealth patches (mask `navigator.webdriver`, plugin spoofing, language
  spoofing) — applied in both headed and headless modes.
- Cloudflare detection with loop avoidance.

### What blocks headless today

1. **`_manual_download_fallback` is disabled in headless mode**
   (browser_fetch.py:853-857). When Cloudflare or a captcha blocks
   automation, the headless run abandons the PMID instead of waiting for
   human intervention. Acceptable trade-off for overnight: those DOIs
   would just stay in the "manual review" bucket.

2. **Default profile path conflicts with Brett's daily Chrome.**
   Without `--profile-path`, the script reaches for
   `~/Library/Application Support/Google/Chrome`, which is locked while
   Brett's main Chrome is open. An overnight headless run should NOT
   share this profile.

3. **First-time auth.** Vanderbilt OpenAthens + most publisher SSO flows
   need an interactive login on first contact with a profile. A fresh
   profile launched in headless mode cannot complete that. Hence the
   one-time interactive setup below.

### Setup recipe for overnight headless

This is the runbook Brett should follow once. After this, the headless
overnight loop is `gvf browser-fetch ... --headless --profile-path ...`.

```bash
# 1. Create a dedicated profile dir for GVF (does NOT touch daily Chrome)
mkdir -p "$HOME/Library/Application Support/GVF-Chrome-Profile"

# 2. INTERACTIVE one-time auth — log into every publisher Brett needs
.venv/bin/python -m cli.browser_fetch \
  results/KCNH2/20260503_161203/pmc_fulltext/paywalled_missing.csv \
  --use-profile \
  --profile-path "$HOME/Library/Application Support/GVF-Chrome-Profile" \
  --pmids 10220144  # any single Wiley DOI to drive the SSO flow

# In the visible browser:
#   - Click through Vanderbilt OpenAthens login.
#   - Visit Wiley, Elsevier ScienceDirect, Springer, Oxford Academic at
#     least once each so each publisher drops a session cookie.
#   - When all publishers show "VU subscriber" badging, close the browser.
# The cookies persist in the GVF-Chrome-Profile dir.

# 3. OVERNIGHT — same command, but headless, full PMID list, with --use-claude
nohup .venv/bin/python -m cli.browser_fetch \
  results/KCNH2/20260503_161203/pmc_fulltext/paywalled_missing.csv \
  --headless \
  --use-profile \
  --profile-path "$HOME/Library/Application Support/GVF-Chrome-Profile" \
  --use-claude \
  > overnight_paywall_fetch.log 2>&1 &
```

### Expected hit rate

`paywalled_missing.csv` from the May-3 run: **490 distinct PMIDs**, breakdown
by publisher prefix:

| Publisher | Count | Headless after SSO? |
|-----------|------:|---------------------|
| AHA / LWW (10.1161) | 42 | yes (subscription works without captcha) |
| Wiley society (10.1111) | 31 | yes (after WILEY TDM key renewal *or* SSO cookie) |
| Elsevier (10.1016) | 40 | yes via SSO; or via Elsevier API if `ELSEVIER_API_KEY` valid |
| Springer (10.1007) | 21 | yes (Springer API + SSO fallback) |
| Wiley (10.1002) | 17 | yes after SSO |
| Oxford (10.1093) | 11 | yes after SSO |
| Karger (10.1159) | 7 | **no — Cloudflare blocks** (known; documented in CLAUDE.md) |
| Other long-tail | ~25 | mixed |
| No-DOI / failed-resolution | 720 *(rows, fewer distinct PMIDs)* | requires manual triage |

**Realistic estimate: 60-75% of the 490 paywalled PMIDs (≈300-370) become
fetchable** in one overnight run, modulo Cloudflare-walled publishers
(Karger and ~30 long-tail) which need manual fetch.

### Recommendation

Do the one-time SSO setup tonight, then schedule the headless run as a
nightly cron entry against `paywalled_missing.csv`. Each new
extraction cycle regenerates the file, so the queue self-shrinks as
papers are recovered.

If `--use-claude` is too slow or noisy, drop it — the deterministic
publisher selectors hit the common cases without LLM assistance.
