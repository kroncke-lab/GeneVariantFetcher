## Grok4.6 differential (no loops, ≤4 new cells + 1 confirmatory)

The prior data already has an anomaly you haven't flagged: tiny prompts timed out at 55s while large paper prompts *sometimes* succeeded. That's inverted from what a naive "big-payload-is-slow" story predicts, and it's the strongest clue you have — don't discard it before designing new cells.

Vary two factors orthogonally, one call each, same deployment/api-key/directHTTP/zero-retry/55s deadline:

| Cell | Payload | Effort |
|---|---|---|
| 1 | trivial ping | low/none |
| 2 | trivial ping | medium |
| 3 | large paper | low |
| 4 | large paper | medium |

Read pattern, not single results:
- **1 succeeds, 2 fails at ~55s** → stall is effort-driven token generation, independent of input size. Points at output/reasoning-cap behavior (e.g., reasoning tokens not bounded the way you expect, or this SKU is slow specifically at higher effort), not endpoint/auth.
- **3 succeeds (replicates prior "sometimes")** while **2 fails** → effort, not size, is the discriminating variable — undercuts any "large payload" theory you might otherwise reach for.
- **All four fail identically at ~55s** regardless of size/effort → look at endpoint/auth/network (wrong deployment string, TLS/handshake, throttling), not model compute.
- **Azure metrics split matters**: 499s = client canceled before server finished (server may have been about to return 200 late = slow provider, not broken endpoint); 200s = model actually answered within window on some calls but not others (rules out hard endpoint failure, points at variable latency or concurrent-capacity contention).

**Retry amplification is a distinct hypothesis you cannot confirm with the zero-retry design above** — by construction you've removed the variable. Add one confirmatory cell only: re-run cell 2 or 4 with SDK default retries enabled, single call, and inspect Azure request logs for multiple overlapping upstream calls in that window. If you see duplicate/concurrent requests consuming the capacity=50 quota, that's amplification; if you see one clean request that simply runs long, it isn't. Don't generalize from more than this one cell — it's a check, not a study.

## Astra expanded test — gaps

1. **No repeats per cell.** Single draws per paper × effort level can't separate stochastic model variance from an effort effect. Acceptable given budget, but the report should say "single-draw, not variance-estimated" rather than imply a stable low-vs-medium delta.
2. **Ambiguous failure bucketing.** You need the operational/content-failure split written into the contract *before* grading — otherwise a malformed-JSON or refused response could get scored as a wrong roster instead of an operational failure, silently inflating or deflating the effort comparison.
3. **Aggregate-only / abstract-only cases need an explicit correct-output rule** — the contract should state the expected answer is "decline/flag as non-extractable," not a best-guess roster. Test whether each arm correctly abstains; a model that fabricates a plausible-looking roster there is a distinct, and worse, failure mode than a miscount on a resolvable case.
4. **Blind grading.** Whoever adjudicates C/A/U against source/reference shouldn't know which effort arm produced the output — otherwise expectation bias (medium effort "should" be better) contaminates judgment calls on ambiguous rows.

## Budget check
Rough worst-case: 7 papers × 2 effort levels on Astra (~5k in / 4096 out) ≈ $3.5–3.6; four Grok diagnostic cells + one retry-confirm cell (small in / 4096 out cap) ≈ $0.10–0.15. Total ≈ $3.7, leaving thin headroom under your $4.90 hard limit once failed-call "unknown reserve" usage is counted. Recommend explicitly capping output tokens on the two low-value cases (count-free functional section, abstract-only) since they don't need full 4096 to produce a correct "no roster" answer, freeing margin for the retry-confirm cell.
