**Do not assign a Grok root cause yet.** Azure `200` plus `499` only shows the client aborted a request the resource later counted as success. That is compatible with deadline, slowness, retry stacking, or a late 200 after a dead client. It does not pick among cap, auth, provider, and retries.

## Bounded Grok differential (four tiny calls, then stop)

All: direct HTTP, **zero retries**, api-key, same deployment, log `x-request-id`, client start/end, Azure status/latency. Prompt ≈ “Reply with ping.” Freeze this matrix before any call.

| ID | Surface | `reasoning_effort` | Output cap | Client deadline |
|----|---------|--------------------|------------|-----------------|
| G1 | Chat | omitted | 16 | 55s |
| G2 | Chat | low (the value previously dropped) | 16 | 55s |
| G3 | Chat | omitted | 16 | 120s |
| G4 | Responses (the other path you actually use) | omitted | 16 | 55s |

**Pre-registered reads (no fifth call):**

- **Auth:** all 401/403. **Endpoint:** 404/route error on G1 xor G4; G1 and G4 both 200 ⇒ path/auth not the tiny-timeout story.
- **Retry amplification:** already removed. If G1 is now 200 inside 55s, prior 499s were client/SDK stacking, not the model. If G1 still 499 and Azure 200 **same request-id after abort**, it is deadline vs server finish, not retries.
- **Slow provider:** G1 499/timeout and G3 200; join request-ids. Zero retries already, so this is TTFB/completion vs 55s, not amplification.
- **Output/reasoning cap:** G1 200, G2 timeout, empty, or `length` stop. Same 16-cap, only effort differs. Do **not** add a 4096+effort probe unless G2 fails and budget remains; that call is larger and confounds the tiny health question.

Large historical successes vs tiny failures are **not** a fifth cell. Record them as unmatched size/route, not as proof.

If G1–G4 still fail, **do not spend Astra budget on Grok papers.**

## Expanded Astra test: gaps

Seven packets × low × medium at reserved 4096, 120s, two prices, is not a differential. It is a spend. `$4.90` with Astra `$10/$50` per M, full output reserved, failed usage **unknown**, and `$81` “unknown retained” **unusable**, cannot underwrite 14 completions. Worst-case output alone is in the same order as the hard limit before source bytes.

**Confounds you are about to bake in:** 55s→120s and SDK→HTTP change two variables versus the Grok incident. Low vs medium is **not** “identical except effort” if reasoning tokens compete with the 4096 cap; medium will look worse on roster completeness when it is actually a length stop. Structurally parsed tables must not be gold-shaped; if the parser was tuned on reference rows, the model is not “transcribing source.”

Gold integer match **must not accept.** A tractable contract, all required:

1. Operational class first: timeout, 4xx/5xx, parse fail, `length` stop, schema fail (not a count).
2. Every person row binds a canonical source row/span + provenance; no IDs minted from totals.
3. C/A/U vs **source text** and vs **frozen reference**, separately; additional and missing both count.
4. Duplicate IDs / co-carrier collisions.
5. Positive-value coverage on source-positive fields, not whole-roster accuracy.
6. Aggregate-only / abstract-only / count-free: **accept iff roster empty** and published totals stay in the aggregate bucket.
7. Frozen endpoint: carrier ≠ diagnostic phenotype; asymptomatic ≠ phenotype-negative; later symptom freedom ≠ baseline-negative; controls ≠ carriers. Suggestive ECG ≠ Unaffected even if no footnote says “healthy.” Derived counts follow the table-defined endpoint.

Keep the production literal validator frozen. No promotion from this calibration.

## Highest-value fixed comparison

After G1–G4: **two Astra packets, low only, 4096, 120s, one in-flight, abort if usage unknown or remaining budget < reserved output.**

- 15-person multi-variant roster/footnotes (identity, footnotes, co-carriers, C/A/U).
- Aggregate-only missing table (must not invent persons).

That pair tests extraction vs fabrication. If both operationally succeed and budget remains, add **one** medium clone of the 15-person packet (effort × cap), then the 8-person diagnostic-vs-symptom packet (endpoint), not the other four. DOC48, MYBPC3, functional, abstract stay frozen packets for a later envelope.

Do not report human-ceiling percentages or population uplift. Report operational fails, binding, identity, source vs reference C/A/U, duplicates, and positive-value coverage on this frozen set only.
