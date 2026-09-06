**1. Main Design Flaw and Budget Repair**
The fatal flaw is mechanical: reserving the *full output cap* for Astra at $50/M tokens will instantly breach your $2.174 remaining budget if you use standard 8K–32K caps. For example, reserving a 32K cap costs $1.60 *per Astra call*. You cannot run 6–8 tasks safely. A secondary flaw is testing Astra against a single-pass Grok when the operational alternative is Grok-with-retry.

**Repair:** You must ruthlessly constrain the output cap to ~1,500 tokens (reserving only $0.075 per Astra call). Instead of feeding entire 20K-token papers, feed only the strictly relevant methodological/results text blocks and tables. Limit the experiment to exactly 6 tasks to guarantee completion without silent reserve releases.

**2. Minimal Paired Arm and Fair Retry Plan**
Run a 3-arm paired design on the 6 tasks:
*   **Arm A:** Grok 4.6 (Single pass).
*   **Arm B:** Grok 4.6 (Retry/Self-Correction).
*   **Arm C:** Astra (Single pass).

**Fair Retry Trigger:** To prevent gold-leakage, Arm B’s retry must be triggered *deterministically by Arm A’s output payload*, not by comparing it to an answer key. Trigger a retry if Arm A yields:
1.  A strict JSON-schema violation.
2.  A `null` for a mandatory endpoint (e.g., missing count).
3.  An internal inconsistency (e.g., subset counts exceeding the stated total).
If triggered, feed Arm A's output back to Grok 4.6 with a generic prompt: *"Your previous extraction had [Schema Error / Null / Math Inconsistency]. Re-read the source to correct omissions, specifically checking for relatives outside tables or implicit aggregates. Output updated JSON."*

**3. Predeclared Gates, Penalties, and Stopping Rules**
*   **Regression Penalty:** Define a Net Utility Score. A false positive (hallucinating/inventing a patient) is operationally dangerous. Score: `(True Positives * 1) - (False Positives * 3) - (Null Misses * 1)`.
*   **Inclusion Gate:** To justify its massive cost premium (5x input/8x output) and prior latency issues, Astra’s Net Utility Score must beat Arm B (Grok + Retry) by at least **+20% absolute margin** across the 6 tasks.
*   **Stopping Rule:** Evaluate sequentially. Stop the experiment entirely if Arm B (Grok + Retry) achieves a perfect Net Utility Score on the first 3 tasks, or if Astra experiences >1 operational failure (timeout/blank output). The remaining budget is saved, and Astra is disqualified.

**4. Outcomes Justifying Protocol Policies**
*   **Default Astra:** Justified *only* if Astra vastly dominates the Net Utility Score on all tasks without operational failures, and if these derivations actually pass the unchanged production literal validator. (Highly unlikely given prior evidence).
*   **Selective Astra:** Justified if Arm B (Grok) succeeds on standard table extraction but deterministically fails (returns nulls) on complex edge cases like relatives-outside-tables or implicit aggregates, while Astra consistently solves them. *Policy:* Route to Grok first; route to Astra *only* if Grok returns a flagged failure or null.
*   **No Astra:** Justified if Grok + Retry matches or beats Astra’s Net Utility, if Astra's derivations continue to be rejected by the rigid production validator (rendering its reasoning moot), or if Astra repeats its historical operational cap/timeout failures.
*   **What Remains Unknown:** This is a failure-enriched, 6-task sample. You cannot establish population-level accuracy, exact API failure rates at scale, or how either model will perform on entirely novel paper formats unseen in your calibration set. You are only deciding if Astra possesses a unique capability ceiling worth paying for in edge cases.
