**Adversarial Review of Model-Routing Test Results**

**1. Justified Model Use**
These results **do not justify promoting Astra Medium for full-paper extraction.** Astra degraded identity recall from 88.71% (Grok 4.3 control) to 77.16%, missed 91 baseline identities, and increased carrier absolute error by 51.9%. Operationally, Astra is fragile in this configuration: 6 of its 15 primary calls exhausted the 32k-token reasoning budget without emitting visible output, triggering costly retries.

The results *do* justify using Astra exclusively for **bounded, source-linked patient roster reading**. The compact diagnostic test proved Astra Low could accurately transcribe 48 people from a clinical table in under 50 seconds for $0.10. Using Astra as a literal-only clinical overlay is unjustified; the strict validation contract rejected almost all of Astra's raw additions, yielding exactly one accepted count across 12 papers.

**2. Errors in Denominators and Causal Claims**
*   **Precision as an Illusion:** Astra’s "100% exactness" on supplied affected/unaffected (A/U) values is a misleading causal claim regarding model competence. It is an artifact of massive omission: Astra only supplied 37 affected values against 650 nonzero gold targets.
*   **The Deterministic Mirage:** The strong top-line recall is heavily distorted by deterministic extraction shortcuts shared across just two massive papers (20129283 and 30059973). When these shortcuts are excluded, Astra’s recall plummets to 46.77% (vs Grok’s 95.70%). Claiming model parity based on the overall metric is a false causal attribution.
*   **Mismatched Denominators:** The raw diagnostic lane combines incompatible units. It notes 54 rejections for "no unique retained baseline identity" alongside rejections for specific field criteria. The 54 rejections count *variant proposals*, while the others count *fields*. Adding these together into a single rejection denominator misrepresents the true field-level error rate.

**3. The 75%-of-Human Ceiling**
This experiment **categorically does not measure a 75%-of-human ceiling.** The test is a 12-paper, failure-enriched, opened calibration, explicitly precluding any claims about random population holdouts or blinded human equivalence. Furthermore, the workflow itself is bottlenecked by non-model factors:
*   Missing source text (e.g., absent Table 1 in RYR2 18929323).
*   Unfinished responses masked by strict length limits.
*   A rigid acceptance contract that drops valid new variant proposals because they lack a pre-existing retained baseline identity.
*   Unresolved endpoint ambiguity (e.g., models failing to distinguish "No Dx" table footnotes from ECG-suggestive prose).
These systemic barriers cap extraction success long before the model's intrinsic reasoning ceiling can be accurately measured.

**4. What the Next Test Should Isolate**
The follow-up test must abandon full-paper zero-shot extraction and isolate a **validated patient-record aggregation path**. Specifically, it should measure:
*   **Bounded table/roster reading:** Isolating the extraction of explicit patient rows (with headers, captions, and footnotes included) from prose synthesis.
*   **Deterministic aggregation:** Testing a pipeline that preserves patient IDs, deduplicates index cases across repeated tables, and calculates phenotype sums strictly in code rather than relying on LLM arithmetic.
*   **Error decoupling:** Isolating reader timeouts and budget guards so that a failed full-paper read does not silently block the merging of valid, pre-filtered table hints.

**5. Honesty of Budget and Operational Caveats**
The budget caveats are honest but reveal severe workflow inefficiencies. Reporting a returned proxy cost of $38.58 against a massive $81.46 unknown-usage/retry reserve exposes the hidden financial drain of 32k-token silent failures and SDK retries. Acknowledging that these are API proxies and not reconciled Azure invoices, and excluding CLI costs, is transparent and necessary. Furthermore, the admission that earlier Grok pilots were abandoned unscored because the SDK silently dropped reasoning_effort parameters demonstrates high operational integrity in reporting pipeline failures rather than obfuscating them behind bad scores.
