# Independent validation notes

## Disposition

**Share with caveats.** The run is reproducible and its arithmetic, hashes, evidence
ledger, timing, and API token telemetry pass independent checks. The nominal
precision and recall should not be interpreted as pure model error rates because
several downloaded sources are incomplete or wrong, and the recall-oriented gold
packet omits source-supported variants.

This 20-paper fixed comparison sample is an experimental Codex evaluation. It does
not replace or update the repository's canonical four-gene recall baseline.

## Integrity and arithmetic checks

- The extraction manifest exposed PMIDs, source paths, source hashes, and available
  representations, but no gold values or gold row counts.
- This historical run predates all-representation hashing: its selection records
  hashes for the primary markdown sources, not artifact JSON, PDFs, or figures.
  New runs hash and verify every selected representation before extraction and
  again before locking.
- `predictions.json` was made read-only and SHA-256 locked before `score` opened
  gold values.
- Selection SHA-256:
  `0dd204263b8cebab5815713bff30fa96ac622900a5a4abd90f16c80db0263955`.
- Prediction SHA-256:
  `73210ffde0b9307242d3d3803fe1fa207a5023143889b7e66ee1f5cab2f564dd`.
- Independent recomputation reproduced 79 TP, 15 FP, 38 FN, 117 gold variants,
  94 predictions, every per-gene and per-paper precision/recall/F1 value, all
  count-recall denominators, weighted MAE/RMSE, token sums, and elapsed-time sums.
- `evidence.csv` has one row for each of the 94 locked predictions.

## Matcher audit and sensitivity

The production-style matcher initially scored 72 TP, 22 FP, and 45 FN:
76.6% precision, 61.5% recall, and 68.2% F1. A blinded-to-model-output audit of
unmatched pairs identified seven certain notation equivalents, after which the
audited score became 79 TP, 15 FP, and 38 FN: 84.0% precision, 67.5% recall, and
74.9% F1.

The seven scorer-only corrections cover:

- `ΔKPQ` versus `p.K1505_Q1507del`.
- Verbose `G189ins` insertion notation.
- Protein frameshift forms with `/length` versus `fsX`.
- Optional `c.` prefixes and spacing in cDNA notation.
- Equivalent exon-deletion notation.

The locked predictions were not changed. The source-supported `N4178S` versus gold
`N4168S` pair was deliberately left unmatched because those are different residues.
Use 61.5–67.5% recall and 76.6–84.0% precision as the scorer-sensitivity interval
when comparing this run with results produced by another matcher.

## Nominal false-negative classification

All 38 nominal false negatives have a concrete benchmark or acquisition limitation:

| class | variants | evidence |
|---|---:|---|
| Missing body/table or missing variant nomenclature in the downloaded source | 35 | RYR2 23595086 (18), RYR2 12093772 (12), KCNQ1 10220144 (3), KCNH2 30844837 (1), KCNQ1 27000522 (1) |
| Engineered non-human variant included in gold | 2 | SCN5A 16054936 `Q1476K`; KCNQ1 18567635 `F340W` |
| Source-versus-gold nomenclature conflict | 1 | RYR2 29925740 source `N4178S`; gold `N4168S` |

Details:

- PMID 23595086 is labeled as main text locally but contains metadata and
  references rather than the article body or the variant table.
- PMID 12093772 names only `S2246L` and `G3946S`; the other 12 gold identifiers
  are not present in the recovered text or tables.
- PMID 10220144 contains the title, abstract, and references; the abstract names
  four variants but not the three nominal misses.
- PMID 30844837 resolves to publisher permission guidelines rather than the
  article containing `D609G`.
- PMID 27000522 reports one KCNQ1-positive case in an abstract but provides no
  variant nomenclature, so the extracted unspecified variant cannot match `L266P`.
- PMID 16054936 studies the human familial variant in SCN1A; the SCN5A homolog was
  engineered for transfected-cell functional assays.
- PMID 18567635 uses `F340W` in Xenopus oocytes rather than reporting a human
  carrier.

This classification does not prove perfect reasoning. It does show that nominal
recall is dominated by source and benchmark validity on this sample.

## Extra-prediction audit

The 15 unmatched predictions are not equivalent to 15 hallucinations:

- Nine rows from KCNH2 PMID 15851119 are present in the source's patient table but
  absent from the count-curated gold packet.
- Three rows from SCN5A PMID 23414114 are present in the source table but absent
  from gold.
- SCN5A `D1823D` is source-supported, although its carrier count is indeterminate.
- KCNQ1 PMID 27000522 yields an unspecified variant and is therefore unscorable
  against variant-level gold.
- RYR2 PMID 29925740 explicitly gives `c.12533A>G (p.N4178S)`, while gold gives
  `N4168S`.

Consequently, the reported precision is a conservative "precision against this
recall packet," not an adjudicated clinical false-positive rate.

## Count disagreements

The largest numerical errors also have source-versus-gold disagreement:

- SCN5A PMID 17544529 reports `H558R` in five sudden unexplained death cases;
  the prediction records five carriers/affected, while gold records one.
- KCNH2 PMID 18452873 identifies the insertion in patient 6 with long-QT syndrome
  and atrial fibrillation; the prediction records affected 1/unaffected 0, while
  gold records affected 0/unaffected 1.
- KCNH2 PMID 16043162 describes two unrelated affected Brugada patients with
  `G873S` and `N985S`; the prediction records one carrier/affected for each, while
  gold records zeros.
- RYR2 PMID 12093772 describes two sporadic affected cases for each named variant;
  the prediction records two, while gold records one.

MAE and RMSE remain conditional on the 63 matched assertions with a supplied value.
They must be read with 53.8% count recall; abstentions do not contribute zero error.

## Cross-harness comparison caveat

Claude's fixed 20-paper comparison used the same PMIDs but a different model
interface, prompt, tool loop, token-accounting boundary, and scorer. Its reported
variant recall was 58% for Fable 5, 56% for Sonnet, and 50% for Haiku. Those values
are useful directional context only. The Codex run's exact 217,579 tokens are API
usage for the router and extractor responses; Claude's roughly 0.90–1.14 million
tokens include its agent/tool context and are not an apples-to-apples cost measure.

Claude's separate 48-paper high-carrier stress set was selected using gold carrier
volume and is therefore not a population-blind sample. Its 74% Fable, 73% Sonnet,
and 65% Haiku recall should not be merged with this fixed 20-paper evaluation.
