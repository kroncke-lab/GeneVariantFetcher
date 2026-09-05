# Source recovery and count-preserving validation, 2026-09-05

The recovered SCN5A **20031634** article now produces **11 TP / 2 FP / 2 FN**,
versus **0 / 0 / 13** in its old paper-derived lock. Eight matched variants
retain affected and unaffected counts; all 16 supplied A/U values equal gold.
A shared extraction/migration fix preserves an explicit insertion label that
the model read correctly but validation previously discarded. The final code
also includes the general repository downloader and measurement safeguards.

These are deliberately selected, already-opened calibration papers. They
demonstrate specific mechanisms, not a passing population-level discovery or
confirmation. The accepted headline and original gold/locks are unchanged.

## What was tested

The [22-paper audit](../phenotype_failure_panel_20260905/README.md) ranked 336
completed attempts and mixed partial A/U capture, missing capture, and exact
controls. The [16-paper downloader test](../repository_fallback_20260905/README.md)
downloaded three bodies, including both requested examples and one of 14
additional papers. That third paper already had a corpus body.

This implementation campaign ran 15 fresh gene-paper attempts over nine unique
papers: ten prototype attempts spanning opened cohorts, a four-attempt final
ablation, and one check of the second recovered body. No cached predictions or
databases entered extraction. The final ablation has identical initial text
and asset hashes **and identical actual extraction-rendering hashes** to its
four corresponding prototype attempts.

| Locked run | Attempts | TP | FP | FN | Test cost proxy |
|---|---:|---:|---:|---:|---:|
| `20260905_mechanism10_candidate` — initial prototype | 10 | 348 | 41 | 12 | $2.6282 |
| Same prototype, restricted to the four ablation attempts | 4 | 55 | 11 | 6 | included above |
| `20260905_mechanism4_final` — final code | 4 | 55 | 6 | 6 | $2.0381 |
| `20260905_repository1_final` — second recovered body | 1 | 0 | 0 | 20 | $0.0254 |

The four-paper-attempt comparison keeps TP and FN unchanged while removing
five gold extras. It trades one RYR2 TP for the recovered SCN5A insertion TP.
Gold extras are scored under the fixture's unmatched-prediction rule; they are
not automatically false biological statements. Do not pool these overlapping
runs or compare the ten-attempt total to the four-attempt total.

Fresh-model count output also varied across the three gene attempts for
22677073. Those A/U coverage changes are not credited to the insertion rule.
The two exact controls in the ten-attempt run retain RYR2 W4645R **4/2/2** and
C2277R **8/7/1**. Functional control 32533946 retains identities with null counts.
The high-volume 30059973 case still supplies carrier counts without A/U;
its acquired clinical tables need a separate reading intervention.

## Changes retained

1. **General repository acquisition.** Normal harvesting and paywall recovery
   follow DOI-bound Unpaywall/OpenAlex/HAL copies, verify the attached article
   beyond repository covers, retain the PDF and page provenance, and preserve
   body-only/incomplete-supplement status. Publisher retrieval and any enabled
   browser completion precede the body-only fallback. Nested captions,
   supplements, malformed redirects and DOI render glue are handled explicitly.
2. **Preserve an explicit source insertion.**
   `p.1570-Phe1571insIle` becomes `p.1570_F1571insI` in the shared validator.
   The left residue remains unspecified; the original notation is retained.
   Adjacent positions, known residue tokens and existing gene/range checks are
   required. No cDNA allele, carrier count or phenotype is invented.
3. **Freeze acquisition inputs for reading experiments.** Production-eval
   setup copies and hashes selected text and local assets. Gold-free staging
   fails on changed/missing inputs and cannot substitute prior-run sources or
   re-download a retry-marked body. PubMed metadata remains an ordinary input;
   post-run rebinding still verifies the actual rendering consumed.
4. **Correct diagnostic and reporting failures.** FN diagnosis inspects
   identities in structured extraction rows, not narrative mentions. The old
   R594Q claim was a model omission, not a parser drop. All-miss count figures
   now report undefined agreement and an empty table rather than divide by
   zero; the companion difference figure retains every gold miss.

The insertion mechanism is isolated by an identical-response, zero-model
replay: validation previously kept 12 rows, now keeps 13, and migration retains
the added row with **11 carriers / 6 affected / 4 unaffected**. See
[insertion_mechanism_replay.json](insertion_mechanism_replay.json). The fresh
final run independently retains the same identity and counts through trusted
projection. One unknown carrier explains the difference from gold's carrier
total of ten. Other +1 carrier discrepancies have the same source/gold
denominator issue; source counts were not altered to fit gold.

Two 20031634 identities remain missed. The acquired PDF's `c.396312T.C` glyphs
were read as `c.3963+12T>C`, while gold has `c.3963+2T>C`; the duplicated DNA
segment also lacks a fully matching source-derived normalization. Neither was
silently corrected from gold. Three matched identities still lack trusted A/U.

## Changes rejected and remaining source gap

All three speculative parser additions were removed: source-only substitution
promotion, explicit cDNA-column pairing, and gene-labelled prose case alleles.
The first two had no demonstrated incremental win. The prose rule proposed
two gold identities but also five gold extras across KCNQ1/SCN5A; only one of
those extra RYR2 gold matches survived the prototype's final projection. The
final ablation removes the five extras and retains the same aggregate TP.

SCN5A **25163546** is an explicit negative acquisition result: its recovered
16-page PDF still yields **0 / 0 / 20**, with no supplied counts. The body
references separate supplementary material, while the acquired corpus lacks
its variant roster. A body download is not proof that count-bearing components
were recovered. The run cost about 2.5 cents; further model repeats on these
bytes are unjustified. Component-level supplement recovery, the unconsumed
30059973 clinical tables, and 21302287 patient-ID joins remain in `TASKS.md`.

## Review, integrity, cost and reproduction

Claude, Grok and Agy CLI reviews are retained in this directory. Grok and Agy
completed bounded facts-only retries after their initial file reviews hit
limits. [PLAN.md](PLAN.md) records accepted advice and post-lock decisions;
reviewer suggestions are not treated as established findings. In particular,
no wrong-gene quarantine was relaxed, no missing-source paper was removed from
an end-to-end denominator, and no small calibration gain was called significant.

Tests cost **$4.69165 by the repository's dated public-price proxy**. Claude
reported **$3.07193** for its review. Known combined cost is about **$7.76**;
Grok/Agy CLI billing is unavailable. This is a ceiling-controlled campaign,
not an instruction to consume all $100. [budget.json](budget.json) separates
test costs from consultation and the earlier closed $84.035 campaign.

All 34 original audit input hashes were rechecked. Each new run locks fresh
predictions, source bindings, production trace manifests and completed
gold-free status before scoring. The initial final reader fingerprint is
`cfd7fc6a4256f775bd9c8b91d662467c04ed01264df9ef4da921279edcef34bb`.
The all-miss figure repair was made after extraction and locking; it changes
report rendering only. Raw PDFs, frozen source assets, databases, traces and
runtime archives remain local; compact locks, predictions, reports and source
hash manifests are versioned. The retired prototype launcher is historical
evidence and is expected to refuse the current runtime.

Rebuild the compact scored ledger with:

```sh
.venv/bin/python docs/evidence/recall_campaign_20260905/summarize.py
```

[results.json](results.json) records per-paper identity/count metrics, paired
source hashes and artifact digests. Final offline suite: **2,890 passed**.
The all-miss regression and difference-figure checks also passed 17 tests.
Ruff and formatting checks passed across all 25 changed Python files.

Figures: [final four-attempt difference](../../../benchmarks/codex_paper_eval/runs/20260905_mechanism4_final/figures/gold_difference.png),
[all-miss source check](../../../benchmarks/codex_paper_eval/runs/20260905_repository1_final/figures/gold_difference.png),
[canonical stratified figure](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.png).
These were visually checked. The canonical figure remains evidence from its
registered cohorts; these unregistered calibration runs do not change it.
