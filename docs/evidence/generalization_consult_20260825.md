# Generalization consult — Grok 4.6 xhigh + GPT-5.6 Sol max (2026-08-25)

Goal: raise precision, accuracy, recall, and MAE on gold-120 and gold-150 with
changes that are **general** — they must work on a new gene-disease pair where
nobody knows the answers in advance. Two independent adversarial reviews were
run against a self-contained brief built from the authoritative lock
`benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118`
(554 TP / 283 FP / 78 FN). Sol had read-only repo access and cites `file:line`;
Grok ran advisory-only.

Every claim below was reproduced locally before being recorded. Claims that did
not survive reproduction are in §6.

## 1. Benchmark contamination: gold-derived aliases in the runtime path

`gvf_data/kcnh2_variant_aliases.json` declares its own provenance:

```json
"source": "Gold standard Excel + generated forms"
"total_variants": 856, "total_aliases": 4766
```

It is loaded at runtime, not only in benchmarks: `_load_gene_aliases` /
`_lookup_alias` are consulted for every gene inside `normalize_variant`
(`utils/variant_normalizer.py:1588`), and `gvf_data/README.md:9` describes the
file as "loaded by `utils/variant_normalizer.py`". KCNH2 is one of the four
gold-120 genes, so KCNH2 normalization has access to a 4,766-entry map derived
from the answer key. It is the only gene with such a file.

**Reproduced:** re-matching all 554 locked pairs with `_lookup_alias` stubbed to
return `None` gives **554/554 — zero pairs depend on the alias map at score
time**. The headline is therefore not manufactured by score-time lookup. What
cannot be recovered from already-normalized predictions is the map's influence
during *extraction*.

This is the single clearest violation of the "must work without knowing the
answers" requirement. It does not invalidate the current score, and it is not a
reason to delete the file — an alias map is legitimate infrastructure. The
requirement is that its contents be derived from public notation resources
rather than from the evaluation key, and that the provenance be stated where
the file is loaded.

## 2. Gold-118 is a calibration set, not an unbiased estimate

Both reviewers reached this independently and it is the load-bearing protocol
conclusion. The set has been through many sealed extractions and a long chain of
residual-driven matcher and table patches. Seventy-eight FNs remain on 23
papers, half of them in three papers. A rule whose only motivating examples are
those 23 papers will read as general and still be a test-set fit.

Treat 554/283/78 as a **frozen non-regression lock**, not a performance
estimate. The first defensible generalization number will come from gold-150
holdout, scored once.

## 3. The headline mixes two different tasks

The locked projection includes every source layer, though the exporter already
distinguishes paper readings from ClinVar/PubTator linkage
(`benchmarks/codex_paper_eval/db_to_predictions.py:129`). Sol's replay of the
locked DB:

| Lane | TP / FP / FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: |
| Linkage-assisted (current headline) | 554 / 283 / 78 | 66.19% | 87.66% | 75.43% |
| Paper-derived only | 512 / 125 / 120 | 80.38% | 81.01% | 80.69% |

Linkage FPs by gene: KCNQ1 79, SCN5A 72, KCNH2 7, **RYR2 0**.

This displaces the previous explanation of RYR2's cleanliness. RYR2 looks clean
because it carries no linkage-layer FPs, not because CPVT literature is tidier.
It also weakens "all five top-FP papers are gold-incomplete": 95 of those 162
extras are ClinVar linkage rows, including **all 63** on KCNQ1 PMID 19632626.

Report both lanes, with paper-source as primary for the extraction task and
linkage-assisted as the enrichment product. Do not delete linkage rows.

Two further naming corrections:

- Counted-extra 97.88% is `554/(554+12)` — every matched row in the numerator,
  only count-bearing extras in the denominator. It is a Gate-2 score, not a
  conventional precision. Stated transparently: **2.17 counted extras per 100
  matched rows**. The conventional count-bearing precision is 210/(210+12) =
  **94.60%**.
- Raw 66.19% is precision against an incomplete fixture under unmatched=false.
  Because some gold rows are themselves out of scope (§5), it is not even a
  clean lower bound.

## 4. Conditional MAE is the wrong decision metric

Carrier MAE 0.193 is computed over the 207 rows the pipeline dared to answer.
Non-zero carrier coverage is 43.2%. End-to-end error is already implemented
(`cli/compare_variants.py:2828`) but not foregrounded:

| Field | Conditional MAE | End-to-end MAE |
| --- | ---: | ---: |
| carriers | 0.193 | **1.541** |
| affected | 0.321 | **1.528** |
| unaffected | 0.321 | **0.512** |

Abstaining improves the first and cannot improve the second. Both reviewers
independently proposed an abstention-charging loss; preregister it for gold-150
**before** holdout opens, since the abstention cost is a value judgement:

```text
supplied prediction:  |prediction - gold|
NULL or missed:       max(1, gold)
count-bearing extra:  max(1, predicted count)
```

Do not pool carriers, affected, and unaffected without preregistered weights.
Keep the gold-zero-versus-NULL convention gap out of the number and report it
separately — forcing the pipeline to emit 0 to match gold's convention would
violate the no-inference contract that was already measured and kept.

## 5. Several "false negatives" are scope errors, not extraction failures

Verified against locked source:

- **PMID 15898185** (10 FN) is an *editorial*, not the cohort paper. It mentions
  S1103Y and R1193Q as commentary. Not a whole-paper extraction failure.
- **PMID 22966897** (2 FN, 3 FP) is an editorial compiling prior molecular-
  autopsy studies. Both its gold rows and its predictions are prior-study
  material.
- **RYR2 K4481R** was explicitly reported as unconfirmable and judged an FFPE
  artifact by the source.
- **KCNQ1 A340E** is the mouse ortholog; the paper states the human residue is
  A341E.
- **SCN5A p.F1617del** is, by contrast, explicit current-human evidence and a
  genuine recoverable miss.

A source-blind status ledger (`source_complete`, `abstract_only`,
`supplement_missing`, `editorial/prior-study`, `mixed_human_model`,
`unresolvable`) should be recorded before scoring. Keep end-to-end recall
primary; report source-observable recall only as a secondary diagnostic. Never
"correct" gold by deleting only the misses — editorials contain matched and
extra rows too.

## 6. Hypotheses that did NOT survive reproduction

Recorded so they are not re-proposed.

- **"No quote captured" means fabricated.** False. No-quote rows are 157/218
  (72%) of carrier-bearing rows but only 43.5% of the 23 carrier errors — about
  **3× less** error-prone than quoted rows, because they come from deterministic
  table/figure parsing rather than the LLM text lane.
- **`carriers=1` is a fabricated singleton default.** False. 128 of 218 rows
  predict 1 and only 6 of them err. Most variants genuinely have one carrier.
- **A systematic "+1 for the proband" rule.** Both reviewers flagged this
  independently as the easiest way to overfit; the −1 pattern is real but small
  (10 rows) and heterogeneous. SCN5A 22882672 (8→1) is an abstract-versus-family
  scope problem, not an off-by-one.
- **Adjacent-codon merge** of `p.Ile1660Val c.4978A>G` with
  `p.Ile1659Val c.4975A>G` (SCN5A 29709101). Both reviewers rejected it: these
  are different nucleotides and different codons unless a transcript proves
  otherwise. Merging on ±1 codon is gold-laundering and is the shape of the
  already-rejected ClinVar drop.
- **Table-as-`*Figure: Table*` placeholders explain the count gap.** Measured:
  they explain only **31 rows on 2 papers**. Carrier coverage by capture mode is
  51.5% with real pipe tables, 35.8% with no tables, 0% placeholder-only.
- **`p.R360_Q361dupQKQR` should collapse.** Already deliberately refused with a
  correct reason at `utils/variant_normalizer.py:1086`: a four-residue insertion
  across a two-residue range is a different allele.
- **584 of 789 extras sit on fully-recalled papers.** Historical figure. The
  current lock is **237/283 = 83.7%**.

## 7. Count coverage, decomposed

Carrier coverage is 207/632 (32.75%), but the addressable part is smaller than
it looks. Splitting gold's explicit zeros from real attribution misses:

- Total gap 425 rows; **real non-zero gap 268 rows**.
- SCN5A 32533946 contributes 83 gold assertions and **all are zero** — the
  convention gap, not a miss.
- Concentration: KCNQ1 14678125 (41), RYR2 19926015 (40), RYR2 19398665 (27),
  KCNH2 29650123 (22, acquisition-blocked), KCNQ1 21956039 (14).

The top papers matched most or all variant identities and supplied **zero**
counts, so this is a whole-paper count-lane failure, not per-row inaccuracy.
Inspection of their captured source shows the per-variant carrier table is
frequently absent from the capture entirely (KCNQ1 14678125 is 13.7 KB with no
tables at all), which makes this largely an acquisition problem rather than an
attribution one.

## 8. Ranked, general, $0 code items

Neither reviewer authorized another paid gold-118 extraction, and neither
recommends one. In agreed priority:

1. **HGVS range-separator repair.** `c.2550-2551insTG` does not match
   `c.2550_2551insTG` (reproduced: `matches()` returns False). Repair `N-M` to
   `N_M` only when both sides are absolute cDNA coordinates immediately followed
   by `del|dup|ins`, and never reinterpret an intronic minus. `c.999-424_1338+81del`
   and `c.2550-2A>G` must be unaffected. Note `c.2550_2551insTG` already matches
   `p.Phe851fs c.2550_2551dup`, so this collapses three spellings to one.
2. **Refuse malformed protein ranges.** `p.Gln359_Arg`, `p.Lys362_His` are
   truncated parses, not alleles. Repair only when the same row carries a
   complete, resolvable cDNA edit; otherwise drop the protein token.
3. **Reference-backed 3′ shift** for repeated-residue deletions, using the
   existing accession-bound references (`gvf_data/reference_sequences/manifest.json`,
   `pipeline/reference_validation.py:150`). Fail closed unless the reference
   proves the repeat. Predicted standalone effect: **555/279/77**, raw P 66.55%.
4. **Provenance-lane scorecard** (§3) and **end-to-end count loss** (§4).
5. **Claim-level observation status** (§5) as a projection, never as deletion of
   canonical rows.

Two existing rules need gold-blind revalidation because they were threshold-
selected on cardiac gold: the codon-shadow gate
(`db_to_predictions.py:274` records that the exact-codon rule was chosen because
it removed 55 gold FPs with zero gold TP loss) and the deletion-span acceptance
at `run_eval.py:1741`.

Also noted: `LOCK.json` binds predictions and selection but **not** the scorer,
normalizer, configuration, or reference files, so a lock does not currently pin
the evaluator that produced it.

## 9. What was executed from this consult

- Cross-gene firewall defect in the gold-150 preregistration found, reproduced,
  and fixed. See `benchmarks/curated_extraction_eval/gold150_preregistered_20260825_amended/README.md`.
- `scripts/build_curation_packet.py` gained `--split-lock` and
  `--canonical-source`; the gene-in-hash assignment became `_split_assignments`
  with pin support and exact preservation of per-gene calibration size.
- `scripts/audit_split_firewall.py` added, with tests, to fail on either defect
  class in any future cohort.

Items in §8 are not yet implemented.
