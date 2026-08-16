# Affected / unaffected value precision (agent briefing)

Last reviewed: 2026-08-16.

Current-facing ranking of how to raise **exact-match precision of emitted
affected and unaffected integers** on gold-120. This is not Gate 2
counted-extra identity precision and not carrier MAE.

`TASKS.md` owns what to do next. Dated paper-level tables live in
`benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/`.

If another note says "improve precision" without naming a denominator, ignore
it. Identity precision is [`PRECISION_COST_LEVERS.md`](PRECISION_COST_LEVERS.md).
This file is the phenotype-value denominator.

## 1. The metric

On a matched gold row, when **both** gold and the prediction supply a value for
the field:

```
exact-match precision = exact / supplied
```

MAE is the same pairs. Null vs gold-0 is **not** an error here (gold encodes
"none reported" as 0; the pipeline abstains with NULL). That convention gap is
count-recall, not this precision number.

**Gold never abstains.** All 6971 rows in `*_recall_input.csv` supply
`carriers`, `affected` and `unaffected` — 865 / 1493 / 6126 of them are an
explicit `0`. So `supplied` is decided entirely by the pipeline: every integer
it emits on a matched row enters the denominator, and abstaining on a row costs
this metric nothing.

**That does not make "abstain more" free.** Nulling a bucket that is *more*
accurate than the field average lowers precision. Break-even for a candidate
rule is `destroyed / killed < exact / errors`:

| Field | Current | Break-even destroyed-per-killed |
| --- | --- | ---: |
| affected | 40/51 | 3.64 |
| unaffected | 27/40 | 2.08 |
| carriers | 135/161 | 5.19 |

A rule can therefore raise this number while throwing away correct counts. Such
a rule is rejected below as a **product** regression (count recall), not because
the arithmetic disfavours it. Say which one you mean.

## 2. Where it stands

Free rescore of the locked `20260813_gold120_verticalfix` predictions against
live gold and the current matcher — **no new extraction, no LLM calls**:

```bash
scripts/phenotype_value_precision.py \
  --run-dir benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815 \
  --compare --errors-out /tmp/value_errors.tsv
```

| Field | Unguarded | Guarded (current) |
| --- | --- | --- |
| affected | 51/74 — 68.92%, MAE 0.905 | **40/51 — 78.43%**, MAE 0.549 |
| unaffected | 37/56 — 66.07%, MAE 0.964 | **27/40 — 67.50%**, MAE 0.750 |
| carriers | 135/161 — 83.85%, MAE 0.292 | 135/161 — 83.85%, MAE 0.292 |

Identity on the same rescore: TP 545, FN 88, recall 86.10%, raw precision
40.85%, counted-extra precision **98.55%** (8 counted extras). The guard itself
removes one counted extra by nulling the counts that made it counted.

Locked parent `report.json` was not rewritten. This is not a new Gate 2 lock.

## 3. What the guard does

[`pipeline/phenotype_count_guard.py`](../pipeline/phenotype_count_guard.py),
always-on, gold-free, no PMID/variant/gene special cases. Wired in
`pipeline/steps.py` (both persist sites) and
`harvesting/figure_variant_reader._parse_response`. Not behind
`COUNT_CLASSIFIER_POLICY`.

1. **Family copy.** `affected == carriers >= 2`, `unaffected in {0, None}`, no
   distinct phenotype column → clear affected and the implied `una=0`.
2. **Figure copy.** Same pattern on `source_layer=figure` for any N, plus a
   figure `affected` with no carrier total. Vision copies pedigree symbols onto
   affected (16/16/0, 7/7/0, 4/4/0).
3. **Partition does not close.** All three integers emitted and
   `affected + unaffected != carriers` → clear **affected only**. Kills 3,
   destroys 0. The contract already asks the model to run this self-check;
   this enforces it on the emitted integers. Overflow *and* underfill both
   fire: there is no `unassessed` slot in the schema, so a short partition is
   indistinguishable from a miscount.
4. **Unsourced zero affected.** `affected == 0` with no sourced affected column
   → clear affected. A counted zero is a positive clinical claim ("this family
   is entirely non-penetrant"), not an abstention. Kills 1, destroys 0.
5. **Leave real splits alone.** `affected != carriers` with a closing partition,
   or a positive sourced `unaffected`, is a counted partition.
6. **Leave one-proband 1/1/0 text/table rows alone.** Often true case reports;
   wiping them moves this metric the wrong way.

Rules 3 and 4 are deliberately **field-scoped to `affected`**. On the lock the
companion `unaffected` is exact on two of the three failing triples, so nulling
the pair destroys counted values the paper really did report. `unaffected = 0`
is likewise never refused by rule 4: it is the ordinary single-proband shape
and is 16/20 exact.

## 4. Rules that look attractive and fail

Measured on the same lock. Do **not** re-land these:

| Candidate | What happened |
| --- | --- |
| Clear every implied `una=0` including n=1 case reports | Lost 30 exact zeros; una only 66.1%→66.7% |
| Clear **all** figure affected/unaffected | Affected **fell** (lost real 2/1/1 splits). After the copy guard the remaining figure affected is 4/6 — the rejection reason is now "it deletes the only correct pedigree splits", not the old 68.4% number |
| Clear figure `unaffected` only | Kills 4, destroys 3. Precision-positive, but strips una off correct `2/1/1` rows |
| Clear n=1 `affected=1` unless the quote has status words | Affected **fell** 68.9%→68.0%, lost 49 exacts. Lock quotes are too short to be a status oracle |
| Null on a functional/`in vitro`/resource `source_location` | Kills 6, destroys 5 — and it cannot work: KCNH2 18675227 "Electrophysiologic characterization" holds both an exact `3/1/2` and an error `2/1/1`. A section heading is not a study-type classifier |
| Null carriers on `source_location == "Abstract"` | Kills 6, destroys 17. Abstracts are majority-correct n=1 case reports |
| `unaffected = carriers - affected` | Permanently forbidden ([`pipeline/count_repair.py`](../pipeline/count_repair.py)) |

Note the shape of the trap: every simple predicate you would reach for
(`pred == 1`, `unaffected == 0`, `layer == figure`) is **majority-correct**.

## 5. The remaining errors

11 affected, 13 unaffected, 26 carriers. Classes, with the general lever:

| Class | Examples | Next general lever |
| --- | --- | --- |
| Status default: gold `affected=0`, predicted always 1 | 18675227 R954C/L955V, 21779290 K897T, 30246897 C566R, 29677589 R190W/R594Q | The `1/1/0` case-report default is indistinguishable by shape from these. Needs a status read, not a wider NULL |
| Identity TP against gold `0/0/0` | 29016797 N588K/S631A, 28739325 D242N, 16301357 L1825P | Gold lists the variant with no patients. Curator question, not a split mistake |
| Figure **split** still wrong | SCN5A 15671429 D1275N 23/2 vs gold 7/15 | Pedigree symbol counting (`pipeline/pedigree_extractor.py`); do not refuse every figure split |
| Symptoms ≠ diagnosis | KCNH2 25819988: one sentence holds symptoms 6, diagnosis 9, carriers 13 | Correct action is keep `carriers=13` and **null the split** — not rewrite it. See §6 |
| Abstract "a patient" vs family N | 22882672 1 vs 8, 22222782 ×3, 23917959 | Source selection, not a location rule |
| Off-by-one / partial pedigree | 24667783, 32866913, 25171853, 19406494 | Better read |

## 6. Ranked work (no one-off aliases)

1. **Do not** add per-paper or per-variant phenotype patches.
2. Keep the guard. The next paid extract inherits all four rules for free.
3. Before any verify-after-merge work (`TASKS.md` §3), settle
   `pipeline/claim_verifier.py`. `_apply_count_identity_guard` implements the
   forbidden `unaffected = carriers - affected` and is correctly skipped for
   card-aware calls (`card is not None` returns early), but
   `_apply_consistency_guards` still runs on every card and can write
   `affected = total` with `unaffected = 0`. Moving verification after merge
   without settling that can *manufacture* the pattern the copy guard then
   wipes — and on 25819988 it would land `13/13/0`, which is also wrong.
4. Pedigree / figure: count filled vs empty symbols separately; if you cannot,
   emit carriers only (already the prompt; the extractor still sometimes emits
   a fake split). Free check first: open the four error figures and the four
   exact ones and ask whether a deterministic legend separates them.
5. Deterministic `regex_table` column binder — MAE / count-recall, not this
   denominator, but the other half of "get the number right".
6. Paid gold-120 re-extract only after 3–5, and only as a measurement of those
   general rules. It will not make 15671429 correct by itself.

## 7. House constraints

- One checkout, `main`. No per-variant aliases (`ΔKPQ`, hard-coded exon maps,
  "K897T is never affected").
- Do not rewrite locked `20260813_gold120_verticalfix/report.json`.
- Do not enable `COUNT_RECOVERY_ENABLED` or un-park the per-paper final check
  to chase this number.
- Offline tests: `.venv/bin/python -m pytest tests/unit -q`.
