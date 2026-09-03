# Mixed-gold tranche 01: baseline, root cause, candidate (2026-09-03)

First paired measurement of the 29-tranche mixed-gold suite
(`benchmarks/evaluation_tiers/mixed_gold/`). Both arms extract the same 49
gene-paper attempts (45 PMIDs; KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 23)
from the frozen corpus, gold-free, and are scored on the paper-derived lane
only. Consumption is recorded in `mixed_gold/consumption_log.jsonl`.

## Baseline arm

`runs/20260903_protocol_mixed01_baseline`, frozen protocol commit `506a949c`
(runtime fingerprint `a7de3089…`, 238 files), 5 genes in parallel, 25.6 min.

| lane | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| paper-derived (primary) | 155 | 61 | 87 | 71.8% | 64.0% | 67.7% |
| linkage-assisted (diagnostic) | 209 | 120 | 33 | 63.5% | 86.4% | 73.2% |

Counts (paper-derived): carriers supplied 48 / 242 gold assertions (19.8%),
conditional MAE 0.812; affected 101 / 242 (41.7%), MAE 0.079; unaffected
7 / 242 (2.9%). By gene: KCNH2 17/11/8, KCNQ1 10/9/56, MYBPC3 73/1/1,
RYR2 21/5/2, SCN5A 34/35/20 (TP/FP/FN).

### Where the 87 misses went (`diagnostics/fn_root_cause/`)

`scripts/recall_audit/fn_root_cause.py` walks each missed gold row through
the stages the run left on disk: corpus (source-presence sweep class) → the
text the run staged → the LLM request → the LLM response → the extraction
JSON → the DB row and its source layers → the paper-derived lane.

| leaf | rows | what it means |
|---|---:|---|
| `acquisition` | 70 | absent from every byte on disk: KCNQ1 14678125 (41 rows, tables never captured), the three gene arms of 15176425 (18; garbled PDF text / 3 KB stub / no source), KCNQ1 21131640 (9), SCN5A 12566525 (4), two abstract stubs |
| `unknown_notation` | 12 | indel/frameshift/structural notation the probe cannot search, or figures on disk |
| `parser_dropped` | 1 | SCN5A 16764707 `L860fsX`: model returned it as `source_notation` only; filter dropped it as nameless |
| `condensing` | 1 | MYBPC3 33673806 `p.Glu843_Arg845del`: deterministic route did not recognise the `HGVSp` column |
| `paper_row_lost_to_linkage_origin` | 1 | SCN5A 28469501 `G289S`: extracted from the paper, but notes mentioning ClinVar made its origin `clinvar` |
| `matcher` | 2 | RYR2 17556193 `R2267H`, `S4565R`: gold lists each variant twice; one-to-one pairing leaves the duplicate unmatched |

Two of the twelve "unknown notation" rows were also reachable: SCN5A 16764707
`1795insD` and `1570insI` were in the model's response as `source_notation`
only and were dropped by the same nameless-row filter.

Counts and extras: MYBPC3 33673806's `Case count` column produced 74 affected
values and zero carrier values (gold carries the same integers as carriers);
SCN5A 23414114's ESP genotype-frequency table produced 21 counted extras;
KCNH2 15851119's 9 extras are real table rows the gold does not list.

Tranche 01 is unusually acquisition-bound: the source-presence sweep puts 75
of its 242 gold rows (31%) behind the hard ceiling versus 15.8% suite-wide
(`gold_source_presence_sweep_20260903.md`). The reachable pool for any
reading-protocol change on this tranche is therefore about 13 rows.

## Candidate arm (commit `b56f469f`)

Six gene-agnostic fixes, each traced to a stage above:

1. Publisher access shells whose folded supplements carry the tables are
   usable/reusable corpus sources (RYR2 22222782: 235 supplement table rows
   read from a PubMed abstract in the baseline; 22 of 70 refused corpus full
   texts suite-wide).
2. Corpus reuse prefers a `CLEANED` sibling ≥ 1.5× a worse `FULL_CONTEXT`;
   preprocessing keeps a staged richer rendering.
3. Source-only legacy protein spellings become identities (`L860fsx89` →
   `L860fsX`; `1795insD`, `1570insI` kept as legacy); the projection falls
   back to `legacy_notation`.
4. Router recognises `HGVSp`/`HGVSc` headers and `p.(Glu101Lys)` values.
5. A closed case-series count column sets the carrier total equal to the
   counted cases when there is no carrier and no unaffected/control column.
6. ClinVar/PubTator are a row's origin only through the ingest markers;
   notes mentioning ClinVar no longer move a paper row to the linkage lane.
   The trust gate also types `genotype` count labels as `population_count`.

### Result (`runs/20260903_protocol_mixed01_candidate`, locked 2026-09-03)

| lane / metric | baseline | candidate | delta |
|---|---:|---:|---:|
| paper-derived TP / FP / FN | 155 / 61 / 87 | 157 / 54 / 85 | +2 / −7 / −2 |
| paper-derived recall | 64.05% | 64.88% | **+0.83 pp** (one-sided 95% LB 0.00) |
| paper-derived precision | 71.76% | 74.41% | +2.65 pp (LB +0.16) |
| counted-extra precision | 82.0% | 91.3% | +9.3 pp (34 → 15 counted extras) |
| carriers supplied / 242 | 48 | 125 | +77 (non-zero gold: 44 → 122) |
| carrier conditional MAE / RMSE | 0.812 / 3.96 | 0.104 / 0.48 | |
| affected supplied / 242 | 101 | 81 | −20 |
| unaffected supplied / 242 | 7 | 4 | −3 |
| total tokens | 1.079 M | 1.091 M | +1.1% |

`run_eval.py compare --phase discovery` → **`reject_or_revise_candidate`**:
the precision non-inferiority guard passed, the recall rule did not (observed
+0.83 pp < the preregistered +1.0 pp). The registered rule scores identity
only; the count changes above are outside it and are reported here because
they are the campaign's second objective.

What moved, by paper:

- SCN5A 16764707: `L860fsX` recovered (the model emitted `p.Leu860fs*89`
  this time; fix 3 was not needed for it). `1795insD` / `1570insI` were not
  in the model's response at all in this run — the baseline's emission of
  them as `source_notation` was stochastic.
- SCN5A 28469501: `G289S` recovered by fix 6 (paper row no longer attributed
  to ClinVar).
- RYR2 22222782: read from the 51,680-byte shell-plus-supplements source
  (fix 1) instead of a 2,785-byte abstract: 7/7 identities as before, 4 → 0
  extras, 5 carrier values supplied.
- KCNQ1 23092362: 3 → 0 extras (a run-to-run difference, not attributable to
  a fix).
- MYBPC3 33673806: the `HGVSp` column is now read (fix 4) and the in-frame
  deletion is a TP, but `c.3408C>A` became the paper's one FN: the scorer
  paired the prediction `p.Tyr1136* c.3408C>A` with the gold **deletion**
  `c.3407_3409delACT` at the same codon and left the true gold row unmatched.
  Carriers are supplied on 73/74 rows (fix 5) with MAE 0.027, but affected
  fell 74 → 55: the always-on phenotype guard clears an affected value equal
  to a multi-carrier total unless its count type is `case`, and the
  case-series rule left the type generic (`copied_carriers_onto_affected`).

Reviewer input on the candidate (Grok 4.6 `xhigh`, read-only; Gemini 3.1 Pro
`high`, brief-only): rank the case-series rule and the fold-marker gate as the
two regression risks. Concretely: a bare case-control `Cases` column would
map to affected and, once typed `case`, survive as a per-variant count; the
fold-marker gate accepts ≥ 1,000 characters of fold chrome without any table;
`p.(?)`/`p.(=)` could be unwrapped into identities; and IUPAC nucleotide
letters outside A/C/G/T (`1100delN`) would promote as protein legacy. All
four are addressed in the v2 candidate.

### Root cause of the candidate's remaining misses

Unchanged in kind: 70 acquisition, 12 notation-unknown, the two RYR2
17556193 gold duplicates, and one scorer pairing error (MYBPC3 above).
