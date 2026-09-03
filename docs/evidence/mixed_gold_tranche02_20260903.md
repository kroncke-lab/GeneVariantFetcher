# Mixed-gold tranche 02: v2 discovery (2026-09-03)

Second paired measurement of the mixed-gold suite. Tranche 01 rejected the v1
candidate on the recall rule (+0.83 pp against a +1.0 pp bar; see
`mixed_gold_tranche01_20260903.md`), so under the registry's
`reject_or_revise_candidate` outcome tranche 02 is a **new discovery** for a
revised candidate (v2), not a confirmation. 49 gene-paper attempts (BRCA1 1,
KCNH2 8, KCNQ1 10, RYR2 6, SCN5A 24), 303 gold rows.

## What tranche 02 can show (gold-blind)

From the source-presence sweep (`gold_source_presence_sweep_20260903.md`):
232 of 303 gold rows are present in the body text, 24 sit behind the hard
acquisition ceiling (11 no source, 12 stub, 1 absent substitution) and 47 more
are undecidable (28 non-searchable notation, 19 figure-only). One paper,
SCN5A 30059973 (a manuscript PDF whose supplemental Table 5 is a fixed-width
text table with the HGVS dots stripped: `c5972G>A`, `pArg1991Gln`), carries
185 of the 303 gold rows; the existing fixed-width parser emits 185 rows from
it with 177 dotted protein identities, so the tranche is dominated by how that
one paper's counts and identities survive the downstream stages.

Gold-blind pre-check of the v1 fixes against tranche 02 source files: the
folded-supplement shell rule fires on one attempt (SCN5A 25341504), the
richer-CLEANED rule on one (SCN5A 25650408, a 0-byte FULL_CONTEXT beside a
5.2 KB CLEANED), the frameshift-with-stop-distance spelling appears in four
attempts (three gene arms of 29544605 and KCNQ1 33498651), and one BRCA1 table
carries a case-series count header. Five attempts are abstract-only stubs and
one has no source. The v1 fixes are therefore expected to be nearly inert on
identity here; v2's additions are hardenings, not new recall levers.

## Baseline arm

`runs/20260903_protocol_mixed02_baseline`: the frozen `506a949c` reading
protocol, re-materialised by checking out the nine behaviour-bearing runtime
files from that commit (`db_to_predictions.py`, `content_validation.py`,
`extraction.py`, `source_quality.py`, `steps.py`, `table_router.py`,
`trust_gate.py`, `legacy_notation.py`, `source_layers.py`). Its runtime
fingerprint (`8235a53c…`, 240 files) differs from tranche 01's baseline
(`a7de3089…`, 238) only by the two audit scripts added since
(`gold_source_presence_sweep.py`, `fn_root_cause.py`), which the pipeline does
not import, and by ruff's reformatting of files that were untracked at setup.

Result (`runs/20260903_protocol_mixed02_baseline`, locked 2026-09-03):

| lane | TP | FP | FN | precision | recall |
|---|---:|---:|---:|---:|---:|
| paper-derived (primary) | 268 | 54 | 35 | 83.2% | 88.4% |
| linkage-assisted (diagnostic) | 278 | 150 | 25 | 65.0% | 91.7% |

Counts: carriers supplied 232 / 303 (non-zero gold 231 / 285), conditional
MAE 0.026; affected 35 / 303; unaffected 6 / 303; counted-extra precision
97.5% (7 counted extras). By gene: BRCA1 0/18/1, KCNH2 14/8/2, KCNQ1 18/22/0,
RYR2 23/0/0, SCN5A 213/6/32. SCN5A 30059973 alone: 182 TP / 3 FP / 3 FN with
carriers on every matched row.

Where the 35 misses went: SCN5A 12106943 (16; figure images only),
SCN5A 25650408 (12; 5 KB stub), SCN5A 30059973 (3), BRCA1 18627636 (1, with
18 extras against a one-row non-exhaustive gold), KCNH2 26022375 and 14642688
(1 each), SCN5A 21167004 (1). The reachable pool for a reading change is at
most ~5 rows (1.7 pp), so the +1.0 pp discovery bar is arithmetically within
reach only if nearly every reachable row is recovered.

## Candidate arm (v2)

v1's six fixes plus:

1. Case-series rule types the affected count as `case` so the phenotype guard
   keeps it, and requires a count word beside the people token so a bare
   case-control `Cases` column never fires (Grok/Gemini counterexample).
2. The folded-supplement gate requires table rows or variant tokens in the
   fold, not merely 1,000 characters.
3. Source-only legacy promotion refuses payloads spelled only in A/C/G/T/N.
4. Predicted-protein unwrap accepts only `p.(<amino acid><position>…)`; never
   `p.(?)` or `p.(=)`.
5. Trust gate: only `genotype frequency/count/distribution` labels are
   population counts.
6. Scorer: when both strings carry concrete but different cDNA alleles, the
   codon-position bridges do not fire (tranche 01's MYBPC3 pairing error).

Result: pending.
