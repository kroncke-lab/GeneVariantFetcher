# SCN5A gold adjudications — 2026-09-08

Approved by Brett Kroncke on 2026-09-08 ("I think your proposals are good") on the
proposal in [`../README.md`](../README.md) section 1, and applied by
[`apply_adjudications.py`](apply_adjudications.py) through the repository's
`gold_v2_*` columns of `gene_variant_fetcher_gold_standard/normalized/SCN5A_recall_input.csv`.
Original values are preserved in the legacy columns; the scorers prefer the v2
values whenever `gold_v2_status` is populated, and `excluded_duplicate_current_cohort`
removes a row from every score (`utils/gold_standard.py`). No locked prediction
and no registry answer key was modified: the continuation registry's frozen
`answer_key/` copies predate these adjudications, so tranches 01-05 keep the scores
they were locked with; the four-gene recall suite and any new registry read the
adjudicated gold. Unlike the 2026-08-10 adjudications, no independent CLI audit
was run; the source rows are cited so a reader can check each decision.

## Policy applied

- A paper that prints one pooled count per nucleotide change and lists the
  contributing centres without per-centre counts has one gold row per printed
  nucleotide change; the per-centre split the legacy gold invented is not in the
  paper. Distinct nucleotide changes with the same protein effect stay separate.
- A paper that prints a per-variant phenotype split defines `unaffected` as its
  negative-phenotype column and `affected` as every phenotype-positive carrier;
  disease-specific buckets are recorded in the note.
- A variant whose phenotype the paper never reports per variant keeps its carrier
  count and takes explicit nulls for affected and unaffected; `carriers = affected`
  is not inferred from ascertainment.

## Summary

| status | rows |
| --- | ---: |
| `adjudicated_phenotype_not_reported` | 173 |
| `excluded_duplicate_current_cohort` | 69 |
| `adjudicated_variant_carrier_count` | 47 |
| `adjudicated_source_phenotype_partition` | 12 |

## Decisions with a changed value

| PMID | variant | legacy C/A/U | adjudicated C/A/U | status | source-grounded rationale |
| ---: | --- | --- | --- | --- | --- |
| 20129283 | A1680T | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5038 G>A A1680T = 2 (centres 2, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | A735V | 1/1/0 | 4/4/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2204 C>T A735V = 4 (centres 2, 4, 8, 9); the paper reports no per-centre split, so the 4 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | c.1890+5G>A | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1890 +5 G>A* = 2 (centres 2, 5); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | c.1936delC | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1936delC Q646RfsX5 = 3 (centres 2, 5, 6); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | c.2602delC | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2602delC L868X = 2 (centres 6, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | c.3840+1G>A | 1/1/0 | 6/6/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3840 +1 G>A = 6 (centres 1, 3, 4); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | c.4437+5G>A | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4437 +5 G>A* = 2 (centres 3, 5); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | D1243N | 1/1/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3727 G>A D1243N = 5 (centres 1, 2, 5); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | D1275N | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3823 G>A D1275N = 3 (centres 1, 5); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | E1225K | 1/1/0 | 4/4/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3673 G>A E1225K = 4 (centres 1, 5, 6, 7); the paper reports no per-centre split, so the 4 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | E161K | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 481 G>A E161K = 3 (centres 3, 4); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | E746K | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2236 G>A E746K = 3 (centres 1, 2, 7); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | G1319V | 1/1/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3956 G>T G1319V = 5 (centres 2, 3, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | G1408R | 1/1/0 | 7/7/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4222 G>A G1408R = 7 (centres 1, 4, 5, 7); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | G1743E | 1/1/0 | 6/6/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5228 G>A G1743E = 6 (centres 2, 3); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | G1743R | 1/1/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5227 G>A G1743R = 5 (centres 4, 5, 7, 9); the paper reports no per-centre split, so the 4 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | G752R | 1/1/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2254 G>A G752R = 5 (centres 1, 5); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | I1660V | 4/4/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4978 A>G I1660V = 5 (centres 2, 3, 5, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | K1493X | 2/2/0 | 1/1/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4477 A>T K1493X = 1 (centres 2); gold 2 conflated K1493del (a separate row). approved by Brett Kroncke 2026-09-08. |
| 20129283 | L1393X | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4178 T>A L1393X = 3 (centres 1, 3, 9); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | L868X | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2602delC L868X = 2 (centres 6, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | N927S | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2780 A>G N927S = 3 (centres 3, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | P.K1493DEL | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4477_4479delAAG K1493del = 2 (centres 1, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | P336L | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1007 C>T P336L = 2 (centres 2, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | Q1695X | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5083 C>T Q1695X = 2 (centres 1, 4); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R104Q | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 311 G>A R104Q = 3 (centres 1, 7, 8); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R104W | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 310 C>T R104W = 2 (centres 1, 2); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R121Q | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 362 G>A R121Q = 2 (centres 2, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R1232W | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 3694 C>T R1232W = 3 (centres 1, 2, 9); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R1583C | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4747 C>T R1583C = 2 (centres 1, 2); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R1623X | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4867 C>T R1623X = 2 (centres 2, 3); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R1638X | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4912 C>T R1638X = 3 (centres 2, 3); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R225W | 2/2/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 673 C>T R225W = 3 (centres 1, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R367H | 1/1/0 | 6/6/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1100 G>A R367H = 6 (centres 1, 2, 8, 9); the paper reports no per-centre split, so the 4 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R376H | 1/1/0 | 4/4/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1127 G>A R376H = 4 (centres 3, 4, 8); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R526H | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1577 G>A R526H = 2 (centres 1, 5); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R535X | 1/1/0 | 4/4/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1603 C>T R535X = 4 (centres 1, 2, 4, 5); the paper reports no per-centre split, so the 4 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R878H | 1/1/0 | 5/5/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2633 G>A R878H = 5 (centres 1, 2, 4, 5, 7); the paper reports no per-centre split, so the 5 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R893H | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2678 G>A R893H = 3 (centres 1, 3, 4); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | R965C | 1/1/0 | 3/3/0 | `adjudicated_variant_carrier_count` | Table 4 prints 2893 C>T R965C = 3 (centres 2, 4, 5); the paper reports no per-centre split, so the 3 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | S1672Y | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5015 C>A S1672Y = 2 (centres 1, 4); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | T1620M | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4859 C>T T1620M = 2 (centres 2, 9); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | T1709M | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 5126 C>T T1709M = 2 (centres 1, 8); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | T220I | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 659 C>T T220I = 2 (centres 2, 3); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | T632M | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 1895 C>T T632M = 2 (centres 2, 4); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | V1405M | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 4213 G>A V1405M = 2 (centres 1, 7); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 20129283 | V232I | 1/1/0 | 2/2/0 | `adjudicated_variant_carrier_count` | Table 4 prints 694 G>A V232I = 2 (centres 2, 6); the paper reports no per-centre split, so the 2 one-row-per-centre gold rows are consolidated here. approved by Brett Kroncke 2026-09-08. |
| 30059973 | D1790G | 7/7/0 | 7/2/5 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2590): ECG-negative 5, isolated LQT3 0, BrS 0, PCCD 1, SSS 1, DCM 0, overlap 0; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | D356N | 3/3/0 | 3/2/1 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2591): ECG-negative 1, isolated LQT3 0, BrS 0, PCCD 2, SSS 0, DCM 0, overlap 0; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | E1784K | 69/69/0 | 69/40/29 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2577): ECG-negative 29, isolated LQT3 13, BrS 0, PCCD 17, SSS 0, DCM 0, overlap 10; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | E901K | 7/7/0 | 7/4/3 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2594): ECG-negative 3, isolated LQT3 0, BrS 0, PCCD 1, SSS 0, DCM 0, overlap 3; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | G1319V | 8/8/0 | 8/2/6 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2586): ECG-negative 6, isolated LQT3 0, BrS 0, PCCD 2, SSS 0, DCM 0, overlap 0; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | G1743E | 8/8/0 | 8/2/6 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2579): ECG-negative 6, isolated LQT3 0, BrS 0, PCCD 1, SSS 0, DCM 0, overlap 1; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | I1768V | 9/9/0 | 9/3/6 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2578): ECG-negative 6, isolated LQT3 1, BrS 0, PCCD 1, SSS 0, DCM 0, overlap 1; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | P.Q1507_P1509DEL | 9/9/0 | 9/5/4 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2585): ECG-negative 4, isolated LQT3 5, BrS 0, PCCD 0, SSS 0, DCM 0, overlap 0; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | P.Q646RFSX5 | 7/7/0 | 7/7/0 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2592): ECG-negative 0, isolated LQT3 0, BrS 0, PCCD 5, SSS 0, DCM 0, overlap 2; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | P.Y1795_E1796INSD | 7/7/0 | 7/3/4 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2588): ECG-negative 4, isolated LQT3 0, BrS 1, PCCD 1, SSS 0, DCM 0, overlap 1; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | V1763M | 6/6/0 | 6/4/2 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2595): ECG-negative 2, isolated LQT3 3, BrS 0, PCCD 0, SSS 1, DCM 0, overlap 0; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |
| 30059973 | V411M | 10/10/0 | 10/8/2 | `adjudicated_source_phenotype_partition` | Table 14 (source line 2583): ECG-negative 2, isolated LQT3 5, BrS 0, PCCD 2, SSS 0, DCM 0, overlap 1; affected = every phenotype-positive carrier. approved by Brett Kroncke 2026-09-08. |

## Excluded per-centre duplicates (69 rows, PMID 20129283)

Each is a legacy `1 / 1 / 0` (or partial) row for a variant whose consolidated row above now carries the printed Table 4 count:

A1680T, A735V, A735V, A735V, D1243N, D1243N, D1275N, E1225K, E1225K, E1225K, E161K, E746K, E746K, G1319V, G1408R, G1408R, G1743E, G1743R, G1743R, G1743R, G752R, I1660V, L1393X, L1393X, L868X, N927S, P.K1493DEL, P336L, Q1695X, R104Q, R104Q, R104W, R121Q, R1232W, R1232W, R1583C, R1623X, R1638X, R225W, R367H, R367H, R367H, R376H, R376H, R526H, R535X, R535X, R535X, R878H, R878H, R878H, R878H, R893H, R893H, R965C, R965C, S1672Y, T1620M, T1709M, T220I, T632M, V1405M, V232I, c.1890+5G>A, c.1936delC, c.1936delC, c.2602delC, c.3840+1G>A, c.4437+5G>A

## Explicit-null phenotype (173 rows, PMID 30059973)

Table 5 of Baruteau 2018 lists per-variant occurrences with no phenotype; Table 1 reports 196 of the 442 carriers with a negative ECG phenotype. Carriers stand; affected and unaffected are explicit nulls for:

A1186T, A1221V, A1330D, A1330T, A1428S, A1746T, A242D, A332T, A413T, A425T, A735V, c.3840+1G>A, c.393-1C>T, c.4245+1G>T, c.4437+5G>A, c.611+1G>A, c.703+1G>A, c.934+1G>A, C341Y, C683R, D1275N, D1595N, D1790N, E1053K, E1105X, E1107X, E1225K, E1240Q, E161K, E346X, F1210S, F1460L, F1486L, F355I, F4V, F892I, F93S, G1262S, G1369X, G1406R, G1408R, G1481E, G1481R, G1481V, G1631D, G1743R, G514C, G897E, H1200Y, H1849R, I1334V, I1521K, I1524T, I1660V, I1749N, K1504E, L1222R, L1346P, L1373X, L1501V, L1634P, L1786R, L212P, L224F, L227P, L276P, L736P, L839P, M1487L, M1498T, M1766L, M1793K, M369K, N1325S, N1380K, N1722D, N1774D, N1774H, N406K, N927S, P.E1071GFSX76, P.E1107RFSX24, P.E1165RFSX6, P.F1486DEL, P.F1617DEL, P.F1775LFSX15, P.F851CFSX19, P.F861WFSX90, P.G1015DFSX14, G1031fsX, P.G1748DEL, P.I1570DUP, P.I1758DEL, P.I759FFSX6, P.L1222LFSX7, P.L1302VFS18, P.L1339DEL, P.L729DEL, P.N1380DEL, P.N291TFSX52, P.N841TFSX2, P.P2006LFSX32, P.Q1475NFSX6, P.Q90WFSX14, P.T631VFSX101, P.Y774TFSX28, P1332L, P1332Q, P336L, P717L, Q1000X, Q1059X, Q1118X, Q1475L, Q1491H, Q1695X, Q779K, Q779X, R1023C, R121Q, R1232W, R1583H, R1623Q, R1623X, R1626C, R1632L, R1644H, R1897W, R1944X, R1991Q, R222Q, R225W, R282C, R340W, R367C, R367H, R376C, R43X, R535X, R689C, R814Q, R878C, R893C, R965C, S1074G, S1079F, S1079Y, S1672Y, S1710L, S401P, S941F, S941N, T1209R, T1304M, T370M, T630M, V1747M, V1763L, V1777L, V1777M, V240M, V263I, W1191X, W1345C, W156X, W374G, Y1449C, Y1615X, Y1767C, Y1795C, Y1795N, Y1811N, Y1811X
