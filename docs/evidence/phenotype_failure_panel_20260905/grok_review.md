**Challenge.** The 22-paper / 30-attempt set is a ranked, high-discrepancy slice of 336 already-locked gene/run attempts, not a new benchmark. It is evidence of mixed failure modes, not a warrant for one count architecture.

The prior “unparsed supplements / stale pedigree” story is false on these papers: 30059973’s 70-page manuscript including supplements was already flat PDF text; 26496715’s Extanded tables.doc is already 382 markdown table lines; 33673806 additional file 2 is already converted with explicit `c.3181C>T` 24 and `c.2373dup` 14. 25814417 already has a source-backed patient-row implementation and is still the largest supplied disagreement (gold 73/106 vs source 97/62 +26 unknown, 68 units).

What the inspections actually show is three different problems that the proposal then folds together.

**1. Gold A/U is not one endpoint.** SCN5A 30059973: identity is strong (182 TP / 3 FN / 3 FP), phenotype fields are empty, 445 A+U discrepancy units. Supplemental Table 5 is carrier N; Table 11 is presentation/follow-up events (E1784K: 0 arrests, 9 syncope, 60 asymptomatic at presentation, 10 later cardiac events); Table 14 is ECG/diagnosis overlap (29 negative ECG, 13 LQT3, 17 PCCD, 10 overlap among 69 E1784K carriers). Gold is affected=69, unaffected=0. Those 69 are the carrier denominator, not symptoms and not ECG. Copying N→affected would fit this gold row and be wrong on Tables 11/14. KCNH2 25819988: 13 carriers, 6 symptomatic, 9 LQTS diagnosed; pred 6/7 vs gold 9/4. Gold is closer to diagnosis than symptoms; residual subtraction (13−6=7) is forbidden. MYBPC3 33673806: suspected HCM, presumed affected, diagnoses not verified; singletons got A, multi-person A stayed NULL. RYR2 33315912 is 17/42 (+3 unknown) vs gold 15/47. Gold zero often means unreported, not healthy. An atomic-observation extractor that later collapses to A/U will re-create the error unless the collapse is refused.

**2. “On disk” ≠ usable source, and corpus growth ≠ table recovery.** KCNQ1 14678125: gold118 35/6/9 vs mixed01 0/41/0; both runs had a 3.3k abstract; current 13.8k object is a publisher shell (abstract + references + navigation), no mutation table, 289 discrepancy units. Legacy TPs are mostly linkage, not table reading. SCN5A 20031634 is still abstract-only (13 FN, 108 units); a web index of Table 1 (G1408R 14=4/9/1; total 115=54/55/6 vs gold 60/48) is not a validated PDF (host timed out). SCN5A 25163546: UCL manuscript PDF is new and verified; supplements unlocated; that is retrieval, not count recovery. 26496715 still lacks the body that defines ascertainment, so the converted “No. of patients” column is not an A/U endpoint. Deterministic replay cannot precede acquisition on this subset. “No whole-paper rerun” is only valid where the locked artifact already is the paper.

**3. Joins and identity are not the same job as counts.** MYBPC3 21302287 Table 3b is patient IDs (`c.2340+18C>G` → 278, 251, 268, 194), diagnosed HCM — 278 is an ID, not a count, and the table is not missense-only. RYR2 30403697 Table 1 has two variant columns, 15 patient rows, multiple patients per family, an explicit asymptomatic row, and phenotypically silent parents in prose; family or figure-symbol totals are the wrong aggregate. 19632626: 63 identity FPs plus S209P 6/1 vs gold 0/6 under AF ascertainment. 32533946: 44 FPs, no clinical counts (engineered/reference vs gold-incomplete). Exact controls (W4645R 4/2/2 from an explicit abstract; C2277R 8/7/1 from patient rows) only prove the scorer can accept a fully specified triple. They do not validate atomic grouping on these tables.

The proposed manifest + PDF coordinates + atomic observations is the right *shape* for 21302287 / 30403697 / 30059973 table structure, and the right *refusals* (no NULL-means-healthy, no symptom⇒affected synonym, no residual subtraction, no force-closed A+U partition). It is not justified as the next global pipeline. Measuring acquisition vs capture separately is required; treating converted supplements as the next extraction target is not.

**Top 3 implementation steps**

1. **Triage each of the 30 attempts into stored classes before writing extractors:** (A) no authoritative body/table/supplement; (B) component present (PDF text, converted .doc, additional file) but not consumed as labelled tables; (C) source consumed and A/U still disagree. Only (C) is a count-logic bug. 14678125, 20031634, 25163546 are (A). Do not score counts on (A).

2. **Lock components, not papers.** Manifest body / table / figure / supplement / rendered block with identity. Gate ingest on visual or structural validation of that component (30059973 table images yes; 20031634 snippet no). Preserve page coordinates and orientation only where the artifact is a PDF page. Converted markdown has no such coordinates.

3. **Emit endpoint-labelled atoms and stop.** Fields: gene, variant, subject ID, genotype, clinical status, timepoint, location. Group distinct people. Keep carrier N, symptoms, ECG class, diagnosis, unknown as separate series. Never map them to gold A/U in the extractor. Adjudicate the map per paper, or leave unmatched.

**Regression tests (must stay red if violated)**

- Controls: W4645R 4/2/2; C2277R 8/7/1.
- 30059973: do not emit A=69 from carrier N without an endpoint label; do not synonym Table 11/14 into A/U; 0 arrests is not U.
- 21302287: `278, 251 268, 194` → four IDs, never count 278.
- 30403697: two variant columns; asymptomatic row ≠ U from family size.
- 25819988: 13 / 6 / 9 remain three numbers; pred U must not become 13−6.
- 33673806: no blanket promote of suspected-HCM N; singleton vs multi-person A policy must be one rule.
- 14678125: 13.8k shell scores as no-table, not as recovered source.
- 20031634: indexed 4/9/1 must not ingest until a validated PDF/HTML exists.
- 26496715: “unparsed supplement” is a false alarm (382 lines exist).
- 25814417 / 33315912: unknown 26 / 3 stay unknown; do not force gold partitions.

**Must remain uncertain**

Whether gold A/U is the same endpoint across even these 22 papers (the inspections say it is not). Whether this slice represents the other 306 attempts. Whether flat PDF / converted markdown retains joinable columns. True 20031634 table values; 25163546 supplement contents; where 14678125’s mutation table lives. Whether 25814417, 33315912, 19632626 S209P, and 32533946 are extractor error or gold error. Cost of atoms-then-group vs per-paper adjudication. No phenotype-precision movement is claimed.
