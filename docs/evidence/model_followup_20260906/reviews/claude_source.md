# Audit: Source-Adjudicated Patient-Count Contract

## MYB204 (48 rows, Additional file 1)
- Row count checks out: 48 rowspan-expanded entries. Endpoint split = 37 positive (numeric age)/10 No Dx (negative)/1 "¿" (unknown) — consistent with the 48-person roster.
- **Real trap — contradictory rowspan metadata**: H40 (IV-1) correctly shows only column 0 (Mutation) inherited, Family blank — matches contract. But **H41 (IV-5)** lists inherited columns `[0,1]`, implying Family *was* inherited via real rowspan, while the literal cell is still blank. This contradicts the instruction that both IV-1/IV-5 family cells are blank-not-merged. Anyone trusting the inherited-columns flag for H41 would wrongly treat the H49 family link as grid-native rather than prose-supplied — flag and override with the blank literal value, keep the H49 join sourced to prose only.
- **Prose-only person excluded correctly, but scope mismatch worth flagging**: prose cites "H197-III:4 = R502Q" (suggestive ECG, no HCM), but no H197 family appears anywhere in the 48-row table (only H147/H614 carry R502Q). Contract rightly excludes this person, but the prose's "46 relatives / 24 carriers" framing is a *different, larger* study population than the 48-row FU table — don't reconcile these totals as if same scope.
- No Dx ≠ ECG-normal is empirically confirmed: H1 (No Dx) still shows ECG "LBBB" — a real case where negative-endpoint status coexists with an abnormal ECG finding.

## RYR304 (15 numbered subjects)
- 13 positive/2 negative (subjects #6, #12) verified against Symptoms column and prose exceptions — correct, and correctly scoped as VT/SCA-at-presentation, not global CPVT-affected status.
- Family correctly nulled despite narrative relational language (sibling/cousin of subject 1, subject 5's self-referential "Sister of subject 5" typo). This is a genuine source error; don't use it to infer or repair family linkage.

## RYR254 (8 carriers)
- 7 positive/1 negative (II:9, footnote d: VE alone not diagnostic) verified.
- 3 prior-symptomatic (II:1, II:9, III:9) is a separate historical field, correctly not conflated with the EST diagnostic endpoint.
- **Real trap — narrative/table contradiction**: prose (N1) claims "all patients were asymptomatic" at follow-up, but the table's own "Events during evolution" column lists **Presyncope for the proband (II:1)**. The packet's "all asymptomatic follow-up" framing glosses over this internal source contradiction; the table-literal exception should be preserved, not suppressed to match the summary sentence.

## MYB213 (22 listed patient IDs, Panel a)
- Count of 22 confirmed (footnote c/d suffixes correctly treated as relatedness markers, not extra people: 180c/186c/208c = 3 people, 177d/185d = 2 people).
- 22 listed vs. 19 unrelated reconciles: the 3-person gap is explained by explicitly-stated relatives (D605H "father, son, and nephew"; L1084P pair) — footnote markers are the mechanism for including relatives beyond unrelated probands, per contract.
- **Provenance trap, not silently fixed**: row for c.2827C>T lists protein change "R493X" — inconsistent with standard MYBPC3 numbering for this variant, but per contract this is copied as reported. Downstream consumers must not assume HGVS-correctness of source-literal protein notation.
- All 22 individuals are HCM-positive by cohort-enrollment definition — this packet has no unaffected/asymptomatic roster members, unlike MYB204; don't cross-apply MYB204's mixed affected/unaffected pattern here.

## Remaining acceptance/reconciliation traps
1. Rowspan-metadata vs. literal-cell conflicts (MYB204 H41) must be resolved in favor of literal blanks.
2. Packet-scoped totals must never be compared to whole-paper prose totals (MYB204 46 vs 48; RYR304 13/2 vs global CPVT).
3. Narrative "all-clear" summary claims (RYR254 follow-up) can contradict table-literal exceptions — table wins.
4. Source-reported protein notation errors (MYB213 R493X) must be preserved, not corrected, and flagged for downstream QA.
5. None of these source-validated counts authorize production promotion — audit trail only.
