REVIEW_VERDICT: reject-as-specified. Salvage a narrower reuse of the existing closed-case-series + control-stamp path. Do not ship tier 2. Do not spend the $100 two-arm extract.

The PI asked for affected **precision**. This plan is a **supply/recall** intervention. Conditional exactness is already 92%; you are proposing to copy carriers onto affected. The contract's standing rule is "a wrong number is worse than a missing one." The 477-row pile is mostly a **column-role miss**, not a missing derivation module.

---

**Q1.** Draw the line at the **count column's own noun**, not at paper-level enrollment.

Source-backed (ship): the only per-variant count header is a case-series people header (`No. of patients`, `Probands (n)`, `Case count`, `No. of Subjects`, and yes `No. of unrelatedindividuals`). That is a reading of a labeled count. `table_router._is_case_series_people_header` already implements it, and `_router_count_provenance(..., case_series=True)` already stamps `count_type="case"` so the guard keeps `affected==carriers`. It does **not** fire on the 477 because those columns were bound as **carriers** (`count_idx` set, `aff_idx` None). Kapplinger Table 4's header is `No. of unrelatedindividuals` (wrapped). That string already satisfies the existing predicate (`individual` + `No. of`). Mapping it to the affected role would recover the 239 without a new module.

Manufactured (do not ship): caption names BrS/LQT/CPVT **and** a generic `Number`/`COUNT`/`N`/`individuals` column is treated as affected. That is exactly "enrollment in a disease cohort does not by itself prove affected." Kapplinger is the type specimen: methods accept "possible or definite BrS", "ECG was not always available", and Table 4 is "possible BrS1-associated mutations" in people "referred for BrS genetic testing." Gold's N/N/0 is a lab convention, not a phenotype column.

The frequency-noun rule is **not** a license here. It converts `freq × named within-study N` and then reads the noun. It does not say "any people-count in a disease paper is affected."

The one existing exception is audited patient-row aggregation of **explicit phenotype cells**, plus claim-verifier's closed quote "these N carriers are N target-disease patients." Your module is `affected := carriers` with a caption stamp. That is the pattern `phenotype_count_guard` was written to kill.

**Tier 2: do not ship.** Title/abstract ascertainment plus a garbled `Coding Effect` parenthetical is paper-level NLP. 52 rows (28341781) are not worth a classifier that will fire on every "Probands With X" registry, molecular-autopsy series, and panel-testing paper. Fix the parenthetical `(n)` parser locally or leave them null.

---

**Q2.** Shapes that fool the stated classifier

| Shape | Why it fires | What is true |
|---|---|---|
| Case-control `Variant \| Cases \| Controls` with only `Cases` bound | Path A (`cases`+count word) | Remainder are controls; `affected:=Cases` can be OK **only if** Controls is seen and exclude-on-two-count-columns fires. You did not write that. |
| `Clinical characteristics of patients with KCNH2 mutations` | Path A/B | `n` is the characterized subset, often mixed index+relatives, often symptoms ≠ diagnosis (contract: "N patients carried… is not enough when a smaller symptom subset is reported"). |
| `Rare variants in 5,000 consecutive patients referred for genetic testing` | `patients` + count | Referral ≠ affected. Kapplinger-class. |
| SIDS/SUDY `Variants in 42 autopsy cases`, column `n` | Path A does **not** require the run's disease | Dead infants, mixed diagnoses, not CPVT/BrS affected. |
| BRCA clinic `BRCA1 mutations in 1,000 patients`, `n` | Path A | Genetics-clinic series mix affected probands and tested relatives unless the caption says so. |
| BMPR2 PAH registry `mutations in 80 patients and relatives` if "relatives" is only in methods | Path B | Cascade contamination. |
| MYBPC3 `Case count` already handled; `Number` in a sarcomere catalogue with a disease caption | Path B | Gene-elusive / VUS / relatives. |
| `Allele count` headed `Number` under `2,611 BrS patients` | Path B if `allele` is only in prose | Allele ≠ person. |
| Compilations `Published SCN5A mutations (n reports)` | `n` + disease | Other studies. Attribution check is parked. |
| `n observations` / `occurrences` / `families (n)` as `N` | `N` is in your people list | Contract already forbids occurrences and family-as-people. |
| 28237968 Table 2 | You exclude on `relative`/`symptomatic` in caption — good, must be hard-fail | Implicit-one-carrier rows with a CPVT title would be **tier 2** if caption lookup fails and stores `Table 2`. |

Missing exclusions: `putative`, `possible`/`suspected`/`referred`, `clinic-based`/`consecutive patients`, `autopsy`/`SIDS`/`SUDY`/`post-mortem`, `literature`/`published`/`reported previously`/`this study and`, `occurrence(s)`, `chromosomes`, `hom`/`het`/`genotype count`, `founder`, `yield`/`prevalence`/`positive rate`, `VUS`/`secondary findings`/`incidental`, `panel`/`any indication`, `index and relative` as a status column, `testing center` (this is the 45-row gold split), case **and** control columns together, `characteristics of`.

Too aggressive (will eat the 477):
- **`genotype-positive`**: Kapplinger's introducing sentence is "identified in the 438 genotype-positive cases (Table 4". Caption lookup that searches on `Table 4` **excludes the 239-row prize**.
- **`frequency`**: 27566755 caption is "…and Frequency in 406 LQT3 Patients". That is 47 rows.
- **`allele`**: Table 2 caption is "Control variants found in 2,600 reference alleles". If EXCLUDE runs before CONTROL, you never stamp the 55 control zeros you want.
- **`carriers`**: real headers `No. of mutation carriers` in disease series.
- **`phenotype`**: 28341781 title is "Genotype-Phenotype Correlation…". Kills the paper that needs tier 2.
- **`affected`/`unaffected` as bare substrings** in a 40-line caption window: too many false excludes, and also too many false includes if you flip the test.

---

**Q3.** Leave **affected NULL** on control tables unless the count is **people** and a row/column already names control. Do **not** treat "reference alleles" as N unaffected people.

Table 2 Number=408 for H558R is 408/2600 **alleles** from 1300 volunteers. Gold's 408/0/408 copies that allele integer as people. The contract says population alleles are not carriers. Emitting 408 unaffected humans will match this gold and then wreck BRCA/ESP/ExAC polymorphism tables.

The parser **already** emits `affected=0`, `unaffected=N` via `row_control_like` (`Rare control`/`Polymorphism`). The guard then kills the 0 (`unsourced_zero_affected`) because `affected_count_type` is still the empty affected header, and `_phenotype_is_sourced` only accepts affected types `{case, proband_count}`. Minimal fix: when `row_control_like`, stamp `unaffected_count_type=unaffected_control` **and** teach the guard that `affected==0` is sourced when the paired una is a control count equal to carriers. Do not invent `count_type=case` for a control zero.

If you emit a closed N/0/N on allele tables, you also create a 100% "unaffected" clinical claim that publish still reads raw (contract open blocker). Leave the 0 null; keep `unaffected=N` only if you are willing to defend N as people.

---

**Q4.** Within-arm on/off on frozen extraction JSON + frozen source **is** the valid effect estimate for a deterministic post-processor that only fills nulls. It is better than a paired live extract.

The second paid arm that "only differs by this step" is **wasteful**. Run-to-run LLM variance is why the registry wants a baseline for prompt/model changes. This is not that. One locked extract + on/off ablation is the experiment. Spend the $10 on a **transfer gene** (MYBPC3 already cited in the case-series comment; BRCA1/2 or BMPR2) or a holdout of family/cascade papers, not on re-extracting cardiac gold you can already rescore for $0.

Your +10 pp / 85% pooled conditional gate is sandbagged. Calibration claims +37–42 pp. 85% pooled lets the 45 center-split wrongs plus new leaks pass. Set instead:

- Primary: within-arm Δ exact affected / 1,118 positive-gold rows. Forecast a number from the $0 replay **before** opening; fail if live Δ is < half the replay Δ.
- Secondary (the precision the PI asked for): exact/supplied on the **newly supplied set only** ≥ 90%.
- Hard guardrails: 0 new affected integers on rows where gold `affected != carriers` (the 100-row real-split class, including 28237968); identity TP/FP/FN unchanged; **counted-extra** precision non-decreasing (stamping affected on extras turns them counted); una exact on **positive** gold una non-decreasing; 0 hits on a preregistered family/cascade/allele smoke list.

Reject even if recovery jumps: any 28237968-style N/N/0; allele counts becoming people; a drop in delta-set exactness below 90%; counted extras up; BRCA/MYBPC3 smoke fail; or you had to special-case PMIDs.

---

**Q5.** Next single lever: **stop the identity-only fixed-width fast path (30059973 and kin) so phenotype/count tables are parsed.**

The 477 plan cannot touch rows with no carrier integer. 232 identity misses are a source/matcher program, not one change. The 45 center-splits are a gold granularity bug (Table 4 has `Testing center`; gold stored four `1`s). Chasing them trains the model on the wrong unit. 126 matched-nulls are where counts exist in the paper and the parser returned early; 185 variants on one paper is the densest remaining pocket, deterministic, no LLM, and unblocks carriers **and** affected together.

Do that **instead of** a new `table_cohort_phenotype.py` if you only get one patch: rebind existing case-series headers to the affected role, reconstruct wrapped captions (`No. of`/`patients`, `### Table 4` → nearest real caption, not the "genotype-positive" prose hit), and stamp control zeros. That is the 26496715 / 11854117 / Kapplinger slice without manufacturing a partition.

---

**Q6.** The 55–60% forecast is unproven and likely incompatible with the exclusion list you wrote.

`(203+477)/1118 = 60.8%` assumes every one of the 477 gets a CASE stamp. `genotype-positive` on Kapplinger's Table 4 intro drops 239 rows → ~39%. `frequency` drops 27566755's 47. `allele` before CONTROL drops the 55 control zeros. I did not see a paper-level breakdown of the 477; A–E plus leftovers must be shown from the locked JSON, not from memory.

Also check, on the two locks, with no new extract:

1. Count of the 477 with `carriers_count_type==per_variant_carrier` and stored caption in `{Table N, empty}` vs a real caption.
2. Whether `_is_case_series_people_header` is already true on those carrier headers (if yes, the module is redundant with a role rebind).
3. Table 4 people-column sum vs 438+13 (you asserted 451; verify).
4. How many of 185 fast-path variants at 30059973 gold-match (cannot all sit in the 126).
5. `~87+62` mixes gold-affected=0 rows with gold-affected>0 rows; they are not one error budget.

Minor: "nothing about identity changes" ≠ counted-extra unchanged. Publish still reads raw counts. TASKS.md's next gate is the **patient-row derived-count contract**, not this copy rule.
