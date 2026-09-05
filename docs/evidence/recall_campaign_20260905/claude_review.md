Write is disabled in this session, so here is the review inline.

**Consult status:** incomplete — this session has no Bash tool, so `grok` and `agy` could not be invoked. This is the Claude arm only. §6 says how to finish it.

---

## 1. Changes most likely to raise TP / cut FN, with precision protection

A–C are **$0 to measure on locked predictions**, because they are deterministic post-model stages — the model output is already on disk, so replaying them carries *zero* extraction variance. That property is the most important budget fact in this review.

**A. Stop the VariantFeatures gate deleting correct, gene-attributed variants — and stop it depending on which fields the model happened to emit.**
`scripts/enrich_from_variantfeatures.py:397-402` returns `wrong_gene_residue_mismatch` when a residue matches *any* other known gene at that position; `pipeline/vf_enrichment.py:38` quarantines that class by default. `docs/evidence/mixed_gold_cont120_01_20260903.md:119-126` is decisive: KCNQ1 30462975 `p.G45S c.133G>A` was quarantined in the candidate arm while the baseline's cDNA-free `p.Gly45Ser` was kept as `novel_in_range` and scored a TP. The same gold variant is accepted or destroyed depending on whether the model wrote the cDNA — an FN source *and* a nondeterminism amplifier.

Two separable fixes:
1. *Field-invariance (do first, near-zero precision risk).* Classify from one canonical resolved identity, not from whichever notation fields are populated. If the protein-only form would classify `novel_in_range`, adding a cDNA must not flip it to quarantine unless the cDNA is itself out of range.
2. *Transcript tolerance (measure before shipping).* `docs/evidence/vb_publication_readiness_20260819.md:44-48` already records that `wrong_gene_residue_mismatch` is 504/230/64 rows on BRCA1/BRCA2/BMPR2 while the *confident* cross-gene subset is 171/181/3. Quarantine only when the row is not paper-derived-and-target-attributed, or when the cross-gene match is exclusive; otherwise demote to the report-only `residue_unverified`. This will raise BRCA/BMPR2 live rows and can cost precision there — measure per gene, never pooled, by re-running `enrich_and_quarantine` on **copies** of locked DBs and diffing `<GENE>_vf_false_positives.csv`.

**B. Paper-derived origin must outrank linkage origin at fold time.** One cont120-01 miss was a correct paper extraction reclassified out of the primary lane by linkage-origin attribution (`mixed_gold_cont120_01_20260903.md:128-131`). Shared resolver is `pipeline/variant_identity.py`; lane assignment is in the provenance-lane projector (TASKS.md:248-258). When two rows fold with differing origins, the merged row keeps the paper origin. *Uncertainty, flagged not invented:* I did not open the projector — confirm the actual lane write site before editing rather than assuming it sits in `variant_identity.py`.

**C. Promote source-only single-letter substitutions.** `utils/legacy_notation.py:71-111` already covers `<Ref><pos>fs[X<n>]` and `<pos>(ins|del)<AA>`, firing only when `cdna_notation`/`protein_notation`/`legacy_notation` are all empty (`:92-94`) — that emptiness check is the precision guard and it is already tested. Add a third shape for the KCNQ1 26481773 `R594Q` class: `^(?P<ref>[ACDEFGHIKLMNPQRSTVWY])(?P<pos>\d{1,5})(?P<alt>[ACDEFGHIKLMNPQRSTVWY*X])$` → `protein_notation`. Refuse `ref == alt`, and gate the position against the gene's known protein length rather than accepting any 1–5 digits — that is the whole precision story here. `pipeline/prompts.py:131-134` already forbids bare single amino-acid letters, so this only reaches rows the model deliberately parked in `source_notation`. Small yield, but it removes a known arm-to-arm flip.

**D. Bare table cDNA** (`14803 G>A` → `c.14803G>A`): scope to a cell under a nucleotide/cDNA-change **column header** using the existing header/caption-scope machinery in `pipeline/table_router.py`. Never in free prose — an unscoped numeric-plus-`G>A` in prose is how you manufacture extras.

**E. Prose exclusion lists** (`A1136V-RyR2 [3 cases]`): `utils/variant_scanner.py`. Require the trailing gene token to equal the target gene or a known alias, and require the variant to validate.

**F. Acquisition — still the largest pool, and it is not "reading."** 70/87 tranche-01 and 72/123 cont120-01 misses are acquisition. Two defects are already reproducible (`phenotype_failure_panel_20260905/README.md:81-91`): `harvesting/supplement_reference_parser.py::parse_supplement_references` sees only `Table S14` in the rich 30059973 context and misses the repeated **Online Table** references, and its gap logic compares table counts to file counts — different units; and the two HTML objects under 25163546's supplements are the EHJ-Supplements journal page and OUP advertising, accepted as supplements today. Validate parent DOI/title, container linkage, MIME/magic bytes and payload. Measure this at $0 with `gold_source_presence_sweep.py` and `fn_root_cause.py` — **bytes recovered is not recall**, and folding it into a reading arm makes both unmeasurable.

---

## 2. Design for a fresh $100 at ~120 attempts/tranche

The binding constraint is resolution, not budget. At ~$10–12/arm you can buy ~8 arms, and the campaign's own conclusion (summary §8, final paragraph) is that a single paired arm at n=49–120 cannot separate a reading change from temperature-0 provider variance (±2 TP, ±5–7 count values). Three more paired arms of the same shape produce three more `reject_or_revise_candidate` rows and no knowledge.

**Phase 0 — $0, before anything is opened.**

| Step | Measures | Tooling |
|---|---|---|
| 0a | A, B, C, scorer bridges, on **locked predictions** | replay deterministic stages on DB copies; re-match with `cli/compare_variants.py` |
| 0b | per-gene precision cost of A(2) | `enrich_and_quarantine` on DB copies, diff `<GENE>_vf_false_positives.csv` |
| 0c | acquisition delta (repository fallback + F) | `gold_source_presence_sweep.py` before/after on the cont120-01 FN worklist |
| 0d | count-side guard precision | `scripts/phenotype_value_precision.py --compare`, quoting `kills / destroys` |

**Gate:** if A+B+C cannot show a positive locked-prediction delta at $0, they will not show one at n=120 against ±2 TP of noise, and no paid arm is justified. D/E/F cannot be fully measured this way and carry into Phase 1.

**Phase 1 — paid, replicated, ≈$66.** On one unopened 120-attempt tranche: 2 × baseline replicates + 2 × candidate replicates ≈ $44. This is exactly the change the campaign asked for — replicates *estimate* within-arm variance instead of assuming it away, so the paired difference is tested against a measured SD. Then 1 baseline + 1 candidate on a second unopened tranche ≈ $22 as confirmation. Leaves ≈$34 reserve.

*Cheaper alternative worth checking first:* if the active provider exposes a sampling seed, fixing it collapses the replicate requirement and halves the cost. Check `utils/llm_utils.py` and `config/settings.py`. **I did not verify this — open question, not a fact.**

**Separating acquisition from reading is now mandatory**, because the repository fallback is a live network path *inside* the extraction run. If the candidate arm gets new bytes and the baseline does not, you are measuring acquisition and calling it reading. Run both arms against a pinned identical corpus snapshot (`--no-source-recovery` with `--pmid-file`, per CLAUDE.md's calibrated-measurement guidance), record per-attempt source hashes, and report acquisition separately as the 0c number.

**Preregister before opening.** Primary: paper-derived identity recall, observed delta ≥ +1.0 pp *and* lower bound > 0 against the replicate-derived SD, paired by paper. Guards: precision LB ≥ −2 pp, observed recall delta ≥ 0. Secondary: e2e carrier MAE per `mixed_gold/secondary_endpoints.json`. Never loosen the identity rule; never re-adjudicate a burned tranche. Do **not** ship a second cumulative multi-fix bundle — the campaign already showed a 6–8-fix bundle at this n is unresolvable.

---

## 3. Pending repository fallback: gaps worth fixing before freezing

The design is sound (two independent identity witnesses, PDF magic bytes, bounded budget, retained PDF + hash). These are specific defects.

**P1 — `supplement_surface_status` is unconditionally `unavailable`, even when supplements are present.** `harvesting/repository_pdf_fallback.py:474` writes the marker and `:502` sets the artifact key regardless of the `supplement_markdown` argument. Both callers can pass real supplement text: `harvesting/orchestrator.py:1823-1825` and `scripts/fetch_paywalled.py:1303-1312`. Consequence: `pipeline/source_quality.py:210 is_reusable_fulltext_source` returns False permanently, so corpus reuse refuses it and **every later run re-fetches a paper it already has** — cost, latency, and a source-nondeterminism channel between arms. It also mislabels the artifact for the component-status work TASKS.md:152-162 makes next. The README claims "preservation of existing supplements" in its acceptance list — check whether that test asserts the *status* or only that the text is present.

**P2 — `main_text` artifact record is replaced, not merged.** `write_repository_source:486-501` does `artifact.update({... "main_text": {...}})`. On the content-validation-failure path the artifact was already `.save()`d at `orchestrator.py:1822` carrying `figure_captions`/`table_captions`/`supplement_descriptions` from `record_main_text` (`:1803-1809`); those are nested inside `main_text` and are lost. Relatedly the repository text discards the captions block `_build_unified_content` produced (`:1797-1801`) — captions are deliberately placed contiguous with the body for the scanner and LLM. Merge into `main_text`, and carry captions through the way `supplement_markdown` is carried.

**P3 — `KeyError` on a redirect without `Location` escapes to the caller.** `_fetch:352` does `response.headers["Location"]`. `KeyError` is in neither `_candidate`'s except tuple (`:456`) nor `recover`'s (`:147`), and `orchestrator.py:1899-1901` does not wrap the call — a malformed 3xx crashes the harvest for that paper. Use `.get("Location")`. Also add an explicit hop limit: the redirect recursion terminates only via `MAX_ASSET_REQUESTS`, so a loop burns the paper's whole 14-request budget before any real candidate is tried.

**P4 — a fifth boundary-free DOI regex, against a landed decision.** `article_identity:101` uses `re.findall(r"\b10\.\d{4,9}/[^\s<>\"]+")` with an ad-hoc `.rstrip(");,")`. The project consolidated exactly this into `utils/doi.py` because four modules each carried a boundary-free copy and a glued-on next word produced an unresolvable identifier that failed a real validation run (TASKS.md:762-767). Reintroducing it in a *rejection* path is worse: a header DOI rendered as `...853374Received` won't equal `doi`, so the function returns **"article header has a conflicting DOI"** for the correct paper. Use `utils/doi.py`.

**P5 — ordering: repository now preempts the browser/EZproxy path.** `_download_with_tier35_fallback` (`orchestrator.py:1866-1895`) returns immediately on a body-only repository success (`:1874-1876`), so the browser path — the one that can return the publisher version *with supplements* — is never tried. Given that the missing gold is concentrated in table-body capture on identifiable papers, silently downgrading a would-be publisher fetch to a body-only PDF is a plausible recall regression. Try repository *after* browser, or continue to browser when the result is body-only and keep the richer one. **Design question, not a proven regression** — the 14-paper panel used no browser comparator, so nothing here measures it.

**P6 — cover-page detection is an enumerated phrase list.** `article_identity:83-91` hardcodes four repository boilerplates (HAL, Wiley, DSpace/Pure, DTU/Pure). Same open-vocabulary shape as `_gene_symbol_tokens` (TASKS.md:843-852), where extending the ignore-list demonstrably did not converge. It is load-bearing because covers *do* print the title. Bounded improvement: widen `pages[:3]` to ~5 (a 4-page cover currently produces a silent false rejection) and require the accepted page to carry an author/affiliation/abstract cue instead of enumerating more cover phrases.

**P7 — flat page text only, and no supplement-gap check.** Conversion is `page.get_text(sort=True)` (`:419`); nothing routes the retained PDF through the table-aware converter. This matters because `harvesting/format_converters.py::pdf_to_markdown` already returns on the first successful text converter *before* its table-aware fallback (phenotype panel README), so nonempty text is not evidence of faithful tables. Separately, `_download_repository_pmid` returns before `check_supplement_gap` (`orchestrator.py:1840-1849` runs only on the PMC path), so a recovered body advertising "Online Table 1" raises no gap warning — precisely the signal the next task needs. The README is honest that table fidelity is unproven; the code should record that status.

**P8 — per-run state on the instance.** `self.deadline`/`self.asset_requests`/`self.visited` are set in `recover()` (`:122-124`), not `__init__`. Single-use, not thread-safe. Current callers are safe (`evaluate_repository_fallback.py:44` builds one per paper even at `--workers 3`), but `orchestrator.py:1899-1901` hands it a shared session. Move the per-run state into a local context, or document the contract.

**P9 — minor.** `_public_url:73` accepts any hostname containing a dot without resolving it (per-hop re-checking limits exposure; write the residual down). `_metadata:202-203` *aborts* recovery on a supplied-DOI/Europe-PMC conflict; preferring the Europe PMC value and recording the conflict would recover papers whose supplied DOI is simply wrong.

---

## 4. Leakage, nondeterminism, budget pitfalls

**Leakage.** `gvf_data/kcnh2_variant_aliases.json` is gold-derived and loaded at runtime by `_lookup_alias` (`utils/variant_normalizer.py:1588`) for every gene — every new arm must run `--gold-free-run` (TASKS.md:231-247). More urgently: the 22-paper phenotype panel, the 16-paper repository panel, and the cont120-01 problem-paper worklist were all selected from *opened* locks and read directly. They are permanently calibration evidence. **Before opening any confirmation tranche, diff its PMID membership against those three lists** — the cont120 tranches were built to exclude the 90 consumed PMIDs, not these. Opening or scoring burns a tranche, so finish all of Phase 0 and freeze the candidate first.

**Nondeterminism.** ±2 TP and ±5–7 count values per tranche at temperature 0 is the same size as every effect sought so far — that is the entire argument for replicates. Fixing the VF gate's field-dependence (A.1) *reduces* measurement noise as a side effect, which is an independent reason to do it before the paid arms. New this cycle: the repository fallback does live network I/O inside the extraction run, so arms run days apart can legitimately read different bytes — pin the corpus snapshot and record per-attempt source hashes. P1 above makes recovered sources permanently non-reusable, turning every run into a re-fetch: a nondeterminism source disguised as a caching bug.

**Budget/process.** ~8 arms available; §2 spends ≈$66 on 6 and holds ≈$34 — don't spend the reserve on a second protocol variant. Use the same public-price proxy as the $84.035 ledger or the numbers aren't comparable. **Never commit while a baseline is frozen** — a docs commit during cont120-01 swept nine files back to `506a949c` and forced a candidate re-lock; check `git status` and stage explicitly rather than a bare `git commit`. Corpus preflight `test -L corpus && test -d corpus` before any corpus-touching job.

---

## 5. Recommended order

1. **$0** — Phase 0a–0d on locked predictions and DB copies. Gate as above.
2. **$0** — fix P1–P4 and P7's status recording; re-run `tests/unit`. Decide P5's ordering explicitly and record the reason.
3. **$0** — Phase 0c acquisition measurement on the cont120-01 FN worklist, reported separately from recall.
4. Freeze the candidate; verify tranche membership against the three contamination lists; preregister the §2 rule.
5. **≈$44** — 2+2 replicated arms, one unopened tranche, frozen matched source.
6. **≈$22** — confirmation pair, second unopened tranche, identical candidate runtime.
7. Append `docs/PROTOCOL_CHANGELOG.md` (one row, linking the gold-difference figure), update `RECALL_STATUS.md`, append `RECALL_HISTORY.md`.

## 6. Completing the consult

Run grok and agy as **bounded, facts-only** reviews of a self-contained brief, not repo-reading tasks, and reproduce any paper-level claim before accepting it — the 2026-09-05 panel had to reject four of agy's unverified paper-level diagnoses outright. Ask both specifically about: (a) whether A(2) transcript tolerance is safe without a BRCA re-measurement; (b) 2+2 replicates vs. a provider seed; (c) the P5 ordering question.

**Verification for any implementation:** `.venv/bin/python -m pytest tests/unit -q` after steps 1–2; locked-prediction replays write to **copies** only; any scored arm regenerates the stratified figure and the companion `figures/gold_difference.*` per `docs/PHENOTYPE_COUNT_FIGURE_POLICY.md`.
