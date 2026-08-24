# GVF Handoff Tasks

Last reviewed: 2026-08-22.

This is the only active GVF checklist. Current measurements and caveats live in
[`docs/RECALL_STATUS.md`](docs/RECALL_STATUS.md); completed benchmark history
lives in [`docs/RECALL_HISTORY.md`](docs/RECALL_HISTORY.md); protocol changes
are append-only in
[`docs/PROTOCOL_CHANGELOG.md`](docs/PROTOCOL_CHANGELOG.md). How to raise
gold-120 **identity** precision and cut cost (which denominator, the 10 counted
extras, Sol/Grok/Kimi shares, $0 vs paid) lives in
[`docs/PRECISION_COST_LEVERS.md`](docs/PRECISION_COST_LEVERS.md). How to raise
**affected/unaffected value** exact-match precision (not Gate 2; no one-off
aliases) lives in
[`docs/AFFECTED_UNAFFECTED_PRECISION.md`](docs/AFFECTED_UNAFFECTED_PRECISION.md).
Dated reports and benchmark run directories are evidence, not competing plans.

Execute the acceptance gates in order. Do not publish a new headline, update the
dashboard, enable a default-off recovery stage, or advance to a larger cohort
until the preceding gate passes and its artifact is recorded.

## 1. Re-establish the scientific baseline

The paired live cardiac rescore and Gate 1 review were recorded on 2026-08-13.
The lead approved advancement from Gate 1 without changing the published
headline. The requested 100--150-paper scale applies to comparison against the
already manually curated cardiac gold standard. Gate 2 is therefore the fixed,
gold-value-blinded `gold_120` sample: originally 30 source-available,
count-eligible papers per cardiac gene (seed 2026081301). Live membership is
118 attempts / 114 unique PMIDs after the non-genetics KCNH2 PMIDs 10086972
and 14642689 were quarantined. It does not widen the experimental genes.

The last completed Gate 2 **passed** under the 2026-08-20 no-inference
contract. The extraction-blinded lock is
`benchmarks/codex_paper_eval/runs/20260820_gold119_noinference`: counted-extra
precision is 97.51%, above the 77.3% floor; the distinct count-bearing-only
diagnostic is 93.30%. Variant recall is 86.57% (548/633), up from 84.09% on the
accepted 2026-08-13 run. Carrier coverage is 30.49%, but conditional carrier
MAE regressed to 0.425 from 0.308. Affected and unaffected coverage deliberately
fell to 7.58% and 5.06% because the pipeline now refuses to infer partitions
from diagnoses, enrollment totals, or arithmetic. This is an identity
recall/precision pass, not a new penetrance or phenotype-count baseline. The
run retained all 42 PMID 26746457 classification identities with null counts
and recorded circuit-breaker empties rather than silently dropping them.

The clean candidate-local replacement is locked at
`benchmarks/codex_paper_eval/runs/20260821_candidate_local_gold119`. It was run
gold-free over the then-live 119-attempt manifest and scored only after a
hash lock: 546 TP / 290 FP / 87 FN (86.26% recall, 65.31% raw precision),
98.03% counted-extra precision, 94.44% count-bearing-only precision, and
carrier MAE 0.255. All four production traces are write-time verified and
bound into the lock (578 LLM calls / 2.552M tokens). The blind circuit breaker
correctly returned no KCNH2 result for unrelated PMID 14642689; PubMed
adjudication then quarantined that record from the live tier without rewriting
the lock or underlying gold snapshot. New Gate 2 runs therefore use 118
attempts, while this immutable run remains the accepted evidence for the code
change.

- [ ] **Act on the 2026-08-24 blind proof result, then stage the reviewer
      surface.** The proof run is DONE and locked at
      `benchmarks/codex_paper_eval/runs/20260824_postfix_gold118` (118
      attempts, clean tree `e4910d9`, 42 min, $11.32 list-price proxy).
      Identity is flat versus the 2026-08-21 lock (TP 546 = 546, FP 290→284,
      recall 86.26%→86.39%); count coverage rose (count-bearing rows
      187→208) and matched-row MAE rose with it. Full comparison and the
      source-mix control are in `docs/RECALL_HISTORY.md` 2026-08-24.
      Decisions this leaves open:
      (a) adjudicate the two single affected values off by 7 (KCNH2
      22338672, KCNH2 16029385) — they carry 14 of the +17 error units,
      while RYR2 30403697 and KCNQ1 26496715 each gained 5 exact values;
      (b) decide whether the coverage/accuracy tradeoff is accepted as
      shipped, or whether extraction should emit fewer marginal counts;
      (c) claim verification grew 188→232 calls because more rows carry
      counts, absorbing the vision lane's −22% saving — batching a paper's
      cards into one call is the standing lever if cost needs to fall.
      Only then stage the reviewer surface.
- [ ] **(superseded, kept for the staging half) Stage the 50/50/50 review
      set.** The batch itself (cruft sweep, prompt-caching split +
      verifier-effort fix, figure-vision cost controls, acquisition
      fail-open closures, eight generic notation/scanner FN fixes) is
      recorded in `docs/PROTOCOL_CHANGELOG.md` 2026-08-22 rows with
      evidence in `docs/evidence/strategy_consult_20260822.md`. Remaining,
      in order: (1) Sol + Grok adversarial review of the full diff, fix
      findings, commit; (2) paid blind proof run
      (`setup_production_eval.py create --run-id 20260822_postfix_gold118`,
      118 attempts — the manifest already excludes 14642689 — then
      `run_extraction.sh` and `lock_and_score.sh`), compared to the
      2026-08-21 candidate-local lock (546/633, 98.03% counted-extra,
      carrier MAE 0.255, ~$9.77) on recall, counted-extra precision, MAE,
      and $ (vision usage is now captured); budget ceiling $100, expected
      ~$6-10; (3) only if the cardiac metrics hold or improve, stage the
      50/50/50 BRCA1/BRCA2/BMPR2 candidates into the login-gated Variant
      Browser reviewer surface (NOT public publish), re-running the
      contamination audit (canine PMID, wrong-gene rows, paper membership,
      fixed-code residue QC; the 2026-08-22 pre-audit found only one
      curator-visible BRCA1 boundary flag, p.Asp1864Tyrfs at the stop
      position) against the exact DB imported, and deciding BRCA2 final5
      vs the staged 3-PMID refresh rebuild explicitly.
- [ ] Finish Gate 3, `reviewer_545`: 545 attempts across 506
      unique PMIDs in the 11 populated reviewer workspaces. The first approved
      experimental strata remain BMPR2 50, BRCA1 50, and corrected BRCA2 50;
      do not expand them to 100 papers per gene. **The complete blind extraction
      and replay process is now staged for exactly 50/50/50 papers.** Pinned
      manifest membership, staged JSON membership, and final DB paper membership
      agree exactly for every gene. PMID 19944633 (the canine BRCA2 paper) is
      absent from all three surfaces. The final local candidate DBs are BRCA1
      `BRCA1.final6_20260821.db`, BRCA2 `BRCA2.final5_20260821.db`, and BMPR2
      `BMPR2.final5_20260821.db`; they contain 3,582 / 722 / 260 live variants
      after mandatory VariantFeatures enrichment and high-confidence
      out-of-range/wrong-gene quarantine. Structural audit: zero placeholder
      titles, wrong-gene live rows, nameless identities, negative counts,
      impossible unquarantined partitions, or duplicate penetrance strata.
      Typed family-count evidence remains raw: 168 / 154 / 3 links carry family
      provenance, and 147 / 142 / 3 of those links carry a raw value. The trust
      gate masks those values from carrier-facing totals where the role is
      `family_count`.

      The refresh path now has an independent, default-on source-evidence gate:
      selection and matching use source identity/provenance only, never gold.
      It fills only a unique identity with a null or identical untyped value,
      restores the backup on absent/ambiguous/conflicting untyped evidence, never
      aliases BIC digits to cDNA, and never matches by codon position. A shared
      ref+position+alt rule folds a truncated `V340G` representation into
      `p.Val340Glyfs*6` across regex/table merge, replay, and same-cDNA SQLite
      migration; protein-only-to-cDNA matching stays fail-closed because
      streaming uniqueness is order-dependent. Grok 4.6 `xhigh` found and then
      cleared the prior P0 replay/migration blockers; its final live-path finding
      (regex had retained the rich frameshift only in `source_notation`) is now a
      regression test. Details and remaining explicit conflicting-protein groups:
      `docs/evidence/brca_bmpr2_full150_audit_20260821.md`.

      **Still HOLD for public publication.** Source-ground the collaborator
      precision sample, re-review the 111 detached BRCA2 subjects, and run a
      disjoint pre-registered cardiac/SCN5A validation before public promotion.
      The trusted Variant Browser import/publish path is fail-closed, but no
      public annotations were changed. The requested Claude desktop adversarial
      pass remains blocked because macOS is locked; Grok CLI review is complete.
- [ ] After an accepted rescore, update `docs/RECALL_STATUS.md`, append
      `docs/RECALL_HISTORY.md` and `docs/PROTOCOL_CHANGELOG.md`, and regenerate
      the dashboard. Until then, the public dashboard remains an archived
      pre-correction snapshot.
- [ ] **Gold-120 refinement follow-through.** Ranking and “do not”
      constraints: [`docs/PRECISION_COST_LEVERS.md`](docs/PRECISION_COST_LEVERS.md).
      Paper-level evidence:
      `benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/`
      (`PRECISION_AND_COST.md`, `NOTES.md`, `remaining_fn.tsv`). Diagnostic
      rescore (locked predictions, live gold + current matcher, not a new
      lock): recall 545/633 (86.10%), raw precision 40.85%, **counted-extra
      98.55% (the Gate 2 / intended precision number)**, carrier MAE 0.292.
      Was 98.20% before the 2026-08-16 deletion-span bridge and the
      document-level gene-attribution pass closed two of the ten extras.
      Remaining, in this order — **do not start the paid 119-paper re-extract
      to chase counted-extra precision**:
      1. Curate the **8 remaining counted extras** (only Gate 2 FP surface
         left, and now all curator work rather than code). List and classes in
         `PRECISION_AND_COST.md` §2a and
         `docs/GOLD_CURATION_QUEUE_2026-08-14.md` §8: two synonymous counted
         rows, four identity TPs against gold `0/0/0`, and the
         `GOLD_GAP_REAL_VARIANT` identities. Do **not** drop synonymous rows
         as a class — gold carries them with real counts (KCNQ1 `A344A`, 23
         carriers on 30758498). The other-gene `L187P` and the `c.693delCA`
         near-twin are closed in code (2026-08-16).
      2. Remaining gold-queue items (`docs/GOLD_CURATION_QUEUE_2026-08-14.md`):
         locate the intended source for quarantined 14642689 before any
         restoration, adjudicate 11 editorial rows, and resolve compound
         `W248F + L347R`. Do **not**
         collapse the 19216760/24394973 exon-3 family pairs (two families).
         Do **not** delete the 780 identity-only extras to raise raw 40.82% —
         584 sit on papers that already matched every gold row. If raw
         precision is the target, expand gold (26746457 first; 24
         source-confirmed).
      3. Scorer-only generic bridges. `c.693delCA` ↔ `c.692_693delCA` is
         **done** in both scorers and in twin identity (endpoint match on
         identical deleted bases; ambiguity refused; single-base deletions
         never bridge). RYR2 `c.169-198_273+820del` ↔ EXON 3 is still open.
         Port the remaining 2026-08-14 matcher rules (arrows, stop/fs,
         splice/translation, twin-merge) from `run_eval.py` to
         `cli/compare_variants.py`. Structural alleles already match via
         `utils/structural_alleles.py`; the rescore added 0 structural TPs
         (those identities were not extracted). No per-variant aliases.
      4. **DONE 2026-08-22:** $0 Sol audit of the 244 locked figure-vision
         calls: 229 (93.9%) added no identity the text extract already had
         (replicated 92.2% on the 08-21 lock); 47 were byte-identical
         re-reads. Redundant-call cost ≈ $2.34/pass. Countermeasures shipped
         (SHA dedup + fail-open caption triage + cap-retry); see
         `docs/evidence/strategy_consult_20260822.md` §A1.
      5. **DONE 2026-08-21:** paid blind re-extract over the then-live
         119-attempt tier. The candidate-local lock records 546/633 recall and
         98.03% counted-extra precision. The unrelated 14642689 paper was
         quarantined afterward, so future live runs contain 118 attempts.
      6. Targeted acquisition: caption-stub table bodies (RYR2
         19926015, KCNQ1 14678125, 31520628, 24667783), abstract-only
         stratum from the 2026-08-14 analysis.
         **CORRECTION 2026-08-17 — KCNH2 29650123 is NOT a figure-OCR item.**
         This entry, and the "20 FNs / source not readable as text (TIFF
         tables) / needs figure OCR" row in
         `runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/NOTES.md`,
         were inferred from "`mmc1.docx` converts to text with zero variant
         tokens" without opening the images. The three embedded TIFFs were
         rendered and inspected: they are a stacked bar chart of patients per
         decade of presentation and two hazard-ratio-vs-baseline-QTc curves
         (restricted cubic spline + linear). There is **no mutation table in
         any image**. `document.xml` holds only "Online Table 1.
         Characteristics of the study population". Gold expects 22 rows for
         this PMID (all `carriers=1`); the acquired `_FULL_CONTEXT.md` contains
         exactly one of them (G628S). The mutation list was **never acquired** —
         this is an acquisition gap, not an OCR gap, and figure vision on those
         TIFFs would return zero variants. Reclassify the 20 FNs and re-derive
         the "figure OCR → ~91%" rung of the gold-120 ceiling ladder, which
         inherited the same error.
      7. Deterministic `regex_table` count-column binder (26496715
         wrapped-header; 114 gold assertions). That paper has 0 counted
         extras — this is MAE / count-recall, not Gate 2. Short-circuit
         refuse already landed; binder is not done. Only after the binder:
         Luna-max post-layer count-ambiguity cards.
         **PARTLY LANDED 2026-08-17 — the diagnosis was wrong and the real
         cause was worse.** The blocker was never a wrapped header alone.
         `26496715_FULL_CONTEXT.md` (a `.doc` supplement converted to text)
         contains **zero** `|---|` separator rows, and
         `enumerate_markdown_tables` required one, so all three mutation tables
         were invisible — the router literally never saw them and the paper
         scored "no usable variant tables". Two deterministic fixes landed in
         `pipeline/table_router.py`:
         (a) a separator-less ("borderless") pipe-block branch that also
         rejoins a header wrapped across physical lines (`No. of` + `patient` +
         `s`), claiming only runs that contain no separator so bordered tables
         are still owned by the existing path and never emitted twice; and
         (b) a caption-scope exemption for the row-level off-target-gene guard.
         (b) was necessary because `_gene_symbol_tokens` accepts any 3–12
         character uppercase token, so `Missense`, `N-Terminal`, `CNBD`, `DII`,
         `DI-S6`, `PAC`, `FRAME`, `SHIFT`, `SITE` and even the raw nucleotide
         run `CGGGGCGAC` all read as "some other gene" and
         `_row_has_off_target_gene_without_target` deleted **every row**. Any
         variant table with a mutation-type or protein-domain column but no
         mapped Gene column lost all of its rows. That vocabulary is open-ended
         and an ignore-list does not converge on it (extending the list by name
         still left 33/99 rows dead), so the guard is now skipped when the
         caption already scopes the table to the target gene — a stronger
         signal, and the misrouted-panel case it was written for (unscoped
         caption) still runs it.
         **$0 measurement, no model calls, source already on disk:** tables
         visible to the router on 26496715 went **0 → 5**, and gold rows
         recovered with an identity *and* an exactly-correct carrier count went
         **0/99 → 87/99 (87.9%)** — KCNH2 48/53, KCNQ1 39/42, SCN5A 0/4, and
         100% of supplied counts exact. Scored with
         `utils/variant_normalizer.normalize_to_single_letter`, not
         `cli/compare_variants.py`, so some of the 12 misses may be matcher
         rather than extraction gaps. 12 new pinning tests; offline suite
         1978 passed / 1 skipped. Still open: the SCN5A arm (0/4), an
         end-to-end re-extract to confirm the short-circuit refuse stops firing
         and the counts actually ship, and a corpus sweep — **132 of 22,148
         `_FULL_CONTEXT.md` files are in this borderless class** (6,846 have
         separators), several with 77 pipe rows.
      8. Residual extraction misses: RYR2 28798025 G1885E (in source, still
         missed), KCNQ1 31293497 A590T (absent from acquirable source —
         gold provenance).

- [ ] **Affected/unaffected value precision (active goal).** Metric and
      rejected wider-NULL rules:
      [`docs/AFFECTED_UNAFFECTED_PRECISION.md`](docs/AFFECTED_UNAFFECTED_PRECISION.md).
      Exact-match among supplied values on the locked gold-120 predictions,
      not Gate 2. Baseline 51/74 affected (68.9%) and 37/56 una (66.1%).
      The always-on guard (`pipeline/phenotype_count_guard.py`) raises that to
      **40/51 affected (78.4%)** and 27/40 una (67.5%) on the same lock, and
      applies on the next extract. Four rules: family copy, figure copy, a
      non-closing partition (`affected + unaffected != carriers`), and an
      unsourced `affected = 0`. The last two are field-scoped to `affected`
      and destroyed **0** exact rows. Measure any new candidate for free with
      `scripts/phenotype_value_precision.py --compare` and quote
      `kills / destroys`, never a bare percentage — every obvious predicate
      (`pred == 1`, `unaffected == 0`, `layer == figure`) is majority-correct.
      Remaining work, in this order — **no per-variant / per-paper aliases**:
      1. Pedigree / figure: count filled vs empty symbols separately; if
         the image is not a counted split, emit carriers only. Largest
         leftover MAE is SCN5A 15671429 (figure 23/2 vs gold 7/15).
      2. Move compact count-semantics verification after merge (already
         §3 below). Target: symptoms ⊂ phenotype (KCNH2 25819988 6
         “cardiac symptoms” vs gold 9/4) and off-by-one splits. **Prerequisite:**
         settle `pipeline/claim_verifier.py` first. `_apply_count_identity_guard`
         implements the forbidden `unaffected = carriers - affected` and is
         correctly skipped when a card is supplied, but
         `_apply_consistency_guards` still runs on every card and can write
         `affected = total` with `unaffected = 0`. Verifying after the merge
         without settling that can manufacture the exact pattern the phenotype
         guard then wipes — on 25819988 it would land `13/13/0`, also wrong.
         The correct action there is to keep `carriers=13` and null the split.
         Free first step: count how many locked cards that rewrite touched.
      3. Do **not** refuse every n=1 `affected=1` from a short evidence
         quote — that dropped 49 exact case reports and lowered
         affected precision. Prompt already forbids defaulting undescribed
         status; enforce with a real quote, not the prediction stub.
      4. Paid re-extract only after 1–2, as a measurement of those
         general rules. It does not by itself fix 15671429.

The three manifests in
[`benchmarks/evaluation_tiers/`](benchmarks/evaluation_tiers/README.md) are the
only rollout cohorts. The curated extraction fixture and dated benchmark runs
remain specialized or historical evidence, not extra gates.

## 2. Recover missing source before adding inference

- [ ] Fetch the 44 never-acquired gold-120 FN rows concentrated in 7 papers
      (routes per paper in `docs/evidence/strategy_consult_20260822.md` §A2:
      JACC/Elsevier online tables for 29650123, Elsevier AJOG supplements for
      31520628 after the mis-bound CDC PDFs were gated, Wiley TDM for
      14678125, AHA/EZproxy for 27114410 and 15466642, Karger/EZproxy for
      17971661, Table 2 body for 19926015), then targeted
      `scripts/refresh_run_db.py` on those PMIDs — not a full re-extract.
- [ ] Complete the one-time EZproxy profile bootstrap with
      `.venv/bin/python scripts/ezproxy_relogin.py --bootstrap`, then rerun
      paywall recovery and `refresh_run_db.py` for the BMPR2 stub backlog.
- [ ] Validate automatic headless EZproxy refresh against the next real session
      expiry; unit tests currently cover only the simulated browser path.
- [ ] Rebuild the two source cards that still lack decisive evidence: RYR2 PMID
      28237968 `c.13352del` and KCNH2 PMID 10862094 `c.526C>T`. Review the other
      seven untested variants in RYR2 PMID 33606749.
- [ ] Work the current source-audit ranking, beginning with blocked SCN5A
      supplements such as PMID 29325976. Use manual PDF acquisition only after
      the current audit confirms the artifact is still missing, then ingest the
      real PDF through `harvesting/format_converters.py`.
- [ ] Re-run the Karger and Sage/Liebert access probes only when credentials or
      publisher access change; the last recorded attempts remained blocked.

Source acquisition is still the dominant recall opportunity. Do not replace a
missing paper or table with speculative count arithmetic.

## 3. Measure count semantics and recovery

- [ ] Join `claim_verification_field_changes` on the locked runs to
      gold-matched exactness ($0): the verify layer changed 80.4% of its
      decisions' emitted values, but changed-vs-CORRECT is unmeasured.
      Required before any verifier threshold, batching, or model change.
- [ ] Move compact count-semantics verification after all extraction and
      recovery layers so it evaluates the final winning observation rather than
      only the early `ExpertExtractor` result.
- [ ] Route only count-bearing multi-cohort or large-count ambiguities. Keep
      broad missing-slot recovery off; its 162-gap probe grounded no additions.
- [ ] Recreate clean DB copies from
      `benchmarks/codex_paper_eval/runs/20260726_fixed48_production`, run the
      hardened carrier-only count-recovery path through the current
      `run_eval.py` validation/lock workflow, and inspect rejections by gene.
- [ ] Treat RYR2 as its own acceptance stratum. Resolve curator questions for
      PMID 29622001 `c.453delC` (24 versus gold 0) and `p.L552S` (73 versus 74)
      before classifying them as recovery defects.
- [ ] Enable `COUNT_RECOVERY_ENABLED` only if the clean measurement improves
      count recall without unacceptable MAE, attribution, or RYR2 regression.
      Keep `COUNT_RECOVERY_FIELDS=carriers`; affected/unaffected widening needs
      a separate benchmark.
- [ ] Resolve the two active BRCA2 gold-scope limitations: family-versus-carrier
      semantics for PMID 26833046 and whether PMID 26848529 is exhaustive
      paper-level gold.
- [ ] Remove the known PMID 33013630 population-annotation contamination from
      canonical data through a reversible, backed-up refresh; its gnomAD and
      predictor columns are not patient carriers.

## 4. Finish the trust and fleet boundary

- [ ] Move report, publish, and other downstream count consumers to the trusted
      projection while keeping raw observations visible for audit. Report
      trusted-tier precision and quarantine rates.
- [ ] Promote count role/evidence type to a first-class persisted axis and have
      the trust gate distinguish patient, cohort-total, control, and population
      observations.
- [ ] Calibrate the trust gate separately on cardiac, BRCA, and one cold-gene
      stratum with `scripts/precision_sample.py`; do not use a pooled result to
      hide a weak stratum.
- [ ] Fail closed for unknown-gene position validation, expand BRCA notation by
      class, and build a BRCA silver standard before claiming generalization.
- [ ] Design a corroboration path that promotes `figure_count_unverified`
      fills in the trusted projection: the masked fills measured 90/91 exact
      (`pipeline/count_repair.py`), far above every destroyed/killed
      break-even. Measure before changing the mask.
- [ ] Variant Browser publish path: one blocked label currently aborts the
      whole publish run with no unblock/skip-and-continue path, and a bare
      legacy label can still merge into a genuine HGVS identity via the
      historical alias keys when both exist. Fix both fail-closed.
- [ ] Fold destructive legacy guards into the trust record so raw values remain
      recoverable, then define and schedule a self-hosted nightly regression
      gate with per-stratum acceptance metrics.
- [ ] Create or import a source-reconciled KCNE1 per-PMID gold input before
      making KCNE1 recall claims.

## 4b. Measurement-artifact integrity (opened 2026-08-17, all $0)

These are defects in the *map*, not the pipeline. They matter because the whole
lever ranking is derived from them.

- [ ] **Regenerate the disagreement artifacts.** `recall_metrics/current/` is
      internally inconsistent: `summary.json` / `report.md` are 2026-08-08 and
      reproduce the published headline (5546/6833 rows, 2596/3010 unique, MAE
      0.6133), but `paper_disagreement_report.csv` in the same "current"
      directory is **2026-07-21** and sums to 2918 matched / 3915 missing rows —
      42.7% row recall, not 81.2%. Root cause is recorded in the summary itself:
      `"disagreement_artifacts_skipped": true`. Consequence: the **failure split
      published in `docs/RECALL_STATUS.md` cannot be regenerated from the
      artifact it cites.** Its classes (568 / 250 / 248 / 184 / 82 / 36 / 8 =
      1,376) do not appear in the on-disk CSV at all, whose classes are
      dominated by `source_missing_or_stub` at 1,251 papers with
      `context_bytes=0` — an unbound run. The 72.8%-acquisition framing that
      anchors the strategy rests on this table, so re-derive it before ranking
      any more work off it. Note also 1,376 ≠ the 1,287 row-recall gap.
- [ ] **Reconcile the per-layer precision table with its artifact.**
      `docs/RECALL_STATUS.md` publishes `regex_table` at **929** counted extras
      / **57.5%** counted precision; `recall_metrics/current/summary.json`
      (2026-08-08) says **899** / **61.60%**, with `matched_db` 1442 vs the
      doc's 1256. Aggregate drifts too: summary has
      `counted_extra_on_gold_pmids` 1637 and
      `precision_vs_counted_gold_pmids` 0.7721, versus the doc's 1,629 and
      77.3% (and a separate prose claim of 1,631). The per-layer column also
      sums to 1,660, not 1,629/1,631/1,637. `regex_table` is the single
      diagnostic that decides where count-attribution work starts, so a 4-point
      drift on it is load-bearing. Grok's 2026-08-17 review additionally flags
      that all of these predate the vertical-table fix and so may be fossils.
- [ ] **Close the open-vocabulary hole in `_gene_symbol_tokens`.** The
      2026-08-17 caption-scope exemption fixes scoped tables only. A table whose
      caption names **no** gene still runs the row-level guard, and that guard
      still reads `CNBD`, `DII`, `PAC`, `FRAME`, `SHIFT`, `SITE` and raw
      nucleotide runs as off-target gene symbols. Note `_caption_gene_scope`
      resolves against the 13-gene built-in roster only, so an off-roster
      symbol (LMNA, TTN, BRCA2, ENG) produces an **empty** scope and its table
      reads as unscoped. The durable fix is column-level: treat a token as
      gene evidence only when its whole column looks like gene names. Do not
      extend the ignore-list further — it demonstrably does not converge.
- [ ] **Remove the three agent worktrees.** `git worktree list` shows
      `.claude/worktrees/code-review-max-e46bda`,
      `gold-standard-precision-ad637f`, and
      `review-works-improvements-inefficiencies-1cad7f` on `claude/*` branches.
      The handoff policy is one checkout on `main` with no worktrees; branch
      sprawl was previously flagged as a major failure. Confirm nothing is
      unmerged before pruning.

## 4c. Cross-gene attribution: re-extract and republish (opened 2026-08-18)

Code is fixed (`cc86a5e`, `eb6666b`; see `docs/PROTOCOL_CHANGELOG.md`). What
remains is data. Contamination is measured per paper in
`docs/evidence/crossgene_contamination_20260818.{md,json}`: BRCA1 259/7,346
(3.5%, 15 papers), BRCA2 113/1,426 (7.9%, 6 papers), **BMPR2 0/5,838 (clean)**.

Both contaminated run dirs already carry a `QUARANTINE.md` (gitignored, so a
tracked companion lives in `docs/evidence/`); nothing was deleted and `corpus/`
source is unaffected.

- [ ] **Re-extract BRCA2 50 + BRCA1 50 on the fixed code, then BMPR2 50 as a
      control** (`results/vb_reextract_20260818/`, publication OFF). BMPR2 is
      expected to show no change; if it does, the fix is not inert on clean data
      and that is a finding in itself. All 150 papers are cached, so ~$24 and no
      fetching.
- [ ] **Review the old-vs-new diff before touching Variant Browser.** Decided
      2026-08-18: no detach, no banner, no republish until the per-paper diff is
      read. The wrong calls are live on collaborator staging until then.
- [ ] **Then decide VB isolation and replacement in one step.** PMID 26848529
      (44/124 wrong) is one of only two lead-approved BRCA2 collaborator papers,
      so the BRCA2 gold provenance needs review, not just the staged evidence.
- [ ] **Measure the four-gene headline and gold-120 Gate 2** against this change
      set. Unmeasured. The change set both closes leaks and adds rejections, so
      precision and recall can move either way.
- [ ] **Decide whether to clean stored rows in the canonical DBs.** They hold
      1,783 rows (BRCA1) and 479 rows (BRCA2) on the affected PMIDs, and
      `BRCA1.db` stores BRCA2-only `c.5291C>G` for 21232165. The fix prevents
      new occurrences but does not retroactively remove stored rows.

## 5. Engineering handoff follow-ups

- [ ] `harvesting/migrate_to_sqlite.py` repair gaps: the fabricated-c repair
      proof strips only spaces (`REPLACE(TRIM(...),' ','')`, missing
      tabs/newlines) and repairs only the first matching row (`LIMIT 1`),
      leaving later duplicates as distinct fabricated identities.
- [ ] Enable run-scoped LLM tracing for standalone extraction, discovery,
      targeted landing, figure extraction, and recall-audit entry points.
- [ ] Give `scripts/smoke_azure_models.py` a normalized, scoped connectivity
      decision event.
- [ ] Move reusable metadata-backfill, dashboard trust-reader, and EZproxy
      self-heal code out of `scripts.*` so installed-package behavior matches a
      source checkout.
- [ ] Add accepted-link tests for clinical triage, extraction-priority triage,
      claim verification, final check, source-grounded summary, and synonym
      relevance.
- [ ] Perform desktop and mobile interaction/visual QA on a completed production
      trace report in an approved browser surface.
- [ ] Extend matcher coverage for indels, isoform offsets, compound variants,
      frameshifts, splice/IVS notation, and structural/CNV events; keep structural
      rows out of precision claims until they are matchable or explicitly
      excluded.
- [ ] Investigate the residual `regex_table` count-role omissions and column
      confusion without reintroducing source-free phenotype zeros.
- [ ] Establish the Friday recall rerun/compare cadence only after the new
      explicit-zero baseline is accepted.
## Deliberate decisions and non-goals

- **One checkout, one branch.** `main` is the handoff branch. Do not create local
  feature branches or experimental worktrees unless the repository owner
  explicitly changes this policy.
- **No `unassessed_count` schema.** The partial experiment was reviewed and
  retired during the 2026-08-12 consolidation. GVF records observed carrier,
  affected, unaffected, and uncertain counts; an unobserved phenotype partition
  remains SQL NULL. It must not be inferred as a residual. A future proposal
  would need end-to-end aggregation, migration, trust/final-check/guard support,
  adjudication and Variant_Browser compatibility, provenance versioning, and a
  paired benchmark before changing the schema.
- **Count recovery stays default off** until the clean carrier-only evaluation
  passes.
- **The per-paper final check and composer stay default off**; their measured
  cost/latency did not justify default-on use.
- **Historical artifacts are immutable.** Correct them with an erratum or a new
  run; do not rewrite locked predictions, historical manifests, or append-only
  recall/protocol history.
