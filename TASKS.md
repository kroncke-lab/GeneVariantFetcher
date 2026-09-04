# GVF Handoff Tasks

Last reviewed: 2026-09-03.

This is the only active GVF checklist. Current measurements and caveats live in
[`docs/RECALL_STATUS.md`](docs/RECALL_STATUS.md); completed benchmark history
lives in [`docs/RECALL_HISTORY.md`](docs/RECALL_HISTORY.md); protocol changes
are append-only in
[`docs/PROTOCOL_CHANGELOG.md`](docs/PROTOCOL_CHANGELOG.md). How to raise
gold-120 **identity** precision and cut cost (which denominator, which *lane*,
the 12 counted extras, Sol/Grok/Kimi shares, $0 vs paid) lives in
[`docs/PRECISION_COST_LEVERS.md`](docs/PRECISION_COST_LEVERS.md). How to raise
**affected/unaffected value** exact-match precision (not Gate 2; no one-off
aliases) lives in
[`docs/AFFECTED_UNAFFECTED_PRECISION.md`](docs/AFFECTED_UNAFFECTED_PRECISION.md).
Dated reports and benchmark run directories are evidence, not competing plans.

Execute the acceptance gates in order. Do not publish a new headline, update the
dashboard, enable a default-off recovery stage, or advance to a larger cohort
until the preceding gate passes and its artifact is recorded.

## Active $100 improvement goal

- [ ] **Mixed-gold tranche campaign (started 2026-09-03; Brett's /goal: 3-4
      tranches, recall up and MAE down, papers only, grok 4.6 + agy advice).**
      **Cohort-size correction 2026-09-03:** Brett specified that each tranche
      should be approximately the 118-attempt Gate-2 size. The already-scored
      49-attempt tranches 01–02 remain immutable calibration evidence; unopened
      03–29 are abandoned in `mixed_gold/abandonment_log.jsonl`. The replacement
      `mixed_gold_continuation_120/` excludes the 90 consumed PMIDs and partitions
      the remaining 1,324 attempts into 11 article-atomic tranches of 120–121.
      Arm schedule (both reviewers agree): every tranche pairs baseline =
      frozen `506a949c` vs candidate = cumulative fixes; confirmation is only
      the *identical* candidate runtime on the next unopened tranche.
      Acquisition ceiling measured first (`docs/evidence/
      gold_source_presence_sweep_20260903.md`: 15.8% hard / 28.7% wide of
      runnable gold rows; supplements are NOT where the missing gold is).
      Tranche 01 baseline locked: paper-derived 155/61/87 (recall 64.0%,
      precision 71.8%); 70/87 misses acquisition, 5 reachable
      (`scripts/recall_audit/fn_root_cause.py`). Candidate `b56f469f` (six
      gene-agnostic fixes, `docs/evidence/mixed_gold_tranche01_20260903.md`)
      scored: 157/54/85, recall +0.83 pp (**below the +1.0 pp discovery
      bar → `reject_or_revise_candidate`**), precision +2.65 pp,
      counted-extra precision 82→91%, carriers supplied 48→125 (conditional
      MAE 0.812→0.104), affected 101→81 (guard interaction, fixed in v2).
      Tranche 01 cost by the public-price proxy: baseline $4.369 + candidate
      $4.206 = **$8.575** (registry estimate was $10.25 paired); ledger
      $44.519 + $8.575 = **$53.094 used / $46.906 remaining** before tranche
      02; tranche 02 baseline $4.555 (1.124 M tokens) → **$57.649 used /
      $42.351 remaining** before the tranche 02 candidate. Tranche 02 = v2 discovery: baseline 268/54/35
      (recall 88.45%), candidate 267/55/36 (−0.33 pp;
      `reject_or_revise_candidate`; secondary count rule also not passed).
      **Every arm difference was model run-to-run variance** (protein
      notation or carrier value present in one arm's extraction row and
      absent in the other's, temperature 0) — a paired tranche with ~5
      reachable rows cannot resolve a reading change from provider
      nondeterminism (`docs/evidence/mixed_gold_tranche02_20260903.md`).
      Tranche 01's secondary endpoint as a locked diagnostic only: e2e
      carrier MAE 2.68 → 1.95, coverage 31% → 80%, all criteria held.
      Ledger after tranche 02 (candidate $4.406): **$62.055 used / $37.945
      remaining**. The 49-attempt tranche 03 was NOT opened: per the
      cohort-size correction above, the campaign continues on
      `mixed_gold_continuation_120/` tranche 01 (120 attempts, 384 gold rows,
      hard ceiling 87 rows) as a v2 discovery under both rules; baseline arm
      `runs/20260903_protocol_cont120_01_baseline` locked: **263/138/121**,
      recall 68.49%, precision 65.59%, carriers 123/384 (MAE 0.236), $10.339
      (2.70 M tokens, 68 min) → ledger **$72.394 used / $27.606 remaining**;
      candidate arm (v2 runtime, commit 8cf78b39 after restoring the nine
      files that a docs commit made while frozen had swept back to
      `506a949c`): **261/141/123**, recall −0.52 pp, precision −0.66 pp
      (guard fails at the bound), carriers 123→130, e2e carrier MAE
      1.393→1.359 (UB +0.08); both rules reject. $11.641 → **final ledger
      $84.035 used / $15.965 remaining**; no further arm fits. Campaign
      summary and ranked next levers (acquisition worklist first; VF
      false-positive gate and linkage row ownership; source-only
      substitution promotion; prose exclusion lists; bare table cDNA;
      replicate arms before scoring again):
      `docs/evidence/mixed_gold_campaign_summary_20260903.md`.
      `docs/evidence/mixed_gold_cont120_01_20260903.md`.
      **Secondary count endpoint preregistered 2026-09-03** (before the
      tranche 02 candidate lock; commits a47e3a5a + 8ef743a4;
      `docs/evidence/mixed_gold_count_endpoint_preregistration_20260903.md`):
      end-to-end carrier MAE (miss/abstention = full error) with observed
      delta ≤ −0.05 and bootstrap upper bound < 0, coverage-on-matched
      non-decreasing, and a hard identity guard (recall LB ≥ −1 pp, precision
      LB ≥ −2 pp, observed recall delta ≥ 0). Identity stays primary;
      tranche 01/02 identity verdicts stand. Rule parameters go in
      `mixed_gold/secondary_endpoints.json` (sibling of the registry, whose
      digest every scored arm binds). Both reviewers (agy, grok) judged this
      legitimate given the measured acquisition ceiling, provided it never
      re-adjudicates burned tranches and never loosens the identity rule.
- [x] **Close the repeated source-positive phenotype-count loss across all
      stages and rerun the full cardiac cohort (2026-09-03).** Code-owned
      provenance is now non-forgeable, refresh preserves it, exact patient-row
      and source-bound partitions survive verification/trust, and dedicated
      carrier denominators outrank the legacy patient mirror. Fresh gold-blind
      lock `20260902_false_zero_recovery_gold118` completed 118/118 attempts;
      required cases survived as W4645R 4/2/2, C2277R 8/7/1, G357S
      185/97/62, P2328S 62/17/42, and Y652X 6/6/NULL. The figure and its raw
      CSV/JSON are generated automatically at score time. Full suite: 2,736
      passed. The run validates the repair but is not promoted because the
      fresh stochastic extraction has lower identity recall and carrier/affected
      supply than the accepted headline.
- [x] Consolidate every worthwhile local branch/worktree change onto `main`;
      preserve the stale pre-merge stash as historical recovery material.
- [x] Freeze the accepted 118-attempt headline and run two new gold-blind,
      pre-lock validations. Neither passed all promotion gates; do not select
      the best metric from different stochastic runs.
- [x] Merge the source-backed scanner/table/structural fixes and targeted-refresh
      usage telemetry. Full result:
      `docs/evidence/gold118_source_recovery_lock_20260824.md`.
- [x] Repair AHA collapsed cohort tables, bind the upgraded PMID 15466642 source,
      and complete one final sealed validation. The source-bound lock passed
      every promotion gate and is now authoritative: 554 TP / 283 FP / 78 FN,
      87.66% recall, 66.19% raw precision, 94.59% count-bearing-only precision,
      and carrier MAE 0.193 over 207 supplied rows. Full result:
      `docs/evidence/gold118_aha_table_lock_20260824.md`.
- [x] Preregister the exact BRCA1/BRCA2/BMPR2 150-paper evaluation as 30
      calibration + 20 holdout papers per gene, with frozen source hashes,
      exhaustive-paper precision, explicit `NONE`, and family-count exclusion.
- [x] **Close the gold-150 cross-gene holdout firewall (2026-08-25).** The
      original split ranked `sha256("<seed>|<GENE>|<PMID>")`, so a PMID shared by
      two genes was ranked independently under each. Six PMIDs were calibration
      for one gene and holdout for another, putting **5 of the 20 BRCA1 holdout
      papers into BRCA2 calibration**; seven of the eight cross-gene PMIDs were
      also bound to different source bytes per gene. Fixed before any answer key
      existed. **Curate from
      `benchmarks/curated_extraction_eval/gold150_preregistered_20260825_amended/`**;
      the 20260824 tree is preserved and marked superseded. Verify with
      `scripts/audit_split_firewall.py` before curation and again before holdout
      release. Known residual: the holdout roster is tracked in git, so this
      cohort is answer-key blind, not roster-sequestered.
- [~] **DESCOPED 2026-08-25 (Brett): do not curate the BRCA1/BRCA2/BMPR2 150.**
      The 90-calibration + 60-holdout human curation effort is not being run.
      The cardiac **gold-120** set is the human-curated standard to work
      against. The amended packets and the firewall audit remain on disk as
      correct, parked artifacts; the split-lock/canonical-source support in
      `build_curation_packet.py` and `scripts/audit_split_firewall.py` are kept
      because they are general packet infrastructure, not 150-specific. Do not
      spend human curation time here without Brett reopening it.
      Exact-150 P/R/F-score/MAE stay **undefined**, not zero.
- [ ] Obtain the missing KCNH2 PMID 29650123 mutation roster from an author or
      other provenance-valid source. The publisher's only supplement has no
      roster; never infer its 20 missing singleton identities from gold.
- [ ] **Remove gold-derived data from the runtime path.**
      `gvf_data/kcnh2_variant_aliases.json` declares
      `"source": "Gold standard Excel + generated forms"` (856 variants / 4,766
      aliases) and is loaded at runtime by `_lookup_alias`
      (`utils/variant_normalizer.py:1588`) for every gene. KCNH2 is a gold-120
      gene, so KCNH2 normalization has access to the answer key. Re-matching all
      554 locked pairs with the lookup stubbed out gives 554/554, so the current
      score is **not** manufactured by score-time lookup — but extraction-time
      influence cannot be recovered from normalized predictions. Rebuild the map
      from public notation resources, or gate it out of production and keep it as
      a benchmark-only input. Details:
      `docs/evidence/generalization_consult_20260825.md` §1.
      **Evaluation containment landed 2026-09-03:** `--gold-free-run` now disables
      all file-backed alias maps (even a warm in-process cache), records that in
      `RUN_STATUS.json`, and the paper-primary projector refuses a run without
      that proof. Rebuilding/removing the map for ordinary production remains
      open.
- [x] **Report the provenance lanes separately.** New registry-created runs lock
      paper-derived identities as the primary score, ClinVar/PubTator as an
      external audit plus secondary linkage-assisted diagnostic, and ambiguous
      `mixed`/legacy/unknown origins as an unscored attribution-quality lane.
      The historical headline mixes paper
      extraction with ClinVar/PubTator linkage. Locked-DB replay: paper-derived
      only is 512/125/120 (P 80.38%, R 81.01%) versus 554/283/78 (P 66.19%,
      R 87.66%) linkage-assisted. Linkage FPs are KCNQ1 79, SCN5A 72, KCNH2 7,
      RYR2 0 — which is why RYR2 looks clean, and why "the top five FP papers are
      just gold-incomplete" is too strong (95 of those 162 extras are linkage
      rows, including all 63 on KCNQ1 19632626). Do not delete linkage rows.
- [x] **Mix the complete named-variant gold inventory into protocol-regression
      tranches.** `benchmarks/evaluation_tiers/mixed_gold/` inventories 1,534
      gene-paper attempts: 1,422 source-available attempts appear exactly once
      in 29 deterministic, article-atomic tranches; 111 unavailable and one
      quarantined wrong-paper attempt remain explicit. Per-tranche one-arm cost
      is $4.98–$5.30; paired baseline/candidate cost is $9.96–$10.61, or
      $12.45–$13.26 with headroom. No paid extraction was run. Consume in
      registry order; any inspected score burns that tranche for confirmation.
      Setup/scoring enforce baseline-before-candidate, one score per arm, and
      append-only consumption. The paired decision is preregistered and applied
      by `run_eval.py compare`; scoring also verifies the exact setup-pinned
      composite gold digest, locks observed gold-free RUN_STATUS evidence, and
      refuses a registry-managed re-lock without its burn contract. Reviews and
      dispositions: `docs/evidence/mixed_gold_grok_review_20260903.md` and
      `docs/evidence/mixed_gold_claude_review_20260903.md`. Full unit suite:
      **2,751 passed**.
- [ ] **Ship the ranked $0 generality fixes** in
      `docs/evidence/generalization_consult_20260825.md` §8, in order: HGVS
      range-separator repair (`c.2550-2551insTG` vs `c.2550_2551insTG` currently
      does not match), malformed-protein-range refusal, reference-backed 3'
      shift (predicted 555/279/77), the provenance-lane scorecard, and end-to-end
      count loss in place of conditional MAE as the decision metric. Neither
      external reviewer authorized another paid gold-118 extraction.

Budget: **$44.51906260 used / $55.48093740 remaining**. The final AHA
source-bound lock cost $11.24050705 by the established public-price proxy and
is promoted. No further paid gold-118 best-replicate run is authorized. The
remaining envelope is reserved for a calibration-informed exact-150 candidate
or one post-freeze holdout extraction if the existing frozen DB cannot serve as
the candidate.

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

- [x] **Act on the 2026-08-24 blind proof result, then stage the reviewer
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

      **2026-08-24 follow-through:** source adjudication did not support masking
      either dominant KCNH2 value to improve the benchmark: PMID 16029385
      supports the extracted 16-person clinical group while the fixture says 9,
      and PMID 22338672 supports seven K897T torsades carriers while the fixture
      says zero. Preserve both as answer-key review items and score any corrected
      gold as a separately versioned snapshot. A generic grouped-header repair
      was instead replayed on SCN5A PMID 29709101; together with structural-only
      projection it yields the source-backed diagnostic in `RECALL_STATUS.md`
      (548/632 recall matches; carrier MAE 0.226) without mutating the lock.
      The fresh grouped/structural and source-recovery locks are now complete.
      The latter restored accepted recall and lowered carrier MAE to 0.198, but
      failed carrier coverage (197 < 206) and count-bearing-only precision
      (93.396% < 93.694%), so neither lock was promoted. Targeted refresh usage
      telemetry is fixed. The attributable ledger is $33.0688 used / $66.9312
      remaining; do not spend it on another best-replicate attempt.
      Post-lock acquisition diagnostics also re-folded four run-local source
      files (RYR2 19216760/22222782/34661651 and SCN5A 29544605); the immutable
      predictions/report hashes are unchanged, but the source-binding check now
      correctly detects those files as drifted. Restore the exact selection-
      hashed snapshots from archival storage, or create a clean scaffold, before
      treating this run directory as a source-valid parent for another lock.

      **Final source-bound follow-through:** the clean replacement scaffold
      bound the repaired AHA source and completed 118/118 gold-free extractions
      before lock. Its 554/283/78 result passes every identity, precision,
      coverage, and MAE floor and supersedes the postfix headline. The final
      ledger is $44.51906260 used / $55.48093740 remaining; do not spend it on
      another gold-118 replicate. Continue with human exact-150 calibration.
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

      The set is not a complete gold standard: only 3/150 papers overlap the
      approved curated fixture and BMPR2 contributes zero. Before any claim of
      150-paper precision, recall, or MAE, seal a human source-adjudicated
      calibration/holdout design (recommended 60/40 within gene), report each
      gene separately, and retain a never-tuned holdout. LLM review may prepare
      evidence cards but may not define the answer key.
- [ ] After an accepted rescore, update `docs/RECALL_STATUS.md`, append
      `docs/RECALL_HISTORY.md` and `docs/PROTOCOL_CHANGELOG.md`, and regenerate
      the dashboard. Until then, the public dashboard remains an archived
      pre-correction snapshot.
- [ ] **Gold-120 refinement follow-through.** Ranking and “do not”
      constraints: [`docs/PRECISION_COST_LEVERS.md`](docs/PRECISION_COST_LEVERS.md).
      Paper-level evidence:
      `benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/`
      (`PRECISION_AND_COST.md`, `NOTES.md`, `remaining_fn.tsv`) — that
      2026-08-15 diagnostic (545/633, raw 40.85%, counted-extra 98.55%, carrier
      MAE 0.292) is **superseded** and must not be quoted as current. The
      authoritative lock is `20260824_aha_table_sourcebound_gold118`: **554 TP /
      283 FP / 78 FN**, recall 87.658%, raw precision 66.189%, counted-extra
      **97.880% (12 counted extras)**, conventional count-bearing 94.595%,
      carrier MAE 0.193. Read
      [`docs/evidence/generalization_consult_20260825.md`](docs/evidence/generalization_consult_20260825.md)
      before proposing any change here: gold-118 is now a calibration set, and
      the headline mixes paper extraction with ClinVar/PubTator linkage
      (paper-only lane is 512/125/120, P 80.38%).
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
         Do **not** delete the identity-only extras to raise raw precision —
         on the current lock **237 of the 283 extras (83.7%)** sit on papers
         that already matched every gold row, and 95 of the 162 on the five
         largest are ClinVar linkage rows rather than paper extractions. If raw
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
- [ ] Decide a policy for the ~15 giant legitimate FULL_CONTEXT folds found by
      the 2026-08-31 corpus audit (largest: APOE 10599054 and 15538542 at
      4.28 GB each — identical sizes, likely the same runaway supplement —
      then LMNA 39468056 505 MB, SCN5A 27063795 412 MB, LMNA 29970176 379 MB,
      MYBPC3 33757590 325 MB; all >2 MB files were scanned and these are NOT
      binary garbage, so `scripts/strip_binary_garbage_blocks.py` must not
      touch them). They predate the new `GVF_FOLD_MAX_TOTAL_CHARS` cap, which
      only prevents new ones. Options: re-fold under the cap (drops real rows
      — needs a per-paper look at what the oversized supplement is), a
      gene-scoped trim, or leave-and-rely-on-the-scout-budget. Do not delete
      or blind-truncate hard-won source.

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

## 4a. Shared pipeline reliability (landed 2026-09-02)

Seven shared, gene-agnostic defects found by auditing one 50-paper BMPR2 cohort.
Full evidence, measurements, and three Grok 4.6 `xhigh` adversarial reviews:
`docs/evidence/shared_pipeline_reliability_20260902.md`. The cardiac gold gate
is bit-identical before and after. This work adds 157 tests; the current working
tree passes all 2692 tests, including a concurrent change set.

- [x] **One canonical source-text policy** (`utils/source_text.py`) used by the
      table router, extraction, the scanner, the normalizer, and the LLM-input
      funnel. Explicit closed substitution tables, deliberately not NFKC. Eight
      real-source tokens went from dropped to parsed; corpus-wide over the eight
      sentinels 711 → 712 router-parsed tokens with zero regressions.
- [x] **Nested archive expansion** (`.zip` only, depth/size/cycle bounded,
      zip-slip refused). Recovered a 14,012-char supplementary PDF holding a
      per-patient variant table that had never reached any extraction route.
      Verified live: that paper's run source grew 82,211 → 96,557 bytes with the
      recovered tables present exactly once.
- [x] **Cross-route identity as a detector plus a spelling-only fold**
      (`pipeline/variant_identity.py`). Candidates are selected by shared
      coordinate digits, so the writer and the read-only detector implement the
      same relation. A field present on one side only is refused, so the
      relation cannot bridge two conflicting rows through a missing value.
      Every database-linkage ingest now uses the shared resolver instead of its
      own raw-string helper.
- [x] **Finalized-list denominators** (`pipeline/source_ledger.py`). Also fixed
      the single-carrier predicate, which read two non-existent schema keys and
      therefore flagged every paper in every run.
- [x] **Typed count roles and zero provenance** persisted on `penetrance_data`,
      additive and nullable. Historical zeros are not rewritten.
- [x] **Protocol fingerprint** covering the text policy, identity rules, ledger
      and archive budget.
- [x] **Re-fold no longer trades real converted text for a converter
      placeholder**, per label rather than all-or-nothing.

A gold-free 50-paper BMPR2 validation run was completed on this protocol
(`results/bmpr2_shared_fixes_20260902/BMPR2/20260902_074142`, fingerprint
`f6787618`). Its own artifacts validate: source ledger 10/10 checks, 50/50 full
text, zero unverified, zero class discrepancies, zero chimeric identities, zero
cross-gene rows, 31/31 archives expanded, and the known source contradictions
still held rather than merged. Foldable spelling duplicates over the same
roster fell 87 → 5. Numbers and caveats: `docs/RECALL_STATUS.md`.

- [x] **Make the count roles load-bearing.** The trust gate now reads the
      persisted role column, falling back to the JSON provenance blob, so a
      declared family count is masked from carrier-facing totals even when that
      blob is missing. Previously the role reached the gate only through a LEFT
      JOIN to `variant_papers`, and a missing or malformed blob silently made a
      family count read as an individual count.
- [x] **One canonical DOI path** (`utils/doi.py`). Four modules each carried the
      same boundary-free pattern, so a DOI printed in rendered page text with
      the next word glued on ("...25052734Submission") became an identifier that
      cannot resolve; that failed the validation run's paywall-fetch stage. The
      trim is structural (capitalised words glued onto a trailing digit), not an
      allowlist, and 13 real DOIs across this corpus's publishers are unchanged.
- [x] **Stop over-splitting on a derived coordinate.** `genomic_position` was
      sparse-refusing in the fold predicate, so two rows with byte-identical
      cDNA *and* protein stayed separate when only one carried a coordinate. It
      is now conflict-only, and the fold carries the coordinate forward. This is
      the class of all five residual foldable pairs in the validation run.

Remaining, in order:

- [ ] **Move the remaining count consumers onto the role columns.** The trust
      gate now reads them, but report, publish, and the scorer still read the
      raw integers, so a family count and an individual count remain
      indistinguishable to a downstream reader.
- [ ] **Re-run the 50-paper BMPR2 cohort under the current fingerprint** if a
      zero-residual artifact is wanted. The coordinate and DOI fixes postdate
      the validation run, so its database still contains the five characterized
      foldable pairs. Extraction output is unaffected by either fix, so this is
      a presentation choice rather than a correctness gate; do not splice a new
      migration onto that run's extractions.
      Command:
      `gvf gvf-run BMPR2 --email <addr> --output ./results/<new> --pmid-file benchmarks/curated_extraction_eval/review_pmids_50_20260824_curated/BMPR2.txt --disease "pulmonary arterial hypertension" --gold-free-run`
- [ ] **Close the mid-suffix DOI glue hole.** The trim is end-anchored, so
      `...34Submissionreceived:6` would still pass through as a bad identifier.
      Fetch-queue robustness only; it does not affect bound source.
- [ ] **Run the nested-archive path on a second gene.** The recall gain is
      measured on exactly one paper. Until a second gene exercises it, the
      repair is correct-by-test but narrow-by-evidence.
- [ ] **Route `scripts/extract_figure_variants.py` and
      `scripts/recall_recovery/merge_v12_db.py` through the shared resolver.**
      Both still raw-insert. The figure path has its own private-to-paper guard
      and contributed 8 variants in the audited database; `merge_v12_db` is off
      the live `gvf-run` path. Neither is urgent; both are identity writers that
      have never seen `fold_decision`.
- [ ] **Work the 51 identity conflicts** the detector surfaced on the audited
      database. These are genuine source contradictions (a nucleotide change
      reported with two different protein consequences, a consequence reported
      against two different nucleotide changes). They need curator adjudication;
      the pipeline must not choose.
- [ ] Retire `_point_alias_candidate`'s sparse/rich completion, or prove it
      safe. It is pre-existing and guarded, but it is the one remaining
      order-dependent fold in migration.
- [ ] Parse intronic `delins`; teach `split_glued_notation`'s callers to emit
      two observations rather than refusing a two-allele cell. Both are yield
      gaps, not correctness risks.

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

- [x] Make phenotype-zero claims source-closed and paper-target aware, then run
      the complete 118-paper blinded cohort. `20260902_zero_claim_guard_gold118`
      completed and locked; 2,726 tests pass. D1790G is 30/30/NULL and G400A is
      1/1/NULL under the shared protocol. The automated recovery figure is
      generated by `lock_and_score.sh`; its CSV preserves raw null separately
      from evaluation-zero filling. The only stored-zero/positive-gold mismatch,
      V141M PMID 28491547, is a documented legacy-gold definition conflict.
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
