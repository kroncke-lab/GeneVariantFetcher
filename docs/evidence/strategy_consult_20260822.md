# Strategy consultation: recall/precision, cost, cruft (2026-08-22)

Advisory evidence only — TASKS.md owns the plan. Sources: GPT-5.6 Sol (codex
exec, xhigh, read-only sandbox), Grok 4.6 (xhigh, plan mode), and five Claude
trace/disk audits. No files were modified by the consultations, no extraction
was run, no gold values were read before the cited locks. Every number below
was either computed from disk this session or is quoted from a named authority
doc. Raw model reports and audit transcripts: session scratchpad
`strategy_consult/` (ephemeral); methods are restated here for reproduction.

## A. Findings all three sources agree on, verified on disk

1. **Figure-vision is the cost problem; claim-verify is not.**
   Trace replay of the locked `20260813_gold120_verticalfix` run
   (llm_traces/trace_index.jsonl + per-call JSON):
   - 229/244 vision calls (93.9%) added **no identity** the text extract
     didn't already have (153 returned empty, 68 all-known, 3 burned the full
     4096-token reasoning cap for nothing, 5 transport failures).
     Redundant-call cost ≈ $2.34 = 80% of the vision lane = 24% of the whole
     $9.774 pass. Replicated on `20260821_candidate_local_gold119`:
     224/243 (92.2%).
   - 47/244 calls (19.3%, $0.60) were re-reads of **byte-identical images**
     (FULL_CONTEXT lists the same figure under two URLs). SHA-256 dedup is a
     zero-risk skip: every duplicate's kept first read returned the same ids.
   - Claim-verify is NOT a no-op lane: 119/148 decisions (80.4%) changed an
     emitted value (86 nulled affected, 29 corrected unaffected, ...);
     run2: 137/188 (72.9%). Cutting verify calls ≠ cutting waste. What is
     unmeasured is whether the changes are *correct* (destroyed/killed vs
     gold) — that join is $0 from `claim_verification_field_changes`.
   - False-skip risk of a text-yield vision gate, honestly measured: novel
     ids live in chromatogram/pedigree images (19184172: R176W, R537W,
     T782fs, E788K, R1068*), so a blanket gate loses 1–4 real ids per run.
     N=5 conservative gate: run1 −$1.64 (56% of lane) losing only 2 probable
     misreads; run2 −$1.30 losing 3108insG (10086971 — likely the real gold
     frameshift) and L187P (18808722 — a cross-gene FP, losing it is a win).
     Ship SHA-dedup + token-cap retry-at-low-effort now; gate the rest on
     caption/pedigree triage, not text yield alone.

2. **The remaining Gate-2 recall gap is concentrated and mostly acquisition.**
   Census of all 88 rows in `remaining_fn.tsv` against the corpus (per-variant
   grep of every FULL_CONTEXT + supplement conversion):
   - Gold-side 13; already-landed-post-lock 11 (21956039 — rescore only).
   - **(a) never acquired: 44 rows in 7 papers**: 29650123 (20; Elsevier
     manifest advertises only the demographics docx — route: JACC publisher
     HTML online tables via EZproxy), 31520628 (7; supplements are caption
     stubs and the on-disk "supplements" are **mis-bound CDC vital-statistics
     PDFs** from refs 19/20 — route: fetch_elsevier_supplements + EZproxy),
     14678125 (6; abstract-only — route: Wiley TDM), 27114410 (6; 0-byte
     KCNQ1 FULL_CONTEXT + abstract-only fallbacks; 2 advertised supplements
     never fetched — route: AHA EZproxy/Europe PMC), 15466642 (3; AHA HTML
     without the per-case Table), 17971661 (1; Karger stub — EZproxy),
     19926015 M81L (1; Table 2 caption stub).
   - **(b) on disk but unparsed: 2** (24667783: textutil flattened a .doc
     table into delimiter-less runs — re-convert table-aware, $0).
   - **(c) extraction failure on parsed text: 18** — the variant strings are
     exact-hit in FULL_CONTEXT: 28798025 (6), 19216760 (4: exon-3 grammar +
     2 missense), 26496715 (2: 3-letter fs normalization), 21288276 (en-dash
     spaced intronic del), 28087622 (ΔKPQ — needs a generic Δ+residues
     grammar bound to nearby explicit coordinates, not an alias), 31293497
     (dup), 10086971 (prose frameshift), 22677073, 19926015 exon-3.
   - So 10 of the 27 rows labeled "residual_extraction" in the 08-15
     diagnostic are actually acquisition gaps, and the $0 pool
     (11 landed + 2 + 18) is **31 rows ≈ +4.9 pts** before any fetch.
     Caveat: counts are vs the 08-15 diagnostic; Grok's read of the 08-21
     lock says several papers already improved (26496715 40TP/2FN,
     21956039 16/16) — regenerate the taxonomy from the new lock first.
   - Path to 92% (583/633 needs +37): 11 landed + 2 (b) + half of (c) +
     ~15 of the 44 fetchable ≈ there, with margin from the rest.

3. **Do not spend on Gate-2 counted-extra precision.** 98.03% with a
   curator-only residue. Raw precision moves only by expanding gold on
   catalogue papers (26746457 first). Unchanged from PRECISION_COST_LEVERS.

## B. Count coverage/value findings (verified)

- Coverage ceilings under no-inference (gold-0 arithmetic on the 08-20
  noinference report): affected ≤64.3%, **unaffected ≤13.0%**, carriers
  ≤74.7%. Unaffected 5.06% is mostly the null-vs-0 convention; actionable
  numbers are non-zero recall: carriers 40.0%, affected 10.1%, una 23.2%.
- Per-layer emission on the same lock: regex_table rows emit counts 20.6%
  carriers but **0.4% affected** — `extraction.py::_table_variant` stores
  `source_table_headers` + `evidence_quote` but never binds count/phenotype
  cells. The identity-only recovery layers (clinvar 294, regex_text 319)
  have no count slot at all; cross-layer fill-only inheritance at merge is
  the $0 mechanism.
- Three highest-yield generic mechanisms, all measurable by parse-only
  replay over cached gold-120 source + `phenotype_value_precision.py`:
  (1) finish the count/phenotype column binder incl. vertical-PDF-text
  layout (21956039) and the regex-sweep cell binding; (2) one-individual-
  per-row roster binder + whole-cell phenotype lexicon (30403697 alone: 21
  carriers/21 affected/8 unaffected explicit); (3) cross-layer count
  inheritance at merge.
- Guard calibration, from AFFECTED_UNAFFECTED_PRECISION §2 arithmetic:
  affected guard 0.92 destroyed/killed vs 3.64 budget (fine); **unaffected
  1.67 vs 2.08 — nearly break-even, and the ~10 destroyed exacts are ~31%
  of achievable coverage**; the N=1 `figure_copied_phenotype` clear was
  never separately measured. `figure_count_unverified` trust mask projects
  away fills measured **90/91 exact** (`pipeline/count_repair.py`) —
  a promote-on-corroboration path restores ~90 correct counts for $0
  extraction. `carrier_guard` nulls affected/unaffected/uncertain alongside
  the offending >100k carriers value; fold into trust record instead.
- **Verified bug (Sol found it):** `pipeline/extraction.py:7583` constructs
  `VariantClaimVerifier(..., reasoning_effort=self.reasoning_effort)` — the
  Grok extractor's effort — while the adjudication path one screen up
  explicitly swaps to `adjudicator_reasoning_effort` because inheriting the
  primary effort "was a silent mismatch". Same mismatch, unfixed, on the
  layer that is ~37% of spend.
- Before any verify-after-merge or model swap: freeze
  `claim_verifier._apply_consistency_guards` (the affected=total/una=0
  trap) and measure verify destroyed/killed from stored field_changes.

## C. Cost plan (evidence-backed, cumulative on the $9.77 pass)

| Lever | Save | Risk |
| --- | --- | --- |
| SHA-dedup images per paper | ~$0.60 | zero (measured) |
| Token-cap vision retry at low effort / accept-empty | ~$0.7–1.0 | zero |
| Caption/pedigree-triaged vision gate (not text-yield-only) | ~$1.0–1.6 | 1–2 ids/run, measure first |
| Prompt caching: lift `_CORE_RULES`+schema into a constant system msg (call sites already accept `system_message` and don't pass one); reorder claim-verify schema above the card (crosses Azure's 1,024-tok floor) | ~$0.6–1.0 (6–10%), ~$2–3.5 on BRCA-shaped runs; also makes CONTINUATION re-sends cache hits | prompt-order QA: task line stays after text |
| Batch claim cards per paper (1 call, N verdict keys) | small on gold-120 (1.3 calls/paper), big on BRCA runs (785→~130 calls) | needs card-id keyed verdicts |
| Pedigree *detection* to cheap model/deterministic pre-filter (currently a full Sol vision call at max_tokens=200 per image) | share of vision lane | low |
| Luna for claim-verify | large | gated: consistency-guard trap + gold-120 rescore first |
| `TABLE_ROUTER_REASONING_EFFORT=low` | ~$0.1 | trivial |

Also fix: the raw Responses-API vision path (`figure_variant_reader.py:433`,
same in figure_text_extractor) captures **no usage** — 0 tokens in every
08-20/21 `figure_variant_read` trace, so the $9.77 proxy is a floor for the
vision slice. Fix capture before trusting any downgrade measurement.

## D. Cruft (method + residue)

Prior executed audit: `docs/evidence/cruft_dead_code_audit_20260820.md`.
Residue found this session:
- DELETE-SAFE (zero references anywhere, superseded): `pipeline/failure_logger.py`
  (+ its test), `scripts/recall_recovery/recover_paywall_oa.py` (violates the
  catalog's own one-section rule).
- ASK/QUARANTINE clusters: deprecated replay pair (replay_cap_trip/
  retry_failed — window condition), adjudication-packet cluster (5 scripts +
  somatic_germline_qc; zero runbook refs but serves no-gold cold starts),
  stage_corpus_figures + fetch_gold_figures (zero-lift experiment),
  gold_source_worklist (superseded by build_acquisition_worklist),
  dedup_db, 9 recall_audit one-off pilots, cli/fetch_manager.py,
  evidence-card pair (dormant B2).
- Import-graph traps confirmed KEEP: run_priority_extraction (shelled out by
  full_coverage), merge_v12_db (shelled by run_all_layers), everything in
  harvesting/ (alive via relative imports).
- Deletion contract (Sol+Grok convergent): reachability from entry points +
  subprocess/doc/runbook refs + flag owner + provenance-hash membership;
  tombstone window; one subsystem per PR with wheel build + `gvf --help`
  snapshot + offline/negative suites. Never delete parked science
  (paper_final_check, count_recovery, count_classifier).

## E. Beliefs to falsify next (both models converged)

1. "Usable-source extraction is ~100%" — false as stated: 28798025 has 6
   exact-hit variants unextracted (largest true extraction failure).
2. "The 08-15 FN taxonomy is current" — it is stale; rebuild from the 08-21
   lock artifacts ($0) before ranking extraction work.
3. "Verification buys its precision" — changed≠correct; join field_changes
   to gold-matched exactness ($0).
4. The four-gene "85% acquisition" split cannot be regenerated from its own
   artifact (`disagreement_artifacts_skipped: true`); regenerate before
   inheriting it into Gate-2 strategy.
