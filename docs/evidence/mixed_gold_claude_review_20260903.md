# Mixed-gold harness: Claude Opus 5 `xhigh` release review

Date: 2026-09-03

This was an independent, read-only repository review after the Grok 4.6 design
review. The local Claude CLI ran Claude Opus 5 with `xhigh` effort and access
only to repository reads/searches. It made no edits, web searches, paid GVF
extractions, or publications. The CLI reported 56 turns, 54,186 output tokens
(32,496 thinking tokens), and a $6.347757 list-price account charge for the
review itself.

## Verified controls

- ClinVar/PubTator-origin identities cannot enter the paper-derived primary
  lane; the stable origin token, rather than later corroborating witnesses,
  selects the lane.
- `mixed`, `manual_or_legacy`, empty, and unknown origins are held in the
  unattributed lane and excluded from both scores.
- Gold-free extraction disables the file-backed alias loader before consulting
  its in-process cache.
- The 1,422 runnable attempts reconcile exactly to 28 tranches of 49 and one of
  50; all 1,534 inventory attempts reconcile to 1,422 included, 111 unavailable,
  and one quarantine. PMID assignment is article-atomic.
- Paired arm order, one-score-per-arm burn semantics, and the cost arithmetic
  were independently recomputed. Suite estimates reconcile to $148.832 one-arm,
  $297.664 paired, and $372.08 paired with 25% headroom.

## Findings and dispositions

| Finding | Disposition |
|---|---|
| No registered paired decision margin or executable comparison | Fixed before any run. The registry now preregisters candidate-minus-baseline deltas, a +1 percentage-point observed recall threshold, one-sided 95% PMID-cluster-bootstrap lower bounds of -1 point for recall and -2 points for precision, 10,000 deterministic resamples, and identical confirmation on the next unopened tranche. `run_eval.py compare` applies the rule and refuses mismatched tiers, arms, lanes, attempts, gold denominators, or registry digests. Each arm binds the active registry in setup and its burn event. These are engineering thresholds, not clinical-effect thresholds. |
| `score` accepted an arbitrary `--gold-root` and could silently fall back | Fixed before any run. `setup.json` now carries the registry's answer-key digest; scoring requires the exact resolved root and recomputes its composite digest. The generated finalize script is also audited for `--paper-primary`, the pinned gold root, and projection-lock-score ordering. |
| Gold-free status recorded CLI intent rather than observed alias-gate state | Fixed. `RUN_STATUS.json` records the observed environment gate; paper-primary projection requires both gold access and aliases disabled. Each production status file and digest is carried into predictions, bound into `LOCK.json`, and revalidated at score. |
| Deleting `setup.json` and re-locking could bypass the tranche ledger | Fixed. A manifest registered under a suite with a consumption log cannot be locked or scored without its setup/burn contract. |
| Per-tranche gold denominators vary (110-593 rows), and rare provenance strata are sparse | Accepted limitation, disclosed in protocol docs. Paired inference remains within the same tranche. Per-provenance numbers are diagnostic with their sample sizes and are not a balanced cross-provenance claim. |
| Dropping a candidate after scoring only its baseline can block progression | Open operational follow-up. Until an explicit append-only abandonment event exists, do not open a baseline arm unless its candidate arm will be completed; never hand-edit the ledger. This does not block the first planned pair. |
| Lane-blind duplicate-notation selection is a latent risk | Open lower-priority follow-up. Claude found zero occurrences in the audited production DBs; the primary/linkage split remains valid for the prepared suite. |
| Non-exhaustive gold makes absolute pooled precision hard to interpret | Accepted limitation. The paired delta remains useful, but absolute pooled precision must not be presented as an exhaustive scientific estimate. |

## Verdict

Claude's initial verdict was conditional GO for `mixed_gold_01` after the first
two blockers. Both blockers landed before any extraction or score, along with
the status-lock and ledger-bypass hardening. The first paired measurement is
therefore ready from a harness perspective. The 111 unavailable sources and
the heavily cardiac, heterogeneous gold composition remain explicit limits;
"mixed" describes the workload, not a balanced provenance cohort.
