# Shared pipeline reliability repair — 2026-09-02

Dated evidence. Not a current-operating-plan document and not a metrics
authority. `docs/RECALL_STATUS.md` owns current metrics; `TASKS.md` owns the
plan.

This records seven shared, gene-agnostic defects found by auditing one
50-paper BMPR2 cohort, the repairs, and the measured before/after. Every
number below was reproduced against a pristine `git archive HEAD` checkout in
a scratch directory. **No historical run output was modified.**

Three independent adversarial reviews by Grok 4.6 at `xhigh` bracket the work:
one on the design before any code was written, one on the actual diff, and one
focused final verdict on the post-run evidence. All are summarized in §7.

## 1. The defects

| # | Defect | Shared? | Evidence |
| --- | --- | --- | --- |
| D1 | Archives unpacked one level only; a nested archive reached no converter | yes — `.zip` absent from the nested dispatch list *and* from the fold's convertible suffixes | one paper's bundle held an inner archive containing a 493 KB supplementary PDF (14,012 chars) with a per-patient variant table, never read by any route |
| D2 | Four independent, mutually inconsistent text normalizers | yes — router, extraction, normalizer, scanner each cleaned differently | 8 real-source tokens parsed by one route and dropped by another |
| D3 | Cross-route identity keyed on the raw notation string | yes — migration *and* every database-linkage ingest | 134 paper-groups in one 50-paper DB where a single cDNA mapped to >1 `variant_id` |
| D4 | Denominators reported from step intent, not the finalized list | yes — one shared writer | four artifacts in one run disagreed about the same 50 papers |
| D5 | Count roles never persisted; zero provenance unrecorded | yes — schema-level | 14 rows with `unaffected_count=0` and carriers>0, provenance unknowable |
| D6 | Provenance omitted the source-text policy, identity rules and archive budget | yes | two runs could share a fingerprint and read source differently |
| D7 | A re-fold could replace real converted text with a converter placeholder | yes — fold path | found live: a PDF that converted in an older environment placeholders now |

### D4 in detail — the four contradictory artifacts

One run directory stated all of the following about the same 50 papers:

| Artifact | Full text | Abstract-only |
| --- | ---: | ---: |
| per-paper extraction records (authoritative) | 50 | 0 |
| `source_completeness.json` | 0 | 50 |
| `RUN_STATUS.json` `source_integrity` | 50 | 0 |
| `workflow_summary` | 50 (`papers_downloaded=0`) | 0 |

Root cause: `download_fulltext` classifies as abstract-only anything it did not
itself fetch, and the completeness report was written from that step's return
value. On a resume path the download step contributes nothing, so every paper
was labelled abstract-only.

Second, independent bug in the same writer: the "single carrier" QA predicate
read `variant["carriers_total"]` and `variant["total_carriers"]`. Neither key
exists in the extraction schema, so `... or 1` made the expression constant and
**every paper in every run** was flagged. The artifact reported
`single_carrier_papers=50` beside a list of 30.

## 2. The repairs

1. **`utils/source_text.py`** — one canonical character policy for every route.
   Explicit closed substitution tables, deliberately **not** NFKC (which folds
   primes, superscripts and the micro sign, any of which can be meaning-bearing).
   Zero-width characters become a **space, never nothing**, so two alleles glued
   by a joiner split instead of fusing into a token that appears in no source.
   Dashes become ASCII hyphen, **never** underscore: range-versus-offset is
   notation semantics and stays behind the strict adjacency guard in
   `utils/variant_normalizer._repair_hyphen_range`.
   Wired into `pipeline/table_router`, `pipeline/extraction._clean_table_cell`,
   `utils/variant_normalizer._preprocess_cdna_token`, `utils/variant_scanner.scan`,
   and the single LLM-input funnel `_prepare_full_text`.
2. **Router grammar** — `delins` and an intronic range requiring *both*
   endpoints to carry their own explicit offset sign.
3. **Nested archives** — recursion for `.zip` only, bounded on three axes
   (depth 3, 256 MB per member, 1 GB per tree), content-digest dedup so a
   repeated or self-referential member converts once, zip-slip refusal,
   `__MACOSX`/dotfile skip. Fold side expands an on-disk archive only when its
   extraction directory is absent *or empty*, preserving the deliberate
   anti-double-count exclusion.
4. **`pipeline/variant_identity.py`** — `fold_decision` folds only a spelling
   difference in a **fully specified** identity. A field present on one side
   only is refused, so the relation is not transitive through a missing value
   and cannot bridge two rows that conflict with each other. Candidates are
   selected by shared coordinate digits, not raw string equality, so the same
   relation is applied by the writer and by the read-only detector.
   `find_identity_conflicts` produces a worklist and never merges.
5. **Typed count roles** — additive nullable columns on `penetrance_data`
   (`*_role`, `*_zero_provenance`, `count_reconciliation`). Roles are
   constrained to the closed vocabulary; an unrecognised declaration is
   `unknown`, never passed through. `sourced` requires a named column **and** a
   recognised count type. Historical zeros are **not** rewritten: a zero whose
   provenance was never recorded is a third state, and nulling it would delete
   a possibly real complete-penetrance observation while trusting it would
   report a possibly defaulted value as a finding. Reconciliation records a
   contradiction; it never balances the arithmetic.
6. **`pipeline/source_ledger.py`** — denominators from the finalized extraction
   records. Each row carries the declared class **and** a class re-derived from
   the recorded source bytes; disagreement is an explicit discrepancy, and an
   unreadable recorded source is `unverified` rather than being promoted on the
   declaration alone. No list is truncated; every count is `len()` of its list.
7. **Provenance** — the source-text policy, identity rules, ledger and archive
   code join the hashed extractor set, plus a `protocol_fingerprint` combining
   code digest, resolved config and the archive-expansion budget.

## 3. Measured results

### Cardiac gold gate — the human-curated standard

**Bit-identical before and after.** `diff` of the two reports is empty.

| Metric | Value |
| --- | ---: |
| pmids | 72/73 (98.6%) |
| variant rows | 2436/2714 (89.8%) |
| unique variants | 1585/1720 (92.2%) |
| carriers MAE | 0.393 |
| affected MAE | 0.321 |
| unaffected MAE | 0.377 |

### Offline test suite

2519 passing at the pristine baseline. This work contributes **157 tests**
across five new files plus additions to the gene-attribution suite, all
passing. The working-tree total (2692) also covers a concurrent session's
change set; see the concurrency note in §6.

### Sentinels — deterministic replay on cached source, no model calls

| PMID | bytes before | bytes after | delta | unique variant tokens |
| --- | ---: | ---: | ---: | ---: |
| 38473983 | 81,845 | 96,073 | **+14,228** | **+2** |
| 31727138 | 191,721 | 124,584 | −67,137 | 0 |
| 18823550 | 253,320 | 143,504 | −109,816 | 0 |
| 29843651 | 181,591 | 101,077 | −80,514 | 0 |
| 21801371, 28017794, 42001268, 11484688 | — | — | 0 | 0 |

The **+14,228 bytes on 38473983 is the nested-archive repair**: the pristine
baseline yields +0 bytes and +0 tokens for that paper. Its two recovered inner
tables appear exactly once each (0 → 1).

The three **shrinks are pre-existing de-duplication**, not loss: the pristine
baseline produces the identical shrink. Those papers' sources carried the same
supplement twice (a harvest-time archive block plus a fold block) and the
existing basename-keyed exclusion removes the duplicate. Verified: **zero
variant tokens lost** on all three, and zero substantive lines lost on two. The
third loses five PolyPhen/SIFT predictor-score lines from a supplement PDF that
no longer converts in this environment — also identical on baseline, and the
cause of repair D7.

### Text policy on real source tokens (baseline → now)

All eight went from **dropped** to parsed:

```
p.S775<NBSP>N            None -> p.S775N
c.1480<THIN>T>G          None -> c.1480T>G
c.722<NBSP>T>C           None -> c.722T>C
c.734<THIN>T>A           None -> c.734T>A
c.646_647delinsCA        None -> c.646_647delinsCA
c.169-198_273+820del     None -> c.169-198_273+820del
p.Arg491Trp<ZWSP>        None -> p.Arg491Trp
c.354T<FULLWIDTH>G       None -> c.354T>G
```

Corpus-wide over the eight sentinels: router-parsed tokens **711 → 712 (+1)**,
**zero regressions, zero changed normalizations**. Reported as a small real
yield: the character class is fully handled and pinned by test, but only one
instance sits in a reachable position in this particular sentinel set.

### Identity detector on the shipped 50-paper database (read-only, 522 variants)

| Class | Count | Disposition |
| --- | ---: | --- |
| foldable spelling duplicates | 81 | a current-code run no longer creates these |
| identity conflicts | 51 | human adjudication required |
| sparse partial observations | 173 | deliberately not completed |

Duplicate attribution by layer: clinvar 128, regex_table 95, llm_table 76,
pubtator 11, llm_text 6. Database linkage was the largest contributor, and
each linkage ingest carried its own copy of a raw-string `ensure_variant`.

Correctness spot-checks in that output:

- `c.1157A>G p.E386G` == `p.Glu386Gly` → **folded** (one-letter/three-letter).
- `c.1156G>A p.Glu386Gln` vs `p.Glu386Lys` → **refused**, conflict.
- `c.1156G>A p.Glu386Gln` vs `c.1156G>C p.Glu386Gln` → **refused**, conflict.
  (`Glu386Gln` is `c.1156G>C`, so one of these labellings is a source error;
  the pipeline surfaces it rather than choosing.)

### Denominator replay on the run whose artifacts contradicted each other

| Field | Before | After |
| --- | ---: | ---: |
| papers with full text | 0 | 50 |
| papers abstract-only | 50 | 0 |
| source unverified | not reported | 0 |
| class discrepancies | not reported | 0 |
| single-carrier papers | 50 (list of 30) | 26 (list of 26) |

## 4. Fresh 50-paper BMPR2 validation run

See §6. Recorded protocol fingerprint for the run:
`f6787618d6776d279f344b4151406835dbfa519a16a79f3386e19ec12453ac5f`
(`git_sha eb34c363`, `git_dirty true` — the repairs are uncommitted by
instruction, so the code digest rather than the commit identifies the protocol).

## 5. Residual uncertainty

- **One paper carries the nested-archive gain.** No second gene has been run
  through the recursion. The measurement is real but narrow.
- **The recovered inner PDF also contains a table for a different gene.** The
  cross-gene attribution guard must hold on it; this is checked in the fresh run.
- **`scripts/extract_figure_variants.py` and `scripts/recall_recovery/merge_v12_db.py`
  still raw-insert** rather than using the shared resolver. The figure path has
  its own private-to-paper guard and accounted for 8 variants in the audited DB;
  `merge_v12_db` is not on the live `gvf-run` path.
- **`_point_alias_candidate` still performs sparse/rich completion** in
  migration, guarded and pre-existing. It is order-dependent by construction,
  which is why the new spelling fold runs first and refuses sparse pairs.
- **No consumer reads the new count-role columns yet.** They are a recording
  device. Report and publish still read the raw integers, so a family count and
  an individual count remain indistinguishable to a downstream reader until
  those consumers are moved over.
- **Intronic `delins` still does not parse** (the new alternative requires the
  inserted payload at end-of-string). Yield gap, not a correctness risk.
- **`split_glued_notation` is a guard, not yet a splitter.** A cell holding two
  alleles is refused rather than parsed into two observations. Safe; yield left
  on the table.
- **BMPR2 has no gold standard.** Nothing here is tuned to an answer key, and
  no precision/recall/MAE for this cohort is defined.

## 6. Fresh-run results

`results/bmpr2_shared_fixes_20260902/BMPR2/20260902_074142`, fingerprint
`f6787618`, 50/50 papers, 85.8 min, 264 LLM calls, trace integrity
`write_time_verified` with zero integrity errors. Gold access disabled.

**Source ledger — 10/10 checks pass.** 50 finalized of 50 requested, zero
missing; 50 full text, zero abstract-only, zero unverified, zero class
discrepancies; every reported count equals the length of its list. The
independent on-disk scan agrees, and `RUN_STATUS.json` records
`agrees_with_disk_scan: true`.

**Identity worklist — clean on every safety property.** 436 variants; zero
chimeric identities; zero cross-gene rows; the `p.Glu386Gln` / `p.Glu386Lys`
contradiction still held, alongside `c.1156G>A` versus `c.1156G>C`; zero
foldable groups that are not canonical-equal on both cDNA and protein.

Same roster, prior run versus this one:

| | 2026-08-26 | 2026-09-02 |
| --- | ---: | ---: |
| variants | 522 | 436 |
| foldable spelling duplicates | 87 | 5 |
| identity conflicts | 51 | 35 |
| sparse partials | 167 | 91 |

**Archive/fold audit — every event enumerated.** 31 archives on disk, 31
expanded, zero unexpanded. Exactly one nested expansion fired. Zero zip-slip
refusals, depth-limit hits, budget exhaustions, oversized-member skips, or
re-fold refusals. Two placeholder events, both pre-existing and both verified
harmless against the pristine baseline: the existing blocks were 140 and 100
characters (themselves figure stubs), so the per-label guard correctly did not
fire, zero variant tokens were lost on either paper, and no placeholder marker
was written into any run source.

**Independent provenance on the recovered supplement.** PMC lists exactly one
supplement for that paper at 392.7 KB; the on-disk nested archive is 402,075
bytes, a byte-exact match. Its entire declared supplementary material was
therefore invisible to every route before the recursion fix. Extraction's
recorded `source_file` is the folded 96,557-byte context, and its four variants
match the published abstract with zero rows from the second gene's table inside
the same file.

**Typed count roles populated.** 283 `per_variant_carrier`, 4 `family_count`,
2 `unknown`; zero provenance 1 `sourced` and 18 `unknown`; one arithmetic
reconciliation flag. No stored integer changed.

### Two further defects the run exposed, both fixed

1. **Malformed DOI failed the paywall-fetch stage.** The queue carried
   `10.3390/ijms25052734Submission`: a rendered landing page prints the
   identifier with the next word glued on, and four modules each carried the
   same boundary-free pattern. Fixed with one canonical `utils/doi.py`; 13 real
   DOIs across this corpus's publishers are unchanged. No data was affected —
   all 50 papers already had full text, the queue held one row, and the override
   file was empty.
2. **A derived coordinate was over-splitting identities.** `genomic_position`
   was sparse-refusing in the fold predicate, so two rows with byte-identical
   cDNA *and* protein stayed separate when only one carried a coordinate. It is
   now conflict-only and the fold carries the coordinate forward. All five
   residual foldable pairs in the run are this class.

### Fingerprint gap, stated plainly

Both fixes above postdate the run, so current code is fingerprint `450ab138`
while the artifact is `f6787618`. The run's database therefore still contains
those five pairs. Its extractions were **not** re-migrated under the new code,
because pairing an extraction from one fingerprint with a migration from another
is the splice this work is required to avoid.

### Concurrent work in the tree

A separate session was editing the phenotype-count and prompt path between 08:28
and 09:07 (`pipeline/count_classifier.py`, `phenotype_count_guard.py`,
`claim_verifier.py`, `prompts.py`, `benchmarks/codex_paper_eval/run_eval.py`,
`docs/EXTRACTION_CONTRACT.md`, `scripts/refresh_run_db.py`, plus their tests and
a new `pipeline/patient_row_phenotype.py`). That work is preserved untouched.
Two consequences: the working-tree suite total spans both change sets (this
work contributes 157 tests), and the BMPR2 run executed while those files were
being written. Python binds modules at process start, so the run used the code
as of 07:41 and is internally consistent, but its recorded fingerprint must not
be re-derived from the current tree.

## 7. Adversarial review

**Round 1 (design, before any code).** Grok 4.6 `xhigh` judged the original
seven-item plan an over-build and reordered it: denominators first (the only
item that cannot invent an allele or a count), nested-zip second and scoped to
zip-in-zip, "Unicode folding is notation policy, not NFKC", identity as a
**detector not a merger**, and no recoding of historical zeros. It supplied the
eight-part conjunction that `fold_decision` implements, including the
NULL-bridge prohibition. All of that was adopted.

**Round 2 (the actual diff).** Findings and disposition:

| Finding | Verdict | Action |
| --- | --- | --- |
| Zero-width chimera still possible for **unprefixed** residue labels — the split only counted `c.`/`p.` prefixes | valid | `_BARE_NOTATION_RE` added; bare residue labels now count as allele-shaped. Pinned by test. |
| `_preprocess_cdna_token` deletes all whitespace, re-fusing what the policy split — "the single policy is not single" | valid | the normalizer now returns the unjoined value when the split guard fires, so the grammar refuses it. Pinned by test. |
| `_canon` returned `None` on parse failure, so two **unparseable** values compared equal and could fold | valid, a real false merge | `_canon` falls back to the stripped raw text; only an absent value is `None`. Pinned by test. |
| The raw-string SQL candidate pre-filter could never see protein-only spelling duplicates, making the linkage repair largely inert for its own headline case | valid and material | candidates are now selected by shared coordinate digits and compared canonically, shared by writer and detector. Pinned by test. |
| `sourced` zero provenance too cheap to earn | valid | requires a named column **and** a count type in the closed vocabulary; roles outside it record as `unknown`. |
| Declared abstract-only plus a missing file produced no discrepancy | valid | the discrepancy now covers any declared class. |
| `agrees_with_disk_scan` ignored unverified rows and discrepancies | valid | both now count against agreement. |
| Placeholder veto was all-or-nothing and could block the only measured nested-archive gain | valid, the two repairs fought | the veto is now per-label with carry-forward of the existing text. Pinned by test. |
| An empty extraction directory marked an archive permanently expanded | valid | requires a directory that actually holds files. Pinned by test. |
| Fullwidth block mapped wholesale is safe; the U+FF0D re-assert is a no-op | accepted | left as documentation. |
| Dash-family widening keeps the range/offset hazard | accepted as residual | adjacency guard is unchanged and now pinned by a non-adjacent test. |
| Still one large change rather than sliced PRs | accepted | recorded here rather than claimed otherwise. |
| Archive budget uses declared member size | accepted, low risk | per-member cap runs first. |

Round 2 verdict: **HOLD** until a fresh 50-paper run exists and its own ledger,
identity worklist and fold refusals are clean — cardiac gold being identical is
necessary but is not that dataset. That run is §6.

### Round 3 — focused final verdict on the post-run evidence

Given the run's ledger, identity worklist, archive audit, the browser-verified
supplement provenance, the diagnosed stage failure, and the fingerprint gap
stated plainly, the reviewer returned **GO** on requesting expert curation, with
the five residual pairs kept off the curator's queue as a cover-sheet item
rather than adjudication work. It judged the genomic-coordinate relaxation
correct on both attacks it was asked to attempt — one-sided absence cannot be an
assembly clash, and both-present-and-unequal still refuses, so assembly
disagreements surface as conservative conflicts rather than silent merges. It
judged the DOI trim safe for this corpus, since the rule only strips
Titlecase words glued onto a trailing digit and therefore cannot touch the
all-caps or trailing-digit suffixes these publishers use, while noting the
residual mid-suffix glue hole now tracked in `TASKS.md`.

Its operative instruction, adopted: bind any curator decision to canonical
identity (cDNA plus protein, plus the conflict pair) rather than to SQLite row
ids, because a later re-run under the current fingerprint collapses those five
pairs and moves row ids while the 35 conflicts persist.
