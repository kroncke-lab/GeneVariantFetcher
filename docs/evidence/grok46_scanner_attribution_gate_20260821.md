# Grok 4.6 review and scanner attribution gate — 2026-08-21

## Incident

The first production extraction of replacement BRCA2 PMID 26843898 showed a
general protocol failure. The primary full-text extractor returned the seven
BRCA2 variants observed in the study and explicitly excluded literature
overview tables. `merge_scanner_results()` then appended 48 count-null regex
identities, including prior-study BRCA2 rows and BRCA1/PALB2 variants, while
stamping every row `gene_symbol=BRCA2`.

The persisted diagnostic extraction is intentionally not rewritten in place:

`results/vb_full150_20260821/BRCA2/20260821_111529/extractions/BRCA2_PMID_26843898.json`

## External adversarial review

The installed Grok CLI was run read-only with model `grok-4.6`, reasoning
effort `xhigh` (the installed maximum), plan permission mode, subagents and web
search disabled. It was asked to assess the general policy, the motivating
paper, and the locked 119-paper evaluation without proposing PMID allowlists or
gold-dependent production behavior.

The completed review found:

- table-like exclusion alone is insufficient because background prose,
  references, testing panels, and other-gene observations can also carry valid
  HGVS-looking tokens;
- scanner hints and additive merge have different trust roles and must be
  separated;
- target-gene and current-study attribution must both be independently positive
  before a scanner row can be appended; silence must abstain;
- scanner confidence and the 150-candidate cap are not attribution controls;
- skip reasons must be persisted; and
- PMID 26843898, wrong-gene, background, bibliography, table-like, duplicate,
  and locked-cohort replays are required acceptance cases.

Grok measured 193 scanner identities on the collaborator-facing trusted
projection, with 20 exclusive gold matches and 173 extras under its projection
comparison. A later repository-local raw extraction replay counted 428 scanner
rows before the trusted projection; these are different surfaces and should not
be conflated.

A second Grok xhigh review was launched against the actual uncommitted diff. It
ran repository and real-paper probes but did not return a final verdict within
the bounded review window and was interrupted. It made no edits and did not
report a new concrete failing case before interruption.

## Implemented protocol boundary

- Regex candidates remain bounded, non-authoritative prompt hints even when the
  additive merge circuit breaker fires.
- Additive merge normalizes identities and recognizes underspecified indel
  prefixes as already extracted when an authoritative row has the complete
  identity.
- A new scanner row requires target-gene attribution and current-study
  observation language in the same independently supportable mention.
- Table-like, compilation/background, bibliography, negative/unconfirmed, and
  wrong-gene contexts abstain.
- Concatenated gene-variant tokens such as `KCNH2_R176W` are bound to their
  embedded gene before the requested gene can be stamped.
- Scanner rows never acquire patient, penetrance, or individual counts.
- Metadata records policy version, additions, and per-candidate skip reasons.
- Scanner offsets are retained when an exact source mention is available.

No PMID branch, gold identity keep-list, paper-specific alias, or public
publication path was added.

## Exact real-paper regression

The old scanner rows were removed from a copy of the authoritative extraction,
the current scanner was rerun on the actual cleaned full text, and the new merge
was applied with confidence 0.6.

| Signal | Result |
|---|---:|
| Authoritative input variants | 7 |
| Regex candidates | 57 |
| Final variants | 7 |
| Scanner additions | 0 |
| Table-like skips | 43 |
| Already-extracted skips | 11 |
| Wrong-gene skips | 2 |
| Background skips | 1 |

No `C61G`, PALB2 `c.509_510delGA`, truncated background
`c.7913_7917del`, or literature-table row remained.

## Locked 119-paper merge replay

The finalized merge policy was replayed over all 119 locked gene-paper attempts
using each frozen source and its prior non-scanner extraction as the baseline.
Gold was consulted only afterward by the evaluation matcher; it was not read by
the scanner or merge control flow.

| Signal | Result |
|---|---:|
| Gene-paper attempts | 119 |
| Old raw scanner rows | 428 |
| New scanner rows | 1 |
| Exclusive gold matches retained | 1 (`SCN5A`, PMID 19305639, `H558R`) |
| New extras | 0 |

Final skip reasons across the replay:

| Reason | Candidates |
|---|---:|
| already_extracted | 612 |
| table_like | 350 |
| wrong_gene | 186 |
| gene_unattributed | 9 |
| reference_list | 8 |
| study_unattributed | 6 |
| background_mention | 2 |

Several apparent recall losses are scientifically desirable. For example, the
five RYR2 calls in PMID 18285261 are present in the normalized gold file, but
the paper states that repeated isolation and analysis could not confirm them
and that they were likely FF-PET-derived DNA artifacts. Gold membership was not
allowed to override that source evidence.

## Tests and disposition

The targeted scanner and extraction-table suite passed 154 tests. The scanner
unit surface includes current-study acceptance, silence abstention, wrong-gene
and embedded-gene rejection, bibliography and unconfirmed-artifact rejection,
background/compilation rejection, multi-mention recovery, normalized identity
deduplication, the BRCA2 compilation regression, and the over-cap hints/merge
split.

This is a GO for a fresh locked-119 production run and fresh private
BRCA1/BRCA2/BMPR2 extraction. It is not evidence for public annotation
publication. Public publication remains OFF pending the repository acceptance
gates and collaborator/human review described in `TASKS.md`.
