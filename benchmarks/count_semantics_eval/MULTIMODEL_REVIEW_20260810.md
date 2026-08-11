# Independent multi-model review — 2026-08-10

## Reviewers

| Route | Model | Effort | Role |
|---|---|---|---|
| Grok CLI | Grok 4.5 | high (strongest available; `max` was unavailable) | implementation and source-scope review |
| AGY | Gemini 3.1 Pro High | high | implementation and source-scope review |
| Claude CLI | Claude Fable 5 | max | trust-chain and source-scope review |

The first pass reviewed the proposed failure-routing implementation and its
metrics. Grok passed it but identified a second scorer that did not honor v2
rows. AGY passed it subject to handling prefixless three-letter protein aliases.
Claude blocked on the open-vocabulary status field, duplicated gold-loading
logic, incomplete count-field reporting, and the risk that adjudication had
followed predictions rather than sources.

The implementation response was to centralize gold loading, validate an exact
status vocabulary, fail closed on unknown states, honor explicit nulls and
exclusions in both scorers, scope the 64k xhigh allowance to GPT-5.6, add alias
coverage, report all three count fields, and conduct the blind audit below.

## Blind source-only audit

Fifteen evidence cards were randomized as B01–B15. The reviewers received the
prospective count policy and source excerpts, but no locked prediction, legacy
gold, v2 gold, error direction, or other reviewer's answer.

### Headline rows

All three reviewers independently selected the same carrier total:

| Card | Gene / PMID / variant | Grok | AGY | Claude | Decision |
|---|---|---:|---:|---:|---:|
| B01 | KCNQ1 / 33141630 / T224M | 124 | 124 | 124 | 124 |
| B02 | SCN5A / 20470418 / S1103Y | 26 | 26 | 26 | 26 |
| B06 | SCN5A / 26746457 / R1193Q | 19 | 19 | 19 | 19 |
| B08 | SCN5A / 28339995 / D1790G | 30 | 30 | 30 | 30 |
| B14 | KCNH2 / 19160088 / R176W | 16 | 16 | 16 | 16 |

### Control rows

| Card | Gene / PMID / variant | Three-model outcome | Action |
|---|---|---|---|
| B03 | SCN5A / 26921764 / R1193Q | 3 carriers | Retain 3; locked prediction is 2 |
| B10 | SCN5A / 26921764 / M369K | 3 carriers | Retain 3; locked prediction is 2 |
| B04 | RYR2 / 33606749 / T2390I | 2 carriers | Correct 1 → 2 |
| B12 | RYR2 / 33606749 / T1223A | 2 carriers | Correct 1 → 2 |
| B13 | RYR2 / 33606749 / R420Q | 2 carriers | Correct 1 → 2 |
| B07 | KCNH2 / 10841244 / L552S | 42 people | Correct allele copies 44 → 42 people |
| B09 | RYR2 / 25814417 / G357S | 185 carriers | Correct living-only 179 → 185 genotyped carriers |
| B15 | SCN5A / 18451998 / E1784K | 41 carriers | Correct 50 → 41 |
| B05 | RYR2 / 28237968 / c.13352del | Insufficient card evidence | No change |
| B11 | KCNH2 / 10862094 / c.526C>T | Insufficient card evidence | No change |

The two retained gold=3 controls move against the locked prediction, while the
R420Q affected correction increases affected-count error by one. These outcomes
show that the audit was not simply optimizing the answer key toward predictions.

## Phenotype-partition handling

Reviewer agreement was strongest on carrier scope. Some phenotype decisions
varied because the papers use different axes (symptoms, QTc, disease status,
case/control enrollment) or assess only a subset. Final partitions use only the
paper's explicit definition:

- T224M remains 34/54 among 88 followed carriers, with 36 additional carriers
  unpartitioned.
- S1103Y remains the paper-defined 17/9 case/control cohort.
- R1193Q uses the explicit 7/19 phenotype result.
- D1790G uses the 30 enrolled LQT3 patients; carrier scope is independently
  robust even though Claude preferred null phenotype partitions.
- R176W is 16/null/null because the population survey did not systematically
  ascertain symptoms.
- L552S and G357S retain explicit subset partitions that do not sum to all
  carrier people.

## Verdict

After the fixes and blind audit, the review outcome is mergeable. Remaining
source-card omissions are documented as follow-up work and do not affect the
adjudicated rows or the locked 56-paper result.
