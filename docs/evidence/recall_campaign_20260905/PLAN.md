# Recall implementation campaign, 2026-09-05

This is the new approximately $100 authorization, separate from the closed
$84.035 campaign. The goal is fewer paper-derived identity FN and more TP,
while preserving precision and source-supported count semantics. Work stays
on main; validated changes will be committed and pushed.

Before any new model extraction, the candidate contains source-only protein
substitution preservation, explicit cDNA-column pairing, and identity-only
gene-labelled case alleles excluded from figures/tables. Generic “Nucleotide
change” columns are deliberately excluded from the new pairing path: a header
alone does not establish coding coordinates. No phenotype/carrier values are
introduced by these rules. The existing downstream attribution/range, trust,
and study-ownership checks remain active.

The repository fallback is hardened against lost caption metadata, malformed
redirects and DOI render glue. Existing enabled browser recovery precedes the
body-only repository route so publisher supplements still get a chance.
Incomplete supplement status remains a retry signal: presence of one converted
supplement does not establish completeness.

The source-only prototype scan covered continuation tranches 01–03, without
reading new answer keys. The broader initial pairing prototype touched five
opened attempts in 01, zero in 02, and one in 03. The final explicit-cDNA rule
is narrower. This is trigger coverage, not a count-source completeness audit.
The full three-rule bundle has no established population gain yet.

1. Run a fresh, source-frozen, gold-free calibration on the ten attempts in
   `mechanism_panel.tsv` (eight distinct papers, including one multi-gene article),
   including the recovered 20031634, parsing misses, partial count capture,
   two exact count controls, and a functional negative control. Reserve $5.
   These papers are already opened and cannot support confirmation.
2. Trace any changes against the old locked artifacts and replay deterministic
   stages from source/intermediates. Separate recovery of new source bytes from
   parsing effects and from run-to-run model differences. No old prediction or
   database feeds the fresh extraction; old locks remain immutable.
3. Only spend on a full paired tranche after the mechanism evidence justifies
   it. Use matched frozen sources for reading contrasts. For an end-to-end
   acquisition contrast, explicitly name differing source hashes and use the
   same reader in a source-only ablation. Reserve enough for replication.
4. Preserve the existing registered discovery/confirmation thresholds. Check
   overlap with *all* opened locks before calling a tranche a holdout. A
   calibration improvement is not protocol acceptance, and sample bootstrap
   bounds do not measure provider jitter. No accepted headline changes without
   a genuine passing discovery and confirmation.

An explicit frozen-source mode verifies every selected text/asset hash and
prevents live harvesting or prior-run cache substitution. This closes a gap
in the old harness: disabling the later recovery stage did not disable the
initial harvester, especially for body-only retry sources. Both baseline and
candidate use this measurement mode. It does not change ordinary cache policy.

Reviewer dispositions: all three CLIs completed a review (Grok and Agy needed
bounded facts-only retries after their file reviews hit a turn/time limit).
Accepted: scope ownership, no ungrounded cDNA prefix, immutable sources,
mechanism checks before large spend, and replication for small effects.
Rejected: removing missing-source papers from the end-to-end denominator;
calling a different unopened tranche a repeat of the same papers; treating
figure/unsearchable notation as proven absent acquisition; equating a +2 TP
change with statistical significance. Claude's alleged cDNA-triggered VF bug
needs direct reproduction: the classifier also differs between one- and
three-letter protein notation. No VF quarantine relaxation is implemented.

Costs and subsequent decisions are recorded in `budget.json` and the report.

## Decisions after the first lock

The ten-attempt prototype finished at 348 TP / 41 FP / 12 FN. The prose
addition recovered two RYR2 gold identities but added five gold extras across
KCNQ1/SCN5A, so all three speculative parser additions were withdrawn. The
source-only substitution examples already survived other stages; the cDNA
pairing rule had no demonstrated incremental win.

The recovered 20031634 body yielded ten gold identities and exposed a real
validation loss: its structured response emitted `p.1570-Phe1571insIle` with
11 carriers, six affected and four unaffected, but validation discarded it.
The final code accepts this explicitly bounded adjacent-position shorthand
without inventing its omitted left residue. A zero-model replay preserves the
row through extraction and migration. A fresh four-attempt ablation uses the
same source/asset hashes as the prototype, including all three gene attempts
for the prose-rule article. Reserve $3.

Before that ablation finishes, add one final-code check of the other recovered
body, SCN5A 25163546 (previously 0 TP / 20 FN). This directly tests the utility
of the second repository recovery; it is already-opened calibration, not a
new holdout. Reserve $2. No further model run is justified without a specific
unresolved mechanism.
