# PMID 29650123 publisher-source audit — 2026-08-24

## Finding

The current public publisher assets do not expose the 22-row KCNH2 mutation
roster represented in the cardiac gold fixture. The pipeline's 2/22 result is
therefore source-limited; the 20 missing identities must not be reconstructed
from the gold file or inferred as one carrier each.

## Checks

- PubMed identity: PMID 29650123, DOI `10.1016/j.jacc.2018.01.078`.
- The live JACC article DOM exposes one supplementary download:
  `suppl_file/mmc1.docx`.
- That exact DOCX is already stored as
  `corpus/KCNH2/29650123/29650123_supplements/mmc1.docx` and folded into the
  full context.
- Direct DOCX inspection finds Online Table 1 (study-population characteristics
  by diagnosis decade), Online Table 2 (beta-blocker dosage), and Online
  Figures 1–3 (clinical/outcome plots). It contains no variant roster.
- The article page exposed no second supplementary mutation-table asset.

## Consequence

This paper accounts for 20/86 (23.3%) of the remaining source-recovery-lock
false negatives. Recovering a genuine roster would have a theoretical ceiling
of +20/632 = +3.16 recall points, but no improvement is claimable until that
source is obtained and bound. The correct status is
`source/provenance blocked`, not `parser failed` and not `infer carriers=1`.

Public records:

- <https://pubmed.ncbi.nlm.nih.gov/29650123/>
- <https://www.jacc.org/doi/10.1016/j.jacc.2018.01.078>
- <https://www.jacc.org/doi/suppl/10.1016/j.jacc.2018.01.078/suppl_file/mmc1.docx>
