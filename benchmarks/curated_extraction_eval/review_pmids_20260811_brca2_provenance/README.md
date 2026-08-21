# BRCA2 active review cohort after provenance correction

This is the exact Variant Browser BRCA2 review cohort adopted on 2026-08-11.
It starts from the 50-paper operational cohort in `../review_pmids_50/BRCA2.txt`
and excludes papers that fail either collaborator provenance or human clinical
scope:

- 10398279 and 21356067 were already absent from the operational 50-paper
  review cohort.
- 15365993, 18489799, 22655046, and 25802882 were removed from the active
  Variant Browser snapshot.
- 19944633 is an explicitly canine BRCA2 ortholog paper. It remains only in the
  wrong-species negative controls and must never enter the clinical browser.

The two papers reviewed by Nate and lead-approved in Variant Browser—26833046
and 26848529—remain present. The resulting operational review cohort has 45
papers. This file is a destructive-import allowlist: publishing it replaces the
current BRCA2 paper/evidence snapshot with exactly these PMIDs.
