# Twelve-paper review cohort (superseded — frozen)

> **Superseded.** `../review_pmids_50/` was the next historical cohort; the
> current rollout manifests live only in
> [`../../evaluation_tiers/`](../../evaluation_tiers/README.md). This directory
> remains frozen for the July 9, 2026 experiment and is not active.

These PMID lists freeze the 12-paper-per-gene cohort used by the July 9, 2026
review experiment. They are a small operational review/smoke sample, not a gold
standard and not a replacement for the authoritative curated benchmark defined
by `../registry.tsv` and `../pmids/`.

The four cardiac lists are the first 12 entries from their corresponding
curated benchmark lists. The non-cardiac lists include the curated papers plus
additional cold-start papers used in that review run. Keep these files stable
for reproducibility; put a materially different future cohort in a new
versioned directory.
