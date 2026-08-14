# Fixed 146-attempt experimental run: cost and QC

This run keeps the experimental scope fixed at BMPR2 50, BRCA1 50, and BRCA2
46. The separately completed 120-attempt cardiac run is the requested
100--150-paper comparison against manually curated gold; it is not an
experimental-cohort expansion.

## LLM burn

The three genes ran concurrently. Provider time is the sum of traced call
latencies, not wall-clock time. Cost is a public-list-price proxy using $5/M
input and $30/M output for GPT-5.6 Sol, $1.25/M and $2.50/M for Grok 4.3, and
$0.95/M and $4/M for Kimi K2.6. Cached-input discounts and human time are not
included.

| Gene | Attempts | Observed wall | Calls (successful) | Tokens | Summed provider time | Cost proxy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BMPR2 | 50 | 70.7 min | 195 (195) | 1,023,662 | 43.8 min | $4.275 |
| BRCA1 | 50 | 91.6 min | 407 (405) | 1,559,111 | 86.3 min | $9.112 |
| BRCA2 | 46 | 97.8 min | 370 (370) | 1,678,568 | 92.6 min | $10.278 |
| **Total** | **146** | **97.8 min concurrent wall** | **972 (970)** | **4,261,341** | **3.712 h** | **$23.664** |

BMPR2 wall time includes the original workflow plus the targeted two-paper
resume after the extraction-marker circuit-breaker repair. Its final
`RUN_STATUS.json` duration covers only the 518-second resume, so that field must
not be presented as the whole-gene elapsed time.

| Model role | Calls (successful) | Input tokens | Output/reasoning tokens | Provider time | Cost proxy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kimi K2.6 table routing | 68 (68) | 97,956 | 245,540 | 1,176.5 s | $1.075 |
| Grok 4.3 primary extraction | 118 (118) | 1,415,249 | 651,355 | 6,069.2 s | $3.397 |
| GPT-5.6 Sol verification and figure vision | 786 (784) | 1,453,816 | 397,425 | 6,116.4 s | $19.192 |

The later pre-publish BRCA2 count-recovery audit added two GPT-5.6 Sol `high`
calls: 32,226 input and 910 output/reasoning tokens (33,136 total), 10.79 summed
provider-seconds, and a $0.188 proxy at the same $5/M + $30/M convention. Both
calls were write-time verified. They examined 26 count gaps across PMIDs
12942367 and 22382806, grounded zero, and wrote nothing. Including this bounded
audit, total LLM burn is 974 calls, 4,294,477 tokens, about 3.715 summed
provider-hours, and a $23.853 proxy. The fixed extraction itself remains the
972-call / $23.664 result above.

## Gold-free QC and trusted projection

| Gene | Full text | Variant identities | Paper-variant links | Count rows | Trusted | Quarantine |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BMPR2 | 45/50 | 587 | 483 | 252 | 231 | 21 |
| BRCA1 | 50/50 | 6,461 | 7,758 | 2,527 | 2,062 | 465 |
| BRCA2 | 45/46 source-integrity pass | 2,155 | 2,526 | 479 | 416 | 63 |

- The trust gate is `tg5-00fb2bd17402`. Family-count observations remain in
  the raw audit record but their carrier fields are masked. This quarantined
  390 BRCA1 and 46 BRCA2 family-count fields instead of importing families as
  patients.
- BRCA2 source QC flags PMID 14559878 as abstract-only and PMIDs 12942367 and
  22382806 as usable-source, zero-variant refresh targets.
- Report-only somatic/germline QC made no mutations. BRCA1 has 188 somatic and
  747 ambiguous records among 7,226 (12.94% combined); BRCA2 has 0 somatic and
  260 ambiguous among 1,420 (18.31%). These require paper-level review rather
  than an automatic drop.
- Every trace record in all three final manifests is write-time verified. The
  BMPR2 manifest truthfully records the two contributing run IDs.
- The extraction jobs themselves kept publication disabled. After the accepted
  gold-120 gate, the lead explicitly authorized refreshing the fixed
  collaborator-facing queues; that staging-only import is recorded below.

## Approved BRCA2 gold diagnostic

Only PMIDs 26833046 and 26848529 were scored, because only those papers have
lead-approved collaborator gold. The run recovered all 7/7 curated variant
identities. It supplied carrier counts for 3/7 assertions (42.86% coverage),
with carrier MAE 1.333. Raw variant precision and counted-extra precision are
not interpretable as paper-level precision here: the collaborator gold contains
seven selected records and is not an exhaustive inventory of every variant in
these two papers.

## Variant Browser collaborator refresh

On 2026-08-13 the three fixed outputs replaced carrier evidence in the existing
Variant Browser staging snapshots using the permanent 50/50/46 manifests. The
preflight found zero paper additions/removals and exact review order. Live
trusted evidence/individual/fact counts are BMPR2 482/470/3,871, BRCA1
7,260/3,663/27,172, and BRCA2 2,346/591/4,920.

All 111 pre-existing BRCA2 adjudication and gold stable keys were preserved.
Every prior adjudication now requires re-review and is detached from current
evidence; the default adjudication and gold exports both return zero. A guarded
identity collision was resolved by retaining the sole reviewed
`gvf:c.2808_2811delACAA` key and moving five state-free duplicate aliases into
a reversible archive namespace. Public annotations were not changed. Exact
input/export hashes and live source-run IDs are recorded in
`VARIANT_BROWSER_PUBLICATION.md`.

## Reviewer-tier estimate

The exact reviewer tier is 546 gene-paper attempts across 507 unique PMIDs.
Cost scales with attempts because a PMID appearing in two gene workspaces is
processed twice.

- Cardiac gold calibration: $0.0814 and 65.9 provider-seconds per attempt ->
  **$44.47 and 9.99 provider-hours** for 546.
- BRCA-heavy experimental calibration: $0.1621 and 91.5 provider-seconds per
  attempt -> **$88.50 and 13.88 provider-hours** for 546.
- The actual 11-gene roster already has seven measured/scaled strata totaling
  about $39.95. Applying the two calibration rates only to the four unmeasured
  50-paper genes gives a roster estimate of **$56--$72 and 11.0--12.5 summed
  provider-hours**. A practical LLM budget with retry headroom is **$70--$90**.

With three or four gene jobs in parallel, pure inference should fit in roughly
4--6 elapsed hours. Source acquisition, refreshes, integrity checks, and QC make
**10--20 elapsed hours** a safer operating window. Human curation time is
separate. For a generic 500--600 attempt population, the two measured
calibrations span **$40.72--$97.25** and **9.15--15.25 summed provider-hours**;
budget roughly **$50--$120** before human review.
