# BRCA2 sex, cancer endpoint and high-penetrance audit — 2026-09-14

**The saved BRCA2 model is not women-only.** It combines clinical observations
without a sex partition with all-sex gnomAD carriers. A female breast-cancer
analysis should use female clinical counts and exact gnomAD XX carrier counts,
keeping the user's assumption that eligible gnomAD carriers are unaffected.
Changing only the population denominator would leave male cases in the numerator.

The accepted five-gene starting point is committed on local `main` as
`dcad8e7a`, after the preceding correction `c8e34c96`. A fresh fetch showed
these two commits ahead of `origin/main`, with none behind. Public push was
rejected by automatic approval review: it requires explicit authorization to
publish the source-level variant/carrier evidence to public
`kroncke-lab/GeneVariantFetcher`. No push workaround was attempted. This audit
is additional evidence; it does not overwrite that frozen starting point.

## Are any BRCA2 variants genuinely 100% penetrant?

The reviewed clinical sources do not establish a heterozygous BRCA2 variant
that inevitably causes breast or ovarian cancer. In a large prospective study
of female carriers, BRCA2 breast-cancer risk by age 80 was **69% (95% CI 61–77%)**
and ovarian-cancer risk was **17% (11–25%)**. These are separate endpoints, not
percentages to add. The BRCA2 breast analysis observed 157 incident cancers in
1,610 women over 7,913.1 person-years; 157/1,610 is not its lifetime penetrance.
[Kuchenbaecker et al., JAMA 2017, PMID 28632866](https://jamanetwork.com/journals/jama/fullarticle/2632503).

Risk can be higher in selected families: one clinical genetics study estimated
BRCA2 breast-cancer penetrance at 88% by age 80. This is a family-ascertained
gene-level estimate, not evidence that a particular allele has 100% penetrance.
[Evans et al., PMID 18513387](https://pubmed.ncbi.nlm.nih.gov/18513387/).

The current structural plot is **missense only** and describes neighbors after
excluding the target variant. It excludes frameshift/nonsense variants, and
averages different substitutions at a residue. A low or moderate neighborhood
therefore does not establish that every allele at that residue is low risk.
Pathogenicity classification, functional loss, an observed case fraction and
age-specific cancer penetrance are different quantities.

In the frozen 4,444-unit BRCA2 missense input, **139 units have observed A>0,
U=0**, but **119 are single affected observations**. A proper finite Beta prior
keeps the posterior below 1 even when the observed fraction is 1. The sole
own posterior above 90%, E1581D at 92.1%, fails the source check below and must
not be cited as a high-penetrance BRCA2 allele.

## What the local source audit establishes

[The executable audit](audit_local.py) replays both previous observation-level
corrections and reproduces every corrected clinical key total. The canonical
missense union has A=1,858, literature U=1,482 and gnomAD U=1,516,774. These
remain observation totals, not globally deduplicated people.

There are 976 retained missense source observations from 227 papers. Only 196
observations have linked individual records, and those records have not been
reconciled to the aggregate counts. Missing sex metadata cannot be taken as
female. A title screen, including cached abstract titles, flags 37 observations
with 48 assigned affected counts for male/prostate review. That screen is a
queue, not a complete or adjudicated sex partition.

**Confirmed male example:** D2723H currently contributes two affected carriers
from PMID 20927582. The primary report explicitly identifies two unrelated men
with this allele in a male breast-cancer study. Those cases belong in a male
endpoint and cannot contribute to a female numerator.
[Ding et al., primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3059396/).

**New source errors:** the [reviewed source patch](reviewed_source_patch.csv)
records exact affected observations, source paths and line locations:

- E1581D / PMID 17100994 / observation 40000141: remove 48 BRCA2 affected
  assignments. The first table contains the BRCA1 mutations identified in the
  Results text; E1581D belongs to that table. Its row has 48 patient and 3
  normal-control observations, not BRCA2 evidence. A coincidentally compatible
  amino acid at the BRCA2 position does not establish gene identity.
- I3412V / the same paper / observation 40000163: preserve 26 BRCA2 patient
  observations and recover 5 omitted normal-control observations from the
  second table. The final value, 110, counts BIC entries and is not people.
- D2723H / PMID 20927582 / observation 20001234: remove its two affected men
  only from a proposed female endpoint; retain their male breast-cancer evidence.

The first paper is [Han et al., PMID 17100994](https://pubmed.ncbi.nlm.nih.gov/17100994/).
Its locally cached primary text and the male study are hashed in
[summary_input_hashes.json](summary_input_hashes.json). The archived databases
and frozen plots have not been rewritten. These decisions must be applied
through the reviewed source workflow before the next model run; the current
BRCA2 plot remains provisional.

## Exact gnomAD sex counts are available

Six public GRCh38 alleles were queried through the single-variant `gnomad_r4`
API, using joint counts compatible with the existing v4.1.1 population inputs.
[Fetcher](fetch_sex_probe.py), [dated receipts](sex_probe.json) and
[arithmetic checks](build_summary.py) preserve the query, source hashes and
deduplication rule. All six combined-sex carrier counts reproduce the frozen
inputs exactly. This is a six-allele probe, not a full-gene sex recount.

| Variant | All-sex carriers | XX carriers | XY carriers |
|---|---:|---:|---:|
| K2729N | 709 | 347 | 362 |
| D2723H | 21 | 14 | 7 |
| R3052W | 8 | 3 | 5 |
| W2626C | 1 | 0 | 1 |
| N372H | 382,581 | 191,239 | 191,342 |
| I2675V | 1 | 0 | 1 |

For this autosomal gene, carrier observations are **AC minus homozygote count**
within the chosen stratum. AN is the number of sampled chromosomes, not the
number of carriers. Use exact XX counts, not half of AC, AN or the combined
carrier total. XX is inferred chromosomal sex, a practical proxy for the
clinical female stratum; it does not supply cancer status, age or censoring.
See the [gnomAD v4.1 release documentation](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/)
and [sex-label changelog](https://gnomad.broadinstitute.org/news/changelog/).

Each tested joint API response repeated the same overall XX and XY rows twice.
The fetcher removes only identical repetitions by exact stratum ID, rejects
conflicting repeats, and verifies XX+XY AC/AN/homozygotes against the joint
total. Ancestry-by-sex and HGDP/TGP subgroup rows are not added again.

[Denominator-only sensitivity](sex_denominator_sensitivity.csv) holds the old
prior and both clinical counts fixed. R3052W's own posterior changes from
40.2% to 53.4%; D2723H from 19.1% to 24.4%. **These are not female penetrance
estimates**: D2723H still has its two male cases in this diagnostic numerator.
They quantify the denominator effect only. Sex restriction can also remove
clinical cases, so a properly matched refit need not increase every estimate.

## Recommended next run

1. Make **female breast cancer** the first BRCA2 endpoint. Curate female
   germline affected and unaffected observations; keep ovarian cancer and male
   breast/prostate outcomes separate. An eventual breast-or-ovarian composite
   needs an explicitly deduplicated union. Unknown sex/phenotype stays unknown.
2. Apply the reviewed gene and control-count corrections, then resolve the
   remaining germline/somatic and endpoint queue. Freeze person/family/cohort
   ownership, age at assessment, onset and risk-reducing surgery where supplied.
3. Recover exact joint XX counts for all eligible alleles, preserving release,
   QC, AC/AN/homozygote and allele provenance. Under the requested convention,
   all retained XX gnomAD carriers add to beta as unaffected. A population-only
   allele observed solely in XY has no observed female donor; do not give it a
   fabricated female unaffected singleton. Preserve independently observed
   female literature evidence for that allele.
4. Rebuild the observed population/clinical union and fit separate missense
   and nonsense priors. Rerun the unchanged short positive-tail distance kernel,
   local structural frames, same-IDR polymer and variant-only LOO. Keep other
   substitutions at the target residue; compare raw fractions, posterior donors
   and one-prior pooled-count features with support displayed.
5. Check known pathogenic missense examples and truncating variants in their
   appropriate analyses, then assess age-specific calibration on independently
   ascertained outcomes. A simple literature-versus-gnomAD count ratio cannot
   reproduce a survival estimate just by matching sex.

Men are not universally unaffected by BRCA2: male breast-cancer risk is
1.8–7.1% by age 70 for harmful variants, and prostate and pancreatic cancer
risks are also elevated. Conversely, a benign variant does not imply zero
absolute lifetime breast-cancer risk in a woman; population risk is about 13%.
Near-zero scores under the adopted gnomAD-unaffected convention therefore
must not be labeled absolute lifetime risk.
[NCI BRCA fact sheet](https://www.cancer.gov/about-cancer/causes-prevention/genetics/brca-fact-sheet).

## Accepted prior-design starting point

The prior concern remains valid and is independent of the sex correction.
[Reproduced diagnostics](prior_weight_diagnostic.json) show the historical
weight w=1−1/(n+0.01) gives n=1 weight 0.0099 but n=2 weight 0.5025. In GCK,
53.3% of variants are singletons yet receive only 1.57% of fitting weight.
The variance divides the weighted squared deviations by variant count M,
rather than sum of weights. This strengthens shrinkage: GCK's unaffected
singleton posterior is 25.78%. Merely normalizing the weighted MSE reduces it
to 8.61%, but is only a diagnostic because observed-fraction variance also
contains binomial sampling noise.

Next compare a proper gene-by-type beta-binomial marginal-likelihood fit with
the historical method; assess whether a near-zero component plus a higher-risk
component is supported, and whether structural/functional features predict
component membership. Do not lower a prior merely to force a preferred curve.
Use held-out counts/outcomes rather than empirical posteriors as ground truth;
refit priors and feature construction inside validation folds to prevent leakage.
The descriptive full-dataset-prior LOO baseline remains saved. **No new prior
or mixture has been fitted in this audit.**

## Reproduction and validation

Run `audit_local.py`, then `build_summary.py`, using the existing BPE Python
environment with pandas/numpy. `fetch_sex_probe.py` replays cached receipts or
queries the public API when raw files are absent. Local raw API receipts and
archived source databases stay in ignored `results/`; this folder contains
compact evidence only, without individual identifiers or full article text.

Checks cover eight immutable database snapshots, both earlier correction
ledgers, all clinical key totals, the 4,444-unit missense union, six population
count reconciliations, five prior diagnostics and three source decisions.
The 111 repeated original (protein key, PMID, variant ID) keys are preserved at
their original row grain, not blindly deduplicated; cohort ownership remains
unresolved. These are audit checks, not a new extraction benchmark or calibrated
women-only model. [Local checks](local_checks.json) and
[input hashes](input_hashes.json) make the reconstruction inspectable.
