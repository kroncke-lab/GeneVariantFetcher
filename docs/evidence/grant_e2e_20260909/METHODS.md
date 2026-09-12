# Methods (draft for the grant application) — frozen workflow, 2026-09-09

Wording is descriptive of what the frozen code does; numbers are filled in from
the run artifacts once each gene completes (see `analysis/SUMMARY.md`).

## 1. Literature collection

For each gene–disease pair (HNF1A–diabetes, GCK–diabetes/GCK-MODY,
LDLR–hypercholesterolaemia, BRCA2–breast cancer) candidate papers are discovered
from PubMed E-utilities with five gene-anchored keyword queries (variant,
pathogenicity, patient/cohort and a carrier/penetrance/segregation lane), each
constrained to an exact disease phrase in title or abstract, unioned with the
gene's PubMind literature set. Abstracts are fetched and passed through a
two-tier relevance filter: a deterministic scope gate (non-human orthologue
rejection) and a keyword tier, then an LLM tier (Azure `gpt-5.6-luna`, maximum
reasoning) that keeps papers reporting original genotyped human data for the
target condition and fails open when the abstract explicitly names the gene
together with genotyped-cohort signals. Filter decisions and reasons are logged
per paper.

## 2. Full-text acquisition

Full text is harvested in priority order from PubMed Central (BioC XML),
publisher APIs (Elsevier with institutional token, Springer, Wiley), indexed
article copies (Unpaywall, OpenAlex, HAL), and advertised supplementary files;
previously acquired source is reused from the consolidated corpus cache. After
extraction a gold-free source-QC pass classifies every paper (usable full text,
fetchable, supplement-only, blocked) and a recovery loop re-fetches and replays
what it can. Abstract-only papers are extracted from the abstract and flagged.

## 3. Variant and carrier-count extraction

Each paper's text and tables are condensed to high-value data zones, clinical
tables are routed by a cheap model (`Kimi-K2.6`) to a deterministic table parser
that reads variant identities and per-variant carrier / affected / unaffected
counts directly from table rows, and the primary reader (`grok-4.3`) extracts
variants, cohort counts and per-individual records from the remaining text.
Count-bearing claims are verified against the source by a second model
(`gpt-5.6-sol`) and disagreements adjudicated; every count carries fact-level
provenance (table, row, column, quote). A default-on trust gate assigns each
count row a tier (trusted / quarantine) from gold-free rules (impossible
partitions, within-paper outliers without structural provenance, cross-gene
attribution, somatic/germline ambiguity); only trusted rows enter the analysis.
Results are aggregated into a normalized SQLite database per gene.

## 4. From extraction to carrier observations

Observations are variant × paper rows. Within a paper, cohort-level rows take
precedence over per-individual records (to avoid double counting). A row is
usable when affected and unaffected counts are explicit (`explicit`), or when
affected and total carriers are explicit with no uncertain individuals
(`total_derived`). Variants are keyed at the protein level (one-letter
notation; nonsense as `X`, frameshift as `fs`), or by normalized cDNA when no
protein change is given. Counts are pooled across papers per variant; duplicate
notations of one variant within one paper are collapsed to the largest
denominator.

## 5. Variant features and classification

Variant-level features come from the lab's variantFeatures warehouse joined on
the MANE Select transcript's enumerated consequences: AlphaMissense, REVEL,
CADD, BayesDel (mean over alleles producing the same protein change), AlphaFold
pLDDT at the residue, gnomAD allele frequency, and the ClinVar classification
(the allele matching the paper's cDNA when available, otherwise the best
review-status record; allele-level conflicts are flagged). ClinVar labels are
collapsed to P/LP, VUS, Conflicting, B/LB, Other, or no record.

## 6. Bayesian penetrance estimation

The unit of observation is a variant × paper row *r* with *y_r* affected among
*n_r* carriers; single-carrier reports of a variant are merged into one row per
variant and ascertainment class. Because most literature rows are series of
affected carriers, ascertainment is modelled explicitly: each row is either an
affected-ascertained series or an informative observation of penetrance,

  *y_r* ~ π_r · Binomial(*n_r*, *p*_asc) + (1 − π_r) · Binomial(*n_r*, *p_r*),
  logit(*p_r*) = α + **x**_v(r)·**β** + σ·*z*_v(r) + τ·*e_r*,
  logit(π_r) = g0 + g1 · [genotype-first ascertainment recorded for the paper],

with *z_v* ~ N(0, 1) (between-variant), *e_r* ~ N(0, 1) (between-study within
a variant), α ~ N(logit of the pooled rate, 1.5), β ~ N(0, 1) on standardized
columns, σ ~ HalfNormal(1), τ ~ HalfNormal(0.5), g0 ~ N(0, 1.5), g1 ~ N(0, 1)
and *p*_asc ~ Beta(19, 1). Population-reference rows (gnomAD allele counts as
unaffected carriers in the gnomAD-anchored arm) are informative by construction
and carry no study effect. The mixing weight learns how often papers of each
ascertainment type are pure affected series; the between-study effect absorbs
the remaining heterogeneity between cohorts reporting the same variant. The
variant-level penetrance *p_v* = sigmoid(α + **x**_v·**β** + σ·*z_v*) is the
estimand: the probability of disease among carriers not selected for being
affected. Its posterior is the refined penetrance, driven by the variant's
own informative counts when they exist and shrunk toward the feature-implied
prior when the variant is known only from case reports. Fixed-effect sets:
intercept only; truncating class; class + AlphaMissense (primary); class +
ClinVar indicators (a classification-only prior); and both. Models are fitted
by NUTS (PyMC, four chains, 1,000 warm-up and 1,000 sampling draws, target
acceptance 0.95) and accepted when R̂ ≤ 1.01 with no divergences.

## 7. Evaluation

*Held-out variants.* Prior models are compared by five-fold cross-validation by
variant: each fold's model predicts the rows of variants it never saw,
integrating both random effects, scored by per-carrier negative log score and
Brier score, log predictive density, calibration slope, and paired bootstrap
differences against the intercept-only and classification-only priors.

*Held-out observations of known variants (the clinical question).* For each
variant, its observations are split into a set A that informs the posterior
and a set B that is scored. Where the extraction records ascertainment per
paper, B is the variant's genotype-first observations (population screening or
cascade testing) and A its phenotype-first observations plus any gnomAD anchor,
so the test asks whether case-based literature plus features predicts the
penetrance seen when carriers are found by genotype. Where ascertainment is not
recorded (the KCNQ1 control), papers are split alternately by PMID. The variant
posterior given A is computed exactly on a grid (the mixture likelihood per row, study effects integrated by
Gauss–Hermite quadrature, with τ, π and *p*_asc at their posterior means) and compared, on the
same held-out carriers, with the features-only prior, the pooled rate, and the
observed A-rate of the variant's ClinVar class among the other variants; AUC
for "observed penetrance in B ≥ 50%" is reported for the same predictors and
for the ClinVar class itself.

*Classification performance.* AUC for identifying variants whose observed
penetrance is ≥ 50% (variants with ≥ 5 literature carriers), for ClinVar
class alone, AlphaMissense alone and the held-out Bayesian prior, plus
sensitivity, specificity, PPV and NPV for the rule "P/LP ⇒ high penetrance".
The in-sample posterior is not used as a classifier because it contains the
counts that define the label.

*Discordant variants.* Variants with ≥ 5 carriers whose 95% credible interval
lies entirely on the opposite side of 50% from their ClinVar class, listed
with counts and source PMIDs. Synonymous alleles are excluded from all fits.

*Common variants.* Alleles with gnomAD frequency above 0.1% are kept in the
model and treated as a known-outcome check: they carry at most a GWAS-scale
risk, so their estimated penetrance must be near the population risk even
though case-series literature reports them mostly among affected carriers.
Each report tabulates literature fraction versus posterior for these alleles.

*Arms.* Primary (literature carriers, ascertainment mixture), sensitivity
(single-carrier reports removed), gnomAD-anchored with the mixture, gnomAD-
anchored without the mixture (the classic literature-plus-population model with
a between-study effect), and a KCNQ1 positive control compared with the curated
Variant Browser values, which selects the anchored model to report.

## 8. Reproducibility

GeneVariantFetcher is frozen at tag `grant-freeze-20260909`; the analysis lives
in `BayesianPenetranceEstimator/iterations/grant_e2e_20260909/` with fixed
seeds; every stage writes its inputs' hashes and the run's parameters
(`FREEZE.md`, per-gene `fit_summary.json`, `observation_summary.json`,
`feature_join_summary.json`).
