# Case-report disappearance and variant-evidence resurfacing

Created: 2026-05-24
Source: Brett voice note, 2026-05-24.

## Why this belongs in GeneVariantFetcher

Part of the long-term point of GeneVariantFetcher is not just extracting variants from papers that already exist. It should help answer a broader publication-ecosystem question:

> Have clinically informative variant case reports disappeared, and can we build tools or incentives that make those observations visible again?

Brett's intuition is that older variant literature often contained case reports because the observation itself was interesting: a patient, a phenotype, and a variant that plausibly connected to disease. Over time, the field shifted toward large cohorts and center-level aggregation. Those cohort papers are valuable, but the case-report stream may have collapsed. After the big cohort era, many clinically meaningful one-off or small-N observations may now remain trapped locally in clinics, genetic-testing reports, EHRs, or internal case conferences instead of entering public literature.

If true, that is a problem for variant interpretation. Case-level observations can be disproportionately informative for:
- penetrance and expressivity;
- phenotype expansion;
- genotype-first discovery;
- rare variant interpretation;
- prior probability for related variants in the same gene/domain;
- clinical narratives that large cohorts flatten away.

## Analysis the repo should support

### Primary question

For a given gene, disease area, or variant class, how has the literature composition changed over time?

Specifically:
- early case reports / small family reports;
- later larger cohorts / registry papers / center series;
- possible decline of case reports after cohort consolidation;
- whether new individual-level case observations stop appearing even though clinical testing volume presumably increases.

### Hypothesis

A plausible pattern:

1. **Discovery / early clinical era** — case reports and family reports appear because each observation is novel.
2. **Aggregation era** — high-impact cohort papers and multicenter studies dominate.
3. **Post-cohort quiet era** — case reports become less publishable, even though new informative patients still exist.

This may be wrong. The point is to measure it, not assume it.

## Candidate metrics

For each paper discovered by GeneVariantFetcher, try to classify:

### Paper type
- single-patient case report;
- family report / pedigree;
- small case series;
- cohort / registry / multicenter study;
- functional-only paper;
- review / guideline / database paper;
- methods paper.

### Evidence granularity
- named variant(s) with individual-level phenotype;
- variant table with per-carrier details;
- aggregate variant counts only;
- no extractable patient-level variant evidence.

### Time trend
- publication year;
- year of first report for each variant;
- number of new variants reported per year;
- number of papers with individual-level variant-carrier data per year;
- ratio of case reports to cohort papers by year.

### Publication incentives / visibility
- journal;
- journal impact factor or proxy rank;
- article type label from PubMed / publisher when available;
- citation counts if easy to retrieve;
- open-access status;
- whether supplementary tables contain the actual variant-carrier payload.

### Variant-level outcomes
- variants reported only once vs repeatedly;
- variants first reported in case reports then later in cohorts;
- variants that appear in ClinVar / gnomAD / disease databases but have thin case literature;
- genes where case reports disappear while variant submissions continue elsewhere.

## Implementation hooks for GeneVariantFetcher

### 1. Add publication-type classification

For each fetched paper, classify article type using a staged approach:

1. PubMed publication types when available.
2. Title/abstract heuristics: `case report`, `case series`, `family`, `kindred`, `cohort`, `registry`, `multicenter`, `systematic review`.
3. LLM-assisted classification from title/abstract/full text.
4. Confidence score and rationale string.

Suggested fields:

```text
publication_type_primary
publication_type_secondary
publication_type_confidence
publication_type_rationale
is_case_report_like
is_cohort_like
has_individual_level_variant_data
```

### 2. Add per-paper evidence-yield metrics

For each paper:

```text
n_variants_extracted
n_carriers_extracted
n_families_or_kindreds
has_supplement_variant_table
has_patient_level_phenotype
has_functional_data
```

### 3. Add time-series report command

Potential CLI target:

```bash
gvf analyze publication-trends --gene KCNH2 --disease LQTS --out reports/kcnh2_publication_trends.html
```

Outputs:
- stacked bar chart of paper types by year;
- new variants per year;
- individual-level evidence papers per year;
- case-report/cohort ratio over time;
- top journals by period;
- candidate “evidence deserts” where testing likely continues but public case evidence dries up.

### 4. Link to impact-factor / journal-rank proxies

Impact factor itself may be annoying/licensed. Safer proxies:
- journal title;
- PubMed citation count / iCite if available;
- Scimago rank if accessible;
- OpenAlex source metadata;
- Crossref / Semantic Scholar citation counts.

The analysis should not depend on exact JIF. It can ask whether case-level observations migrated from higher-visibility journals to lower-visibility outlets or disappeared entirely.

## Candidate test genes / domains

Good first targets because Brett knows the biology and literature:
- `KCNH2`
- `KCNQ1`
- `SCN5A`
- broader congenital long-QT genes
- cardiomyopathy genes such as `TTN`, `MYH7`, `LMNA` if the pipeline generalizes

## Separate strategy stream: how to stimulate case evidence back into existence

### Core problem

Clinicians and health systems probably still see informative rare variant cases. But the publishing path is too slow, unrewarding, and low-status for many of those observations to become accessible evidence.

The goal is not to resurrect low-quality anecdote. The goal is to create a credible, structured, privacy-safe way for clinically informative observations to re-enter shared knowledge.

### Possible strategy 1 — Advocate/Wake as a scale partner

Pitch to Advocate/Wake:

> You touch enough patients that rare variant observations that look like isolated anecdotes locally become a meaningful evidence stream at system scale. We can turn those observations into structured, shareable variant-phenotype evidence.

Why this is attractive:
- large patient volume;
- many genetic tests and cardiology encounters;
- possible EHR/genomics integration;
- institutional interest in learning-health-system infrastructure;
- research value from converting local clinical insight into public evidence.

### Possible strategy 2 — Structured variant case briefs

Instead of traditional case reports, create a lightweight standardized format:

- gene / variant / transcript;
- phenotype ontology terms;
- age/sex/ancestry where appropriate and privacy-safe;
- testing indication;
- segregation / family history;
- relevant ECG/imaging/labs;
- treatment or outcome if relevant;
- ACMG evidence codes potentially affected;
- consent / IRB / de-identification status;
- links to ClinVar submission or publication.

This could feed:
- ClinVar submissions;
- journal short reports;
- preprint-style case evidence notes;
- internal variant boards;
- public gene-specific evidence dashboards.

### Possible strategy 3 — “Evidence yield” dashboard

Use GeneVariantFetcher to show that case evidence has declined, then offer a dashboard that tracks new internal observations waiting to be converted into shareable evidence.

Dashboard buckets:
- likely publishable case report;
- ClinVar-ready evidence update;
- needs segregation;
- needs phenotyping cleanup;
- already known / low novelty;
- ethically/privacy constrained.

### Possible strategy 4 — Make it academically motivating

People will not do this if it feels like clerical work. Possible incentives:
- coauthored annual gene/disease evidence updates;
- credit for clinicians who contribute structured cases;
- integration with fellows/residents as scholarly output;
- automated first drafts for case briefs;
- institutional “variant evidence commons” branding;
- measurable impact on variant reclassification.

### Possible strategy 5 — Avoid the garbage-case-report trap

The standard must be higher than “we saw a variant and a disease.” Required guardrails:
- clear phenotype quality threshold;
- population frequency check;
- inheritance/segregation where possible;
- prior literature and ClinVar context;
- uncertainty language;
- no overclaiming causality;
- privacy protection.

## Talk/pitch framing

A blunt but useful line:

> The field solved part of the aggregation problem, but may have accidentally killed the public case-observation stream that made variant interpretation possible in the first place.

Better polished version:

> As genomic testing scales, many informative variant observations are now generated inside health systems but never become part of the public evidence base. We need infrastructure that turns local clinical observations into calibrated, shareable variant evidence.

## Next concrete repo task

Add an issue or implementation note for:

1. paper-type classifier;
2. individual-level evidence-yield fields;
3. publication trend report;
4. first test run on KCNH2 / long-QT literature;
5. strategy memo for Advocate/Wake variant-evidence commons.
