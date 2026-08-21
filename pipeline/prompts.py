"""LLM Prompt Templates for Variant Extraction

This module centralizes all LLM prompts used in the extraction pipeline.

**Keep these tight.** Instruction overhead is paid on every paper and competes with
the paper itself for the context window: the pre-2026-07-26 pair spent ~6.8K tokens
of boilerplate per call, which on a median extraction exceeded the source text
actually sent. The rules below are terse and stated ONCE, in ``_CORE_RULES``, shared
by both prompts -- the previous per-prompt copies drifted apart and had to be
re-synced by hand (that is why ``TABLE_ATTRIBUTION_GUIDANCE`` exists as a separate
constant at all).

The authoritative statement of what extraction must capture and must refuse is
``docs/EXTRACTION_CONTRACT.md``. Edit that first, then mirror it here. The schema
fragments are load-bearing -- downstream parsers, the count classifier, and the
evidence-card validator read those exact keys, and
``tests/unit/test_prompts_count_provenance.py`` pins them.
"""

# Threshold for switching to compact extraction mode
HIGH_VARIANT_THRESHOLD = 30

# Canonical table-scope/attribution guidance. Lives here so every extractor shares
# one copy: it was previously carried only inside an ad-hoc eval harness and was
# lost when that harness was not retained. Consumed by
# benchmarks/codex_paper_eval/run_eval.py; pinned by
# tests/unit/test_codex_paper_eval.py. Keep it free of {} so it can be concatenated
# into str.format templates.
TABLE_ATTRIBUTION_GUIDANCE = """Understand each table before extracting from it:
- Read its caption, column headers, footnotes, legends, symbol definitions, and the
  prose introducing it.
- Decide its scope: this study's own observations, a compilation citing other
  studies, or both.
- Decide what each column counts (families, individuals, alleles, probands, cases)
  before mapping to carriers/affected/unaffected.
- Judge every row and count for attributability using whatever provenance signal the
  table offers - reference column, novel-vs-reported label, footnote marker, scoping
  caption. Assume no fixed format or keyword.
- Count only what this study observed; exclude anything credited to another
  publication. If provenance is ambiguous, say so and do not fabricate a count."""

# Every extraction rule, stated once. Brace-free so it concatenates into
# str.format templates.
_CORE_RULES = (
    """TABLES
"""
    + TABLE_ATTRIBUTION_GUIDANCE
    + """

CLINICAL vs FUNCTIONAL - decide before extracting any count
Functional/assay signals ("n cells", "replicates", "n=", current density,
trafficking, fluorescence, patch clamp, massively parallel assay, saturation
mutagenesis, systematic testing of every variant at a position) mean the sample
sizes are ASSAY REPLICATES, not patients: leave every count null, put metrics in
functional_data, note "Functional study - assay data only" in source_location, and
take clinical counts only from a column explicitly labelled patients, carriers,
cases, probands, or the disease. Sanity check: 50+ variants with similar counts and
~100% penetrance means you read an assay table as patients - re-read it.
Non-human subjects (cat/Maine Coon/Ragdoll, dog/canine, mouse, rat, zebrafish, cell
lines) are never human carriers; they belong in functional_data/additional_notes
unless the row explicitly describes human patients.

COUNTS
- Implicit counts are real: "a patient"/"a case"/"an individual"/"the proband" = 1
  carrier; "a healthy individual"/"an asymptomatic carrier" = 1 UNAFFECTED carrier;
  "an affected patient" = 1 AFFECTED carrier.
- affected = has the phenotype/symptoms. unaffected = a carrier the paper
  explicitly calls healthy, asymptomatic, clinically silent, or unaffected. NEVER
  infer unaffected and never default it to 0.
- Case report whose carrier's status is not described: count the carrier, leave the
  affected/unaffected split null - do not default to affected. "A novel mutation"
  with no described human patient is a variant, not a carrier.
- Same variant in several tables/cohorts: SUM independent cohorts, take the LARGER
  when they overlap or one is a subset. Fail toward completeness.
- Never push a study-wide, cohort, family-set, domain, or mutation-class total onto
  individual variants. The same applies to patients.phenotype, demographics, and
  penetrance_percentage ("35% of carriers overall", "median age 61"). Use a count
  only when the row, cell, or sentence is variant-specific; put aggregates in
  additional_notes.
- Never usable as a carrier count: "Total case(s)" when a "Carrier(s)" column
  exists, screening denominators, MAF / allele frequency / allele counts (gnomAD,
  ExAC, TOPMed, 1000 Genomes), "No. of occurrences", row identifiers (patient or
  case no., adult/child number), and bare het/hom/genotype columns. A frequency
  cell "5/200" means 5 carriers. Families are not individuals.

IDENTITY - required, and load-bearing
Emit a variant only when at least one of cdna_notation, protein_notation,
genomic_position, or variant_class + structural_description is real. Return zero
entries rather than a placeholder for: mutation classes ("34 patients with
pore-region mutations"), single nucleotide or amino-acid letters, gene symbols,
nan/NA/NaN, and clinical-significance strings in a notation field. Valid
protein_notation is (p.)?<AA><pos><change>; valid cdna_notation starts "c." plus a
digit; valid genomic_position has a chromosome and a coordinate. Convert 1-letter
to 3-letter amino acids; "X" or "*" at a position is a stop codon. A trailing
asterisk AFTER a name ("D16A*") marks a NOVEL variant, not a stop - say so in
additional_notes. Record source_notation verbatim as the paper wrote it.

PROVENANCE - required whenever any count is populated
count_provenance records WHY a count was assigned:
  carriers_column_label / affected_column_label / unaffected_column_label - the raw
  column header text, or null. Do not invent labels.
  carriers_count_type / affected_count_type / unaffected_count_type - exactly one
  of: per_variant_carrier (the count to actually use) | family_count |
  proband_count | cohort_total | screened_N | case | control | unaffected_control |
  unknown.
Rule: when a count_type is cohort_total, screened_N, or unknown AND the value would
be large (>10x the smallest carrier count in this paper), leave that count NULL and
record the aggregate in additional_notes instead.

fact_provenance: emit one row per extracted variant identity, carrier count,
affected count, unaffected count, and individual affected/unaffected status, with
the most exact location available - fact_type, fact_value, individual_id,
source_location, source_table ("Table S1"), source_row ("row 20"), source_column
("Affected"), source_section ("Methods"), source_paragraph ("paragraph 2"), and a
short verbatim evidence_quote.

Per-observation locators: on every aggregate patients object and every
individual_records entry, fill what the paper provides and leave the rest null -
source_container ("main" or "supplement"; a supplementary table is still a table),
source_kind ("table"/"figure"/"text"/"abstract"), source_ref ("Table 2",
"Figure 3B"), page_label, pdf_page, row_label ("II-1"), row_ordinal, column_ref,
figure_panel, locator_extra. Do not invent page, row, column, or panel values. If
the paper gives only aggregate counts, emit one grouped variant record carrying the
aggregate locator. The SQLite writer computes source_record_id; do not hash it."""
)

CONTINUATION_PROMPT = """You previously extracted variants from this paper but the response was truncated.
You extracted {extracted_count} variants so far. The paper contains approximately {expected_count} variants total.

Previously extracted variants (DO NOT re-extract these):
{extracted_variants_list}

Please continue extracting the REMAINING variants starting AFTER the last one listed above.
Return ONLY the variants you haven't extracted yet in the same JSON format.

TARGET GENE: {gene_symbol}
Paper Title: {title}

Full Text:
{full_text}

Return a JSON object with this structure:
{{
    "continuation": true,
    "variants": [
        ... remaining variants only ...
    ],
    "extraction_metadata": {{
        "continuation_variants_found": integer,
        "notes": "any notes about this continuation"
    }}
}}
"""

COMPACT_EXTRACTION_PROMPT = (
    """You are an expert medical geneticist extracting {gene_symbol} variants from a
paper that contains MANY variants (estimated {estimated_variants}+). Extract ALL of
them in the COMPACT format below: completeness beats detail, so skip
individual_records, functional_data, key_quotes, and age_dependent_penetrance
rather than dropping variants. Do NOT stop early.

A wrong number is worse than a missing one. When support is not explicit, use null.

TARGET GENE: {gene_symbol} - ignore variants in every other gene.
Paper Title: {title}

Full Text:
{full_text}

"""
    + _CORE_RULES
    + """

OUTPUT - JSON only:
{{
    "paper_metadata": {{"pmid": "{pmid}", "title": "{title}",
        "extraction_summary": "Compact extraction of {gene_symbol} variants"}},
    "variants": [
        {{
            "gene_symbol": "{gene_symbol}",
            "cdna_notation": "c.XXX or IVS... or null",
            "protein_notation": "p.XXX or null",
            "source_notation": "verbatim as printed, or null",
            "variant_class": "missense|nonsense|frameshift|inframe_indel|splice|deep_intronic|large_deletion|large_duplication|cnv|exon_deletion|exon_duplication|complex|other or null",
            "structural_description": "e.g. 'deletion of exons 3-5' or null",
            "clinical_significance": "pathogenic|likely_pathogenic|VUS|likely_benign|benign",
            "patients": {{"count": N, "phenotype": "brief",
                "source_container": "main|supplement|null", "source_kind": "table|figure|text|abstract|null",
                "source_ref": "Table X or null", "page_label": null, "pdf_page": null,
                "row_label": null, "row_ordinal": null, "column_ref": null,
                "figure_panel": null, "locator_extra": {{}}}},
            "penetrance_data": {{"total_carriers_observed": N or null,
                "affected_count": N or null, "unaffected_count": N or null}},
            "count_provenance": {{
                "carriers_column_label": "string or null",
                "carriers_count_type": "per_variant_carrier|family_count|proband_count|cohort_total|screened_N|case|control|unaffected_control|unknown",
                "affected_column_label": "string or null", "affected_count_type": "(same enum)",
                "unaffected_column_label": "string or null", "unaffected_count_type": "(same enum)"}},
            "source_location": "Table X or Results",
            "fact_provenance": [
                {{"fact_type": "variant_identity|patient_count|total_carriers_observed|affected_count|unaffected_count|individual_affected_status",
                  "fact_value": "value", "individual_id": "string or null",
                  "source_location": "string", "source_table": "string or null",
                  "source_row": "string or null", "source_column": "string or null",
                  "source_section": "string or null", "source_paragraph": "string or null",
                  "evidence_quote": "short exact quote"}}
            ]
        }}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high|medium|low",
        "study_type": "clinical|functional|mixed (REQUIRED)",
        "study_design": "case_report|case_series|case_control|cohort_population|cohort_biobank|family_segregation|functional_invitro|gwas|review_meta|other (REQUIRED)",
        "study_summary": "1-3 sentences on what this study is",
        "compact_mode": true,
        "notes": "string"
    }}
}}

Structural events (exon deletions, large del/dup, CNVs) may omit point-form
notation when variant_class + structural_description are set.
"""
)

EXTRACTION_PROMPT = (
    """You are an expert medical geneticist extracting genetic variants from a
scientific paper, with emphasis on penetrance (affected vs unaffected carriers).

Extract ALL variants - completeness beats per-variant detail. If the paper holds
many variants, keep minimal fields each rather than dropping any; above ~50
variants, at most one individual_records entry per variant.

A wrong number is worse than a missing one. When support is not explicit, use null.

TARGET GENE: {gene_symbol} - only this gene. If the paper mentions it but reports
no variants, return an empty variants list.
Paper Title: {title}

Full Text (including tables):
{full_text}

"""
    + _CORE_RULES
    + """

ALSO CAPTURE, when the paper provides it
- Individuals: EVERY person carrying the variant (proband, case, patient, subject,
  family member, or a label like II-1 / P1 / Case 2), with age at evaluation, onset,
  and diagnosis; sex; affected status; their own phenotype; reported
  ethnicity/ancestry and geographic origin verbatim; and the exact sentence. A young
  unaffected carrier may simply not be past the risk window - say so rather than
  reclassifying.
- Variant-specific cohort penetrance: total carriers, affected, unaffected,
  uncertain, and age-stratified penetrance, only when explicitly tied to ONE
  variant.
- Functional results, segregation, population frequency, evidence level, and
  demographics (age, sex, ancestry, country of origin - verbatim).

ABSTRACT-ONLY SOURCE (marked "[ABSTRACT ONLY - FULL TEXT NOT AVAILABLE]")
Extract what is there, including a variant notation on its own; apply the same
implicit-count and abstain rules; set extraction_confidence low or medium and note
it. Do not skip extraction just because the information is thin.

OUTPUT - JSON only:
{{
    "paper_metadata": {{"pmid": "{pmid}", "title": "{title}",
        "extraction_summary": "brief summary of what was extracted"}},
    "variants": [
        {{
            "gene_symbol": "string",
            "cdna_notation": "string or null (c. or IVS notation)",
            "protein_notation": "string or null",
            "source_notation": "the variant EXACTLY as written in this paper, verbatim, before normalization (e.g. 'IVS3+2T>G', 'R1443X', '5382insC', 'del exon 13'), or null",
            "genomic_position": "string or null",
            "variant_class": "missense|nonsense|frameshift|inframe_indel|splice|deep_intronic|large_deletion|large_duplication|cnv|exon_deletion|exon_duplication|complex|other or null",
            "structural_description": "e.g. 'deletion of exons 3-5' or null",
            "clinical_significance": "string",
            "patients": {{"count": integer, "demographics": "string", "phenotype": "string",
                "source_container": "main|supplement|null", "source_kind": "table|figure|text|abstract|null",
                "source_ref": "Table X or Figure Y or Results or null",
                "page_label": "string or null", "pdf_page": "integer or null",
                "row_label": "string or null", "row_ordinal": "integer or null",
                "column_ref": "string or null", "figure_panel": "string or null",
                "locator_extra": {{}}}},
            "penetrance_data": {{
                "total_carriers_observed": "integer or null",
                "affected_count": "integer or null",
                "unaffected_count": "integer or null",
                "uncertain_count": "integer or null (unclear status or too young)",
                "penetrance_percentage": "float or null; ONLY when the paper explicitly states a variant-specific percentage; never calculate it from counts",
                "age_dependent_penetrance": [
                    {{"age_range": "e.g. '40-50 years'", "penetrance_percentage": "float",
                      "carriers_in_range": "integer", "affected_in_range": "integer",
                      "evidence_quote": "exact source quote containing the stated percentage",
                      "source_location": "table/figure/section or null"}}
                ]}},
            "count_provenance": {{
                "carriers_column_label": "string or null (raw column header the count came from)",
                "carriers_count_type": "per_variant_carrier|family_count|proband_count|cohort_total|screened_N|case|control|unaffected_control|unknown",
                "affected_column_label": "string or null", "affected_count_type": "(same enum)",
                "unaffected_column_label": "string or null", "unaffected_count_type": "(same enum)"}},
            "individual_records": [
                {{"individual_id": "e.g. 'II-1', 'P1', 'Case_2'",
                  "age_at_evaluation": "integer or null", "age_at_onset": "integer or null",
                  "age_at_diagnosis": "integer or null", "sex": "male|female|other|null",
                  "affected_status": "affected|unaffected|uncertain",
                  "phenotype_details": "string", "ethnicity": "string or null",
                  "geographic_origin": "string or null",
                  "evidence_sentence": "exact sentence from paper",
                  "source_container": "main|supplement|null", "source_kind": "table|figure|text|abstract|null",
                  "source_ref": "string or null", "page_label": "string or null",
                  "pdf_page": "integer or null", "row_label": "string or null",
                  "row_ordinal": "integer or null", "column_ref": "string or null",
                  "figure_panel": "string or null", "locator_extra": {{}}}}
            ],
            "functional_data": {{"summary": "string", "assays": ["list of assays"]}},
            "segregation_data": "string or null",
            "population_frequency": "string or null",
            "evidence_level": "strong|moderate|weak|supporting",
            "source_location": "e.g. 'Table 2, Row 3' or 'Results, paragraph 4'",
            "additional_notes": "string",
            "key_quotes": ["relevant quotes from paper"],
            "fact_provenance": [
                {{"fact_type": "variant_identity|patient_count|total_carriers_observed|affected_count|unaffected_count|penetrance_percentage|individual_affected_status",
                  "fact_value": "value", "individual_id": "string or null",
                  "source_location": "string", "source_table": "string or null",
                  "source_row": "string or null", "source_column": "string or null",
                  "source_section": "string or null", "source_paragraph": "string or null",
                  "evidence_quote": "short exact quote"}}
            ]
        }}
    ],
    "tables_processed": [
        {{"table_name": "e.g. 'Table 1', 'Supplementary Table 3'",
          "table_caption": "string", "variants_extracted": integer}}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high|medium|low",
        "study_type": "clinical|functional|mixed (REQUIRED)",
        "study_design": "case_report|case_series|case_control|cohort_population|cohort_biobank|family_segregation|functional_invitro|gwas|review_meta|other (REQUIRED)",
        "ascertainment": "proband_referral|population_screening|biobank|family_cascade|unknown (drives penetrance interpretation)",
        "cohort_source": "e.g. 'LQTS referral clinic, Japan'; 'UK Biobank' or null",
        "population": "study-level ancestry/geography/founder population or null",
        "study_summary": "1-3 sentences on what this study is",
        "challenges": ["any issues during extraction"],
        "notes": "string"
    }}
}}

Preserve the paper's exact nomenclature. Do not invent data that is not there.
"""
)
