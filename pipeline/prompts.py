"""
LLM Prompt Templates for Variant Extraction

This module centralizes all LLM prompts used in the extraction pipeline,
making them easier to maintain, test, and iterate on.
"""

# Threshold for switching to compact extraction mode
HIGH_VARIANT_THRESHOLD = 30

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

COMPACT_EXTRACTION_PROMPT = """You are an expert medical geneticist extracting genetic variants from a scientific paper.

This paper contains MANY variants (estimated {estimated_variants}+). Use COMPACT output format.

TARGET GENE: {gene_symbol}
Paper Title: {title}

Full Text:
{full_text}

INSTRUCTIONS:
1. Extract ALL {gene_symbol} variants - completeness is critical
2. Use MINIMAL fields per variant to fit all variants in the response
3. Skip detailed fields (individual_records, functional_data, key_quotes, age_dependent_penetrance)

CRITICAL - FUNCTIONAL STUDY DETECTION:
Before extracting, determine if this is a FUNCTIONAL STUDY vs CLINICAL STUDY:

FUNCTIONAL STUDY indicators (DO NOT extract as patient carriers):
- Tables with columns like "number of cells", "n (cells)", "replicates", "n="
- Assay metrics: "trafficking", "current density", "fluorescence", "normalized"
- Study types: "massively parallel assay", "high-throughput", "saturation mutagenesis", "functional characterization"
- Systematic testing of all possible variants at specific positions
- Sample sizes are ASSAY REPLICATES, not patients

If this is a FUNCTIONAL STUDY:
- The "n" or sample size column = cells/replicates tested, NOT patient carriers
- Set penetrance_data counts to null (do not use assay replicates as carrier counts)
- Look for SEPARATE columns labeled "LQTS", "patients", "carriers", "cases", "probands" for actual clinical data
- Only extract penetrance_data from columns explicitly describing patient/carrier counts
- Note in source_location: "Functional study - assay data only"

CLINICAL STUDY indicators (extract normally):
- Case reports, cohort studies, family studies
- Columns like "patients", "carriers", "probands", "cases", "affected", "unaffected"
- Individual patient descriptions with ages, symptoms, phenotypes

RECOGNIZING IMPLICIT COUNTS (for CLINICAL studies only):
- "a patient" / "a case" / "an individual" = 1 carrier
- "a healthy individual" / "an asymptomatic carrier" = 1 UNAFFECTED carrier
- "an affected patient" / "a patient with [disease]" = 1 AFFECTED carrier
- Disease-associated mutation tables: patient counts are typically AFFECTED carriers

CRITICAL — COHORT/SCREENING TABLE CARRIER COUNTS:
Large cohort or genetic screening papers often report the SAME variant across multiple
tables, cohorts, or registries. Capture the TOTAL carrier count across all sources,
not just the first table. If variant X appears in Table 1 with 50 carriers and Table 2
with 78 carriers from a non-overlapping cohort, total = 128. "N patients" or "No. of
patients" columns are CARRIER COUNTS. Frequency columns like "5/200" → 5 carriers.
When same variant in multiple tables: SUM if independent, use LARGER if overlapping.

CRITICAL — DO NOT COPY STUDY-WIDE TOTALS ONTO EACH VARIANT:
If the paper states an aggregate such as "43 carriers, 28 affected and 15 unaffected"
for the whole study, cohort, family set, domain, or mutation class, do NOT assign
43/28/15 to every variant. Use counts only when the row, sentence, or table cell is
variant-specific. If only aggregate counts are available, leave variant-level counts
null and put the aggregate in additional_notes.

This rule also applies to patients.phenotype, patients.demographics, and
penetrance_percentage. Study-wide statements such as "35% of carriers overall"
or "median age 61 years" are not variant-level facts. Do not copy them onto
every variant; mention them only in additional_notes when useful.

REQUIRED — COUNT PROVENANCE (records WHY a count was assigned):
For every count you populate (patients.count, penetrance_data.*), also record
where it came from in a separate `count_provenance` block on the variant:

  count_provenance: {{
      "carriers_column_label":   raw column header text or null
                                 (e.g. "No. of patients", "Carriers", "N")
      "carriers_count_type":     one of the count_type values below
      "affected_column_label":   raw column header text or null
      "affected_count_type":     one of the count_type values below
      "unaffected_column_label": raw column header text or null
      "unaffected_count_type":   one of the count_type values below
  }}

count_type values (use exactly these strings):
  per_variant_carrier — this row's number specifically describes carriers
                        of THIS variant (the count to actually use)
  family_count        — families, not individuals
  proband_count       — probands only, not all carriers
  cohort_total        — study/cohort total ("N studied = 200"); NOT variant-specific
  screened_N          — screening denominator; NOT a carrier count
  case                — case-arm count in a case/control study
  control             — control-arm count
  unaffected_control  — healthy-control count
  unknown             — cannot tell from context

Rule: if count_type is cohort_total, screened_N, or unknown AND the count would
otherwise be large (>10× the smallest carrier count in this paper), leave the
corresponding patients.count / penetrance_data.* NULL — do NOT assign the cohort
total to this row. Record the aggregate in additional_notes instead.

This block is REQUIRED whenever any count is populated. Set fields to null when
you can't identify the source column. Do NOT invent column labels.

REQUIRED — REJECT COHORT-CLASS SUMMARIES (no specific variant identifier):
Do NOT emit a variant entry where ALL THREE of cdna_notation, protein_notation,
and genomic_position would be null. The variant identifier is load-bearing.

REJECT examples:
- "34 patients with pore-region mutations" (mutation class, not a variant)
- "Subjects with N-terminal mutations, 4 affected" (region, not a variant)
- Aggregating multiple variants by domain/class without naming each one

ACCEPT examples (at least one identifier present):
- protein_notation: "p.Arg176Trp" → emit
- cdna_notation: "c.1234G>A" → emit
- genomic_position: "chr7:150646404G>A" → emit

If the paper describes a cohort by mutation class without listing the variants,
return ZERO entries for that cohort description rather than a placeholder.

REQUIRED — REJECT TABLE-CELL ARTIFACTS (single letters, gene symbols, NaN):
Do NOT emit a variant entry where the notation is a fragment of a table cell
without enough structure to be a real variant. These come from genotype call
columns, allele frequency tables, and SNP arrays — they are NOT variants.

REJECT examples (all of these have appeared in past extractions):
- protein_notation: "A" / "T" / "C" / "G" (single nucleotide letter)
- protein_notation: "R" / "K" / "Y" (single amino acid code with no position)
- protein_notation: "nan" / "NaN" / "NA" (parser/Python NaN literal — not a variant)
- protein_notation: "KCNH2" or any gene symbol (the gene name is not a variant ID)
- cdna_notation: "c.nan" / "c.KCNH2" (concatenation artifacts)
- protein_notation: "Likelybenign/benign" / "pathogenic" (clinical_significance leaked
  into the wrong field — put those values in clinical_significance, not in a notation
  field)

A valid protein_notation has the form `(p.)?<AA1-or-3-letter><pos>(_<AA><pos>)?<change>`.
A valid cdna_notation starts with `c.` followed by a digit. A valid genomic_position
contains a chromosome and a coordinate. If you cannot produce one of those forms, do
NOT emit the entry — return zero entries for that row instead.

OUTPUT FORMAT (compact):
Return a JSON object with this structure:
{{
    "paper_metadata": {{
        "pmid": "{pmid}",
        "title": "{title}",
        "extraction_summary": "Compact extraction of {gene_symbol} variants"
    }},
    "variants": [
        {{
            "gene_symbol": "{gene_symbol}",
            "cdna_notation": "c.XXX",
            "protein_notation": "p.XXX",
            "clinical_significance": "pathogenic/likely_pathogenic/VUS/likely_benign/benign",
            "patients": {{"count": N, "phenotype": "brief description"}},
            "penetrance_data": {{
                "total_carriers_observed": N or null,
                "affected_count": N or null,
                "unaffected_count": N or null
            }},
            "count_provenance": {{
                "carriers_column_label": "string or null",
                "carriers_count_type": "per_variant_carrier|family_count|proband_count|cohort_total|screened_N|case|control|unaffected_control|unknown",
                "affected_column_label": "string or null",
                "affected_count_type": "(same enum)",
                "unaffected_column_label": "string or null",
                "unaffected_count_type": "(same enum)"
            }},
            "source_location": "Table X" or "Results"
        }}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high/medium/low",
        "study_type": "clinical/functional/mixed (REQUIRED - indicate if this is a functional assay study)",
        "compact_mode": true,
        "notes": "Compact extraction - detailed fields omitted for completeness"
    }}
}}

CRITICAL: Extract ALL variants. Do NOT stop early. Completeness > detail.
"""

EXTRACTION_PROMPT = """You are an expert medical geneticist and data extraction specialist. Your task is to extract genetic variant information from the provided scientific paper, with special emphasis on penetrance data (affected vs unaffected carriers).

RESPONSE SIZE GUIDANCE:
- You MUST extract ALL variants found in the paper - do not stop early or truncate the list
- For papers with many variants (large tables), use a simplified format per variant:
  * Include: gene_symbol, cdna_notation, protein_notation, clinical_significance, patients.count, source_location
  * Use empty arrays for: individual_records, key_quotes, functional_data.assays, age_dependent_penetrance
  * Use null for: genomic_position, penetrance_data counts (unless explicitly stated), segregation_data, population_frequency
  * Use brief strings for: patients.phenotype, additional_notes
- If there are >50 variants, limit individual_records to at most 1 per variant (or empty list)
- Prioritize completeness over detail: extract ALL variants with minimal fields rather than few variants with full detail

TARGET GENE: {gene_symbol}

CRITICAL: Only extract variants in the gene "{gene_symbol}". Ignore all variants in other genes.

Paper Title: {title}

Full Text (including tables):
{full_text}

EXTRACTION INSTRUCTIONS:
Extract ALL variants in the {gene_symbol} gene mentioned in this paper with the following structured information:

For each variant, provide:
1. Gene Symbol (e.g., "BRCA1", "TP53")
2. Variant Notation:
   - cDNA notation (e.g., "c.1234G>A") - convert from nucleotide position if needed (e.g., "47 A>C" → "c.47A>C")
   - Protein notation (e.g., "p.Arg412His") - use full HGVS notation:
     * Single-letter amino acids should be converted to three-letter (e.g., "D16A" → "p.Asp16Ala")
     * Nonsense/stop mutations: "X" or "*" in protein position means stop codon (e.g., "R412X" → "p.Arg412*")
     * Frameshift: include fs and stop position (e.g., "G24fs+34X" → "p.Gly24fs*58" or "p.Gly24fs")
   - Genomic position (if available)

   IMPORTANT - Novel mutation markers:
   - An asterisk (*) AFTER a mutation name (e.g., "D16A*", "R20G*") typically indicates "novel mutation"
     (first reported in this study). This is NOT a stop codon - note it in additional_notes as "novel mutation".
   - An asterisk (*) IN the protein position (e.g., "R412*", "p.Arg412*") indicates a stop codon/nonsense mutation.

3. Clinical Significance: pathogenic, likely pathogenic, benign, likely benign, VUS, etc.
4. Patient Information:
   - Number of patients/cases
   - Demographics (age, sex, ethnicity if mentioned)
   - Phenotype/disease presentation
5. Functional Data:
   - In vitro studies
   - Functional assays
   - Protein effects
6. Segregation Data: Does it segregate with disease in families?
7. Population Frequency: Allele frequency in databases (gnomAD, ExAC, etc.)
8. Evidence Level: How strong is the evidence for pathogenicity?
9. Additional Notes: Any other relevant clinical or functional information

CRITICAL - FUNCTIONAL STUDY DETECTION:
Before extracting penetrance data, determine if this is a FUNCTIONAL STUDY vs CLINICAL STUDY:

FUNCTIONAL STUDY indicators (DO NOT extract as patient carriers):
- Tables with columns like "number of cells", "n (cells)", "replicates", "n="
- Assay metrics: "trafficking", "current density", "fluorescence", "normalized", "patch clamp"
- Study types: "massively parallel assay", "high-throughput", "saturation mutagenesis", "functional characterization"
- Systematic testing of all possible variants at specific positions (e.g., all SNVs in an exon)
- Sample sizes are ASSAY REPLICATES, not patients

If this is a FUNCTIONAL STUDY:
- The "n" or sample size column = cells/replicates tested, NOT patient carriers
- Set penetrance_data counts to null (do not use assay replicates as carrier counts)
- Look for SEPARATE columns labeled "LQTS", "patients", "carriers", "cases", "probands" for actual clinical data
- Only extract penetrance_data from columns explicitly describing patient/carrier counts
- Extract functional metrics to functional_data.assays instead
- Note in source_location: "Functional study - assay data only"

CLINICAL STUDY indicators (extract penetrance normally):
- Case reports, cohort studies, family studies, registry data
- Columns explicitly labeled "patients", "carriers", "probands", "cases", "affected", "unaffected"
- Individual patient descriptions with ages, symptoms, phenotypes
- Pedigrees with affected/unaffected family members

SANITY CHECK: If you find 50+ variants all with similar carrier counts (e.g., 20-50 each) and ~100% penetrance,
this is likely a functional study where you're mistaking assay replicates for patients. Re-evaluate the table structure.

PENETRANCE DATA EXTRACTION (CRITICAL - for CLINICAL studies only):
For calibrating disease prediction models, extract detailed penetrance information:

RECOGNIZING IMPLICIT COUNTS (IMPORTANT - for CLINICAL studies):
Natural language often implies counts without stating explicit numbers. You MUST recognize these:
- "a patient" / "a case" / "an individual" / "one patient" = 1 carrier
- "a healthy individual" / "an asymptomatic carrier" = 1 UNAFFECTED carrier
- "an affected patient" / "a patient with [disease]" = 1 AFFECTED carrier
- "proband" / "the proband" = 1 carrier (usually affected unless stated otherwise)
- "two patients" / "three families" = 2 or 3 carriers
- For case reports mentioning a mutation, assume at least 1 affected carrier

AFFECTED vs UNAFFECTED:
- "Affected" = has the disease phenotype, symptoms, or clinical manifestations
- "Unaffected" = carrier WITHOUT disease (described as: healthy, asymptomatic, no symptoms,
  clinically silent, normal phenotype, unaffected family member)
- Probands and "patients" in disease studies are typically AFFECTED
- Family members described as "healthy" or "asymptomatic" carriers are UNAFFECTED

A. Individual-Level Records:
   - Extract EVERY individual person mentioned who carries the variant
   - Look for: proband, case, patient, subject, individual, family member, sibling, parent, or designations like II-1, III-2, P1, Case 2, etc.
   - For each individual, capture:
     * Age at evaluation/assessment (if mentioned)
     * Age at disease onset (if mentioned)
     * Age at diagnosis (if mentioned)
     * Sex
     * Affected status: "affected", "unaffected", or "uncertain"
     * Phenotype details specific to that individual
     * Exact sentence where this information appears
   - Use age-appropriate penetrance logic: a young person may be "unaffected" but not past the risk window

B. Variant-Specific Cohort/Aggregate Penetrance Data:
   - Extract cohort statistics only when explicitly tied to one variant
     (e.g., "10 carriers of variant X, 4 affected, 6 unaffected")
   - Total carriers observed for that variant
   - Affected count for that variant
   - Unaffected count for that variant
   - Uncertain/unclear cases for that variant
   - Age-dependent penetrance data if stratified by age groups for that variant
   - Variant-specific penetrance percentages if explicitly stated
   - If the statistic describes the whole study/cohort/table/domain rather than
     one variant, leave the variant-level count and penetrance fields null

C. Age-Dependent Penetrance:
   - Capture penetrance stratified by age ranges if mentioned
   - Note age at which penetrance was assessed
   - For age-dependent diseases, track penetrance by age group

CRITICAL REQUIREMENTS:
- Extract data from BOTH the main text AND all tables
- Pay special attention to supplementary table references
- For tables: int ALL rows with variant data AND individual person data
- If a variant is mentioned multiple times, consolidate the information but preserve individual-level detail
- Include exact quotes for key clinical descriptions
- Note the specific section/table where each variant was found
- Distinguish between individual case reports and cohort study statistics
- When paper states "10 people with variant X, 4 with disease" → extract: total_carriers=10, affected=4, unaffected=6
- When individual cases are described, create individual_records entries for each person
- IMPORTANT: For CLINICAL disease-associated mutation tables (e.g., "LQT2-associated mutations", "pathogenic variants"):
  * Patient counts in these tables represent AFFECTED carriers
  * Extract patient count to BOTH patients.count AND penetrance_data.total_carriers_observed/affected_count
  * Example: Table shows "No. of patients: 1" for a pathogenic variant → patients.count=1, total_carriers_observed=1, affected_count=1
  * BUT: If table has "n (cells)" or assay metrics, it's a FUNCTIONAL study - do NOT use those numbers as patient counts

CRITICAL — COHORT/SCREENING TABLE CARRIER COUNTS:
Large cohort or genetic screening papers often report the SAME variant across multiple
tables, cohorts, or registries. You MUST capture the TOTAL carrier count across all
sources within the paper, not just the first table you encounter.

Common patterns in cohort papers:
- Table 1: "Discovery cohort" (e.g., 50 carriers of R176W)
- Table 2: "Replication cohort" (e.g., 78 more carriers of R176W)
- → Total: 128 carriers. Do NOT report just 50 from the first table.

- "N patients" or "No. of patients" columns are CARRIER COUNTS, not replicates
- "n" in clinical tables (not functional assays) = number of patients/carriers
- Frequency columns (e.g., "5/200") give carrier counts: numerator = carriers, denominator = total screened
- "Families" ≠ "carriers": if "3 families" are reported, look for how many total carriers across those families
- If the paper says "identified in N individuals" or "found in N patients", that is the total carrier count

When the same variant appears in multiple tables or sections:
1. SUM the carrier counts if the cohorts are independent/non-overlapping
2. Use the LARGER number if the cohorts overlap or one is a subset
3. When in doubt, use the largest reported count (fail toward completeness)

CRITICAL — DO NOT COPY STUDY-WIDE TOTALS ONTO EACH VARIANT:
If the paper states an aggregate such as "43 carriers, 28 affected and 15 unaffected"
for the whole study, cohort, family set, domain, or mutation class, do NOT assign
43/28/15 to every variant. Use counts only when the row, sentence, or table cell is
variant-specific. If only aggregate counts are available, leave variant-level counts
null and put the aggregate in additional_notes.

This rule also applies to patients.phenotype, patients.demographics, and
penetrance_percentage. Study-wide statements such as "35% of carriers overall"
or "median age 61 years" are not variant-level facts. Do not copy them onto
every variant; mention them only in additional_notes when useful.

REQUIRED — COUNT PROVENANCE (records WHY a count was assigned):
For every count you populate (patients.count, penetrance_data.*), also record
where it came from in a separate `count_provenance` block on the variant:

  count_provenance: {{
      "carriers_column_label":   raw column header text or null
                                 (e.g. "No. of patients", "Carriers", "N")
      "carriers_count_type":     one of the count_type values below
      "affected_column_label":   raw column header text or null
      "affected_count_type":     one of the count_type values below
      "unaffected_column_label": raw column header text or null
      "unaffected_count_type":   one of the count_type values below
  }}

count_type values (use exactly these strings):
  per_variant_carrier — this row's number specifically describes carriers
                        of THIS variant (the count to actually use)
  family_count        — families, not individuals
  proband_count       — probands only, not all carriers
  cohort_total        — study/cohort total ("N studied = 200"); NOT variant-specific
  screened_N          — screening denominator; NOT a carrier count
  case                — case-arm count in a case/control study
  control             — control-arm count
  unaffected_control  — healthy-control count
  unknown             — cannot tell from context

Rule: if count_type is cohort_total, screened_N, or unknown AND the count would
otherwise be large (>10× the smallest carrier count in this paper), leave the
corresponding patients.count / penetrance_data.* NULL — do NOT assign the cohort
total to this row. Record the aggregate in additional_notes instead.

This block is REQUIRED whenever any count is populated. Set fields to null when
you can't identify the source column. Do NOT invent column labels.

REQUIRED — REJECT COHORT-CLASS SUMMARIES (no specific variant identifier):
Do NOT emit a variant entry where ALL THREE of cdna_notation, protein_notation,
and genomic_position would be null. The variant identifier is load-bearing.

REJECT examples:
- "34 patients with pore-region mutations" (mutation class, not a variant)
- "Subjects with N-terminal mutations, 4 affected" (region, not a variant)
- Aggregating multiple variants by domain/class without naming each one

ACCEPT examples (at least one identifier present):
- protein_notation: "p.Arg176Trp" → emit
- cdna_notation: "c.1234G>A" → emit
- genomic_position: "chr7:150646404G>A" → emit

If the paper describes a cohort by mutation class without listing the variants,
return ZERO entries for that cohort description rather than a placeholder.

REQUIRED — REJECT TABLE-CELL ARTIFACTS (single letters, gene symbols, NaN):
Do NOT emit a variant entry where the notation is a fragment of a table cell
without enough structure to be a real variant. These come from genotype call
columns, allele frequency tables, and SNP arrays — they are NOT variants.

REJECT examples (all of these have appeared in past extractions):
- protein_notation: "A" / "T" / "C" / "G" (single nucleotide letter)
- protein_notation: "R" / "K" / "Y" (single amino acid code with no position)
- protein_notation: "nan" / "NaN" / "NA" (parser/Python NaN literal — not a variant)
- protein_notation: "KCNH2" or any gene symbol (the gene name is not a variant ID)
- cdna_notation: "c.nan" / "c.KCNH2" (concatenation artifacts)
- protein_notation: "Likelybenign/benign" / "pathogenic" (clinical_significance leaked
  into the wrong field — put those values in clinical_significance, not in a notation
  field)

A valid protein_notation has the form `(p.)?<AA1-or-3-letter><pos>(_<AA><pos>)?<change>`.
A valid cdna_notation starts with `c.` followed by a digit. A valid genomic_position
contains a chromosome and a coordinate. If you cannot produce one of those forms, do
NOT emit the entry — return zero entries for that row instead.

OUTPUT FORMAT:
Return a JSON object with this structure:
{{
    "paper_metadata": {{
        "pmid": "{pmid}",
        "title": "{title}",
        "extraction_summary": "Brief summary of what was extracted"
    }},
    "variants": [
        {{
            "gene_symbol": "string",
            "cdna_notation": "string or null",
            "protein_notation": "string or null",
            "genomic_position": "string or null",
            "clinical_significance": "string",
            "patients": {{
                "count": integer,
                "demographics": "string",
                "phenotype": "string"
            }},
            "penetrance_data": {{
                "total_carriers_observed": "integer or null (total people with variant)",
                "affected_count": "integer or null (number with disease)",
                "unaffected_count": "integer or null (number without disease, past risk window)",
                "uncertain_count": "integer or null (unclear status or too young)",
                "penetrance_percentage": "float or null (calculated: affected/total_carriers * 100)",
                "age_dependent_penetrance": [
                    {{
                        "age_range": "string (e.g., '40-50 years')",
                        "penetrance_percentage": "float",
                        "carriers_in_range": "integer",
                        "affected_in_range": "integer"
                    }}
                ]
            }},
            "count_provenance": {{
                "carriers_column_label": "string or null (raw column header the count came from)",
                "carriers_count_type": "per_variant_carrier|family_count|proband_count|cohort_total|screened_N|case|control|unaffected_control|unknown",
                "affected_column_label": "string or null",
                "affected_count_type": "(same enum as carriers_count_type)",
                "unaffected_column_label": "string or null",
                "unaffected_count_type": "(same enum as carriers_count_type)"
            }},
            "individual_records": [
                {{
                    "individual_id": "string (e.g., 'II-1', 'P1', 'Case_2', or generate unique ID)",
                    "age_at_evaluation": "integer or null",
                    "age_at_onset": "integer or null (age when symptoms started)",
                    "age_at_diagnosis": "integer or null (age when diagnosed)",
                    "sex": "string or null (male/female/other)",
                    "affected_status": "string (affected/unaffected/uncertain)",
                    "phenotype_details": "string (disease manifestations for this person)",
                    "ethnicity": "string or null (reported race/ethnicity/ancestry of this individual or their cohort, verbatim, e.g. 'East Asian', 'Caucasian', 'Ashkenazi Jewish')",
                    "geographic_origin": "string or null (reported country/region/population of origin, e.g. 'Japan', 'Northern Italy')",
                    "evidence_sentence": "string (exact sentence from paper)"
                }}
            ],
            "functional_data": {{
                "summary": "string",
                "assays": ["list of assays performed"]
            }},
            "segregation_data": "string or null",
            "population_frequency": "string or null",
            "evidence_level": "string",
            "source_location": "e.g., 'Table 2, Row 3' or 'Results, paragraph 4'",
            "additional_notes": "string",
            "key_quotes": ["relevant quotes from paper"]
        }}
    ],
    "tables_processed": [
        {{
            "table_name": "string (e.g., 'Table 1', 'Supplementary Table 3')",
            "table_caption": "string",
            "variants_extracted": integer
        }}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high/medium/low",
        "study_type": "clinical/functional/mixed (REQUIRED - indicate if this is a functional assay study)",
        "challenges": ["any issues during extraction"],
        "notes": "any additional notes about the extraction process"
    }}
}}

IMPORTANT NOTES:
- ONLY include variants in the {gene_symbol} gene. Do NOT include variants from other genes even if they are mentioned in the paper.
- If the paper mentions the target gene {gene_symbol} but does not report any variants, return an empty variants list.
- Be thorough but accurate - don't invent data not present in the paper
- If a field is not available, use null or an empty string as appropriate
- Preserve exact nomenclature from the paper

ABSTRACT-ONLY EXTRACTION:
When only an abstract is available (indicated by "[ABSTRACT ONLY - FULL TEXT NOT AVAILABLE]"):
- Extract whatever variant and carrier information is available, even if limited
- Recognize implicit counts: "a patient with mutation X" = 1 affected carrier with variant X
- Case reports typically describe at least 1 affected carrier
- "A novel mutation" or "we report a mutation" implies at least 1 carrier
- Even if carrier counts aren't explicit, extract the variant notation if mentioned
- Set extraction_confidence to "low" or "medium" due to limited information
- Note in extraction_metadata.notes that this was abstract-only extraction
- Do NOT skip extraction just because information is limited - extract what's available
"""
