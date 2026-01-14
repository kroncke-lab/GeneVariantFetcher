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
            "source_location": "Table X" or "Results"
        }}
    ],
    "extraction_metadata": {{
        "total_variants_found": integer,
        "extraction_confidence": "high/medium/low",
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

PENETRANCE DATA EXTRACTION (CRITICAL):
For calibrating disease prediction models, extract detailed penetrance information:

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

B. Cohort/Aggregate Penetrance Data:
   - Extract study-level statistics when provided (e.g., "10 carriers, 4 affected, 6 unaffected")
   - Total carriers observed
   - Affected count
   - Unaffected count
   - Uncertain/unclear cases
   - Age-dependent penetrance data if stratified by age groups
   - Penetrance percentages if explicitly stated

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
- IMPORTANT: For disease-associated mutation tables (e.g., "LQT2-associated mutations", "pathogenic variants"):
  * Patient counts in these tables represent AFFECTED carriers
  * Extract patient count to BOTH patients.count AND penetrance_data.total_carriers_observed/affected_count
  * Example: Table shows "No. of patients: 1" for a pathogenic variant → patients.count=1, total_carriers_observed=1, affected_count=1

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
            "individual_records": [
                {{
                    "individual_id": "string (e.g., 'II-1', 'P1', 'Case_2', or generate unique ID)",
                    "age_at_evaluation": "integer or null",
                    "age_at_onset": "integer or null (age when symptoms started)",
                    "age_at_diagnosis": "integer or null (age when diagnosed)",
                    "sex": "string or null (male/female/other)",
                    "affected_status": "string (affected/unaffected/uncertain)",
                    "phenotype_details": "string (disease manifestations for this person)",
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
        "challenges": ["any issues during extraction"],
        "notes": "any additional notes about the extraction process"
    }}
}}

IMPORTANT NOTES:
- ONLY include variants in the {gene_symbol} gene. Do NOT include variants from other genes even if they are mentioned in the paper.
- If the paper mentions the target gene {gene_symbol} but does not report any variants, return an empty variants list.
- If full text is not available and only abstract is provided, note this limitation
- Be thorough but accurate - don't invent data not present in the paper
- If a field is not available, use null or an empty string as appropriate
- Preserve exact nomenclature from the paper
"""
