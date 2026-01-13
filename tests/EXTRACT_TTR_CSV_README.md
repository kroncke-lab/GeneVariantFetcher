# TTR Database CSV Extraction Guide

## Overview

The `extract_ttr_to_csv.py` script extracts variant data from the TTR.db SQLite database and exports it to CSV format.

## What Gets Extracted

### Main Variant CSV (`ttr_variants.csv`)
- **Variant Identity**: cDNA notation, protein notation, genomic position
- **Number of Affected**: Aggregated count from penetrance_data and individual_records
- **Number of Unaffected**: Aggregated count from penetrance_data and individual_records
- **Number of Uncertain**: Count of uncertain cases
- **Notes**: Additional notes from variant_papers table (concatenated across all papers)
- **Clinical Significance**: Pathogenic, benign, VUS, etc.
- **Source Information**: PMIDs and source locations where variant was reported
- **Phenotype Descriptions**: Disease phenotype information

### Individual Records CSV (`ttr_variants_individual_records.csv`)
- **Individual-level data** for each person carrying a variant
- **Onset Information**: Age at onset, age at diagnosis, age at evaluation
- **Affected Status**: affected, unaffected, or uncertain
- **Phenotype Details**: Individual-specific disease manifestations
- **Evidence Sentences**: Exact quotes from papers

## Usage

### Basic Usage

```bash
python3 extract_ttr_to_csv.py --db "/Users/kronckbm/OneDrive-VUMC/Kroncke_Lab/MAGen/hATTR/TTR/TTR.db"
```

This will create:
- `ttr_variants.csv` - Main variant data
- `ttr_variants_individual_records.csv` - Individual-level records with onset data

### Custom Output File

```bash
python3 extract_ttr_to_csv.py \
  --db "/Users/kronckbm/OneDrive-VUMC/Kroncke_Lab/MAGen/hATTR/TTR/TTR.db" \
  --output my_variants.csv
```

This will create:
- `my_variants.csv` - Main variant data
- `my_variants_individual_records.csv` - Individual records

### Skip Individual Records

If you only want the aggregated variant data:

```bash
python3 extract_ttr_to_csv.py \
  --db "/Users/kronckbm/OneDrive-VUMC/Kroncke_Lab/MAGen/hATTR/TTR/TTR.db" \
  --no-individual-records
```

## CSV Column Descriptions

### Main Variant CSV Columns

| Column | Description |
|--------|-------------|
| `variant_id` | Unique variant identifier in database |
| `gene_symbol` | Gene symbol (TTR) |
| `cdna_notation` | cDNA notation (e.g., "c.148G>A") |
| `protein_notation` | Protein notation (e.g., "p.Val30Met") |
| `genomic_position` | Genomic position if available |
| `clinical_significance` | Pathogenic, benign, VUS, etc. |
| `evidence_level` | Evidence level for classification |
| `total_carriers` | Total carriers observed (sum across all papers) |
| `total_affected` | Total affected count (sum across all papers) |
| `total_unaffected` | Total unaffected count (sum across all papers) |
| `total_uncertain` | Total uncertain count (sum across all papers) |
| `individual_affected_count` | Count from individual_records table |
| `individual_unaffected_count` | Count from individual_records table |
| `individual_uncertain_count` | Count from individual_records table |
| `notes` | Additional notes (concatenated with " \| " separator) |
| `source_pmids` | Comma-separated list of PMIDs reporting this variant |
| `source_locations` | Source locations (e.g., "Table 2") concatenated |
| `phenotype_descriptions` | Phenotype descriptions concatenated |
| `paper_count` | Number of papers reporting this variant |

### Individual Records CSV Columns

| Column | Description |
|--------|-------------|
| `record_id` | Unique record identifier |
| `variant_id` | Links to variant in main CSV |
| `cdna_notation` | cDNA notation for this variant |
| `protein_notation` | Protein notation for this variant |
| `individual_id` | Individual identifier (e.g., "II-1", "P1", "Case_2") |
| `age_at_evaluation` | Age when evaluated/assessed |
| `age_at_onset` | **Age at disease onset** (when available) |
| `age_at_diagnosis` | Age at diagnosis |
| `sex` | Male, female, or other |
| `affected_status` | affected, unaffected, or uncertain |
| `phenotype_details` | Individual-specific phenotype details |
| `evidence_sentence` | Exact sentence from paper |
| `pmid` | PubMed ID of source paper |
| `paper_title` | Title of source paper |

## Notes

- The script aggregates data across all papers for each variant
- If a variant appears in multiple papers, counts are summed
- Notes and source locations are concatenated with " \| " separator
- Individual records are sorted by variant, then by affected status, then by age at onset
- Empty fields are represented as empty strings in the CSV

## Troubleshooting

### Database Not Found
```
Error: Database file not found: /path/to/TTR.db
```
Make sure the path to TTR.db is correct. Use absolute paths if needed.

### No Variants Found
```
Warning: No variants found in database.
```
Check that the database contains TTR variants. The script filters for gene_symbol = 'TTR'.

