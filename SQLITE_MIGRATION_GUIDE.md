# SQLite Migration Guide

This guide explains how to migrate your Gene Variant Fetcher data from file-based storage to a SQLite database.

## Overview

The migration system consists of three main components:

1. **`migrate_to_sqlite.py`** - Main migration script that creates the database and imports data
2. **`query_variants_db.py`** - Query helper for exploring and analyzing the database
3. **SQLite Database (`variants.db`)** - Normalized relational database for efficient storage and querying

## Benefits of SQLite Migration

- **Reduced disk usage**: Eliminates redundant directories and scattered JSON files
- **Faster queries**: Indexed lookups and SQL-based filtering
- **Data integrity**: Foreign keys and constraints ensure consistency
- **Scalability**: Easily query across thousands of papers and variants
- **Archival**: Original files compressed into ZIP for backup

## Database Schema

The database uses a normalized schema with the following tables:

### Core Tables

- **`papers`** - Paper metadata (PMID, title, journal, etc.)
- **`variants`** - Unique variants with HGVS notations
- **`variant_papers`** - Many-to-many relationship between variants and papers

### Data Tables

- **`penetrance_data`** - Cohort-level penetrance statistics
- **`age_dependent_penetrance`** - Age-stratified penetrance
- **`individual_records`** - Person-level carrier and affected status
- **`functional_data`** - In vitro and functional assays
- **`phenotypes`** - Patient phenotype descriptions
- **`variant_metadata`** - Segregation and population frequency data

### Tracking Tables

- **`extraction_metadata`** - Extraction quality and confidence
- **`tables_processed`** - Tables processed from each paper

## Usage

### Installation

First, ensure you have Python 3.8+ and install dependencies:

```bash
pip install tabulate  # For query_variants_db.py
```

### Basic Migration

To migrate a data directory to SQLite:

```bash
python migrate_to_sqlite.py --data-dir automated_output/TTR/20251125_114028
```

This will:
1. Create `variants.db` in the current directory
2. Import all JSON files from the `extractions` subdirectory
3. Show migration statistics

### Migration with Cleanup

To migrate AND clean up the file system:

```bash
python migrate_to_sqlite.py \
    --data-dir automated_output/TTR/20251125_114028 \
    --cleanup \
    --delete-pmc-after-archive
```

This will:
1. Migrate data to SQLite
2. Delete all empty directories (especially `*_supplements` folders)
3. Archive `pmc_fulltext` to `pmc_fulltext.zip`
4. Delete original `pmc_fulltext` directory after successful archival

### Custom Database Name

```bash
python migrate_to_sqlite.py \
    --data-dir automated_output/TTR/20251125_114028 \
    --db ttr_variants.db
```

### Alternative Extraction Directory

If your extraction files are in a different subdirectory (e.g., `extractions_rerun`):

```bash
python migrate_to_sqlite.py \
    --data-dir automated_output/TTR/20251125_114028 \
    --extractions-subdir extractions_rerun
```

### Dry Run (Preview Only)

To see what would happen without actually migrating:

```bash
python migrate_to_sqlite.py \
    --data-dir automated_output/TTR/20251125_114028 \
    --cleanup \
    --dry-run
```

## Querying the Database

### Database Statistics

```bash
python query_variants_db.py --stats
```

Example output:
```
=== DATABASE STATISTICS ===
Papers: 142
Variants: 387
Unique Genes: 3
Variants With Penetrance: 156
Individuals Affected: 1234
Individuals Unaffected: 567
```

### List All Genes

```bash
python query_variants_db.py --list-genes
```

### List Variants for a Gene

```bash
python query_variants_db.py --list-variants --gene TTR
```

### Get Variant Details

```bash
python query_variants_db.py --variant "p.Val30Met" --details
```

### Penetrance Summary for a Gene

```bash
python query_variants_db.py --gene TTR --penetrance
```

### Search Variants

Search pathogenic variants with at least 10 carriers:

```bash
python query_variants_db.py \
    --search \
    --clinical-significance pathogenic \
    --min-carriers 10
```

### JSON Output

All queries support JSON output for programmatic use:

```bash
python query_variants_db.py --list-genes --json > genes.json
```

## Advanced SQL Queries

You can also query the database directly using SQLite:

```bash
sqlite3 variants.db
```

### Example Queries

**Find variants with highest penetrance:**

```sql
SELECT
    v.protein_notation,
    v.gene_symbol,
    AVG(pd.penetrance_percentage) as avg_penetrance,
    SUM(pd.total_carriers_observed) as total_carriers
FROM variants v
JOIN penetrance_data pd ON v.variant_id = pd.variant_id
GROUP BY v.variant_id
HAVING total_carriers >= 5
ORDER BY avg_penetrance DESC
LIMIT 10;
```

**Count affected vs unaffected by variant:**

```sql
SELECT
    v.protein_notation,
    ir.affected_status,
    COUNT(*) as count
FROM variants v
JOIN individual_records ir ON v.variant_id = ir.variant_id
WHERE v.gene_symbol = 'TTR'
GROUP BY v.variant_id, ir.affected_status;
```

**Papers with most variants:**

```sql
SELECT
    p.pmid,
    p.title,
    COUNT(DISTINCT vp.variant_id) as variant_count
FROM papers p
JOIN variant_papers vp ON p.pmid = vp.pmid
GROUP BY p.pmid
ORDER BY variant_count DESC
LIMIT 10;
```

## Migration Performance

- **Small datasets** (< 100 papers): ~1-2 minutes
- **Medium datasets** (100-500 papers): ~5-10 minutes
- **Large datasets** (500+ papers): ~15-30 minutes

Migration is I/O bound, so SSD storage provides significant speedup.

## Data Integrity

The migration script includes validation:

- **Unique constraints**: Prevents duplicate variants
- **Foreign keys**: Ensures referential integrity
- **NULL handling**: Flexible matching for incomplete data
- **Transaction safety**: Rollback on errors

## Backup Recommendations

Before migration:

1. **Backup original data**: `cp -r automated_output automated_output.backup`
2. **Use dry-run first**: Test migration without changes
3. **Archive before cleanup**: Ensure ZIP creation succeeds before deleting originals

After migration:

1. **Backup database**: `cp variants.db variants.db.backup`
2. **Keep archived files**: Don't delete `pmc_fulltext.zip`
3. **Version control**: Consider adding to git (if < 100MB)

## Troubleshooting

### "No extraction directory found"

The script searches for `extractions`, `extractions_rerun`, and `extraction` directories. If your directory has a different name:

```bash
python migrate_to_sqlite.py \
    --data-dir /path/to/data \
    --extractions-subdir your_directory_name
```

### "Failed to migrate file"

Check the error log. Common issues:
- Malformed JSON files
- Missing required fields
- Encoding issues

The migration continues even if individual files fail. Review the error summary at the end.

### "Database is locked"

Close any other connections to the database:

```bash
# Find processes using the database
lsof variants.db

# Or just remove lock file
rm variants.db-wal variants.db-shm
```

### Performance Issues

For very large datasets:

1. Use SSD storage
2. Increase SQLite cache: `PRAGMA cache_size = 10000;`
3. Disable synchronous mode during migration: `PRAGMA synchronous = OFF;`

## Integration with Existing Code

You can query the database from your Python code:

```python
import sqlite3

conn = sqlite3.connect('variants.db')
cursor = conn.cursor()

# Get all variants for a gene
cursor.execute("""
    SELECT protein_notation, clinical_significance
    FROM variants
    WHERE gene_symbol = ?
""", ("TTR",))

for row in cursor.fetchall():
    print(f"Variant: {row[0]}, Significance: {row[1]}")

conn.close()
```

Or use the query helper as a library:

```python
from query_variants_db import VariantDatabaseQuery

with VariantDatabaseQuery("variants.db") as db:
    variants = db.list_variants(gene_symbol="TTR")
    for v in variants:
        print(v['protein_notation'])
```

## Future Enhancements

Potential improvements:

- **Web interface**: Flask/Django app for browsing
- **Export functions**: Export to CSV, Excel, or JSON
- **Incremental updates**: Add new extractions without re-migration
- **Multi-database**: Separate databases per gene
- **Full-text search**: SQLite FTS5 for text search

## Support

For issues or questions:

1. Check the error logs in the migration output
2. Review this guide's troubleshooting section
3. Open an issue on GitHub

## Schema Diagram

```
papers (pmid, title, ...)
  ↓
  ├─ variant_papers (variant_id, pmid)
  │    ↓
  │  variants (variant_id, gene_symbol, protein_notation, ...)
  │    ↓
  │    ├─ penetrance_data (variant_id, pmid, ...)
  │    │    ↓
  │    │  age_dependent_penetrance
  │    ├─ individual_records (variant_id, pmid, ...)
  │    ├─ functional_data (variant_id, pmid, ...)
  │    ├─ phenotypes (variant_id, pmid, ...)
  │    └─ variant_metadata (variant_id, pmid, ...)
  │
  ├─ extraction_metadata (pmid, ...)
  └─ tables_processed (pmid, table_name, ...)
```

## License

This migration system is part of the Gene Variant Fetcher project and shares the same license.
