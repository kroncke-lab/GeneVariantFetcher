#!/usr/bin/env python3
"""Run SQLite migration for KCNH2 extractions."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher')

from pathlib import Path
from harvesting.migrate_to_sqlite import create_database_schema, migrate_extraction_directory

extraction_dir = Path('/mnt/temp2/kronckbm/gvf_output/KCNH2/20260202_173749/extractions')
db_path = '/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db'

print(f"Creating database at: {db_path}")
conn = create_database_schema(db_path)

print(f"Migrating from: {extraction_dir}")
stats = migrate_extraction_directory(conn, extraction_dir)

print(f"\nMigration complete!")
print(f"Total files: {stats['total_files']}")
print(f"Successful: {stats['successful']}")
print(f"Failed: {stats['failed']}")

# Get counts
cursor = conn.cursor()
cursor.execute("SELECT COUNT(*) FROM papers")
print(f"Papers: {cursor.fetchone()[0]}")

cursor.execute("SELECT COUNT(*) FROM variants")
print(f"Variants: {cursor.fetchone()[0]}")

cursor.execute("SELECT COUNT(*) FROM individual_records")
print(f"Individual records: {cursor.fetchone()[0]}")

conn.close()
print("Done!")
