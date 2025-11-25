#!/usr/bin/env python3
"""
Query Helper for Gene Variant SQLite Database

Provides convenient functions to query the variants database.

Usage examples:
    python query_variants_db.py --list-genes
    python query_variants_db.py --gene TTR --list-variants
    python query_variants_db.py --variant "p.Val30Met" --details
    python query_variants_db.py --stats
"""

import sqlite3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Any, Optional
from tabulate import tabulate


class VariantDatabaseQuery:
    """Query interface for variant database."""

    def __init__(self, db_path: str):
        """
        Initialize database connection.

        Args:
            db_path: Path to SQLite database
        """
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Enable column access by name

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()

    def get_database_stats(self) -> Dict[str, int]:
        """Get overall database statistics."""
        cursor = self.conn.cursor()

        stats = {}

        # Count records in each table
        tables = [
            'papers', 'variants', 'variant_papers', 'penetrance_data',
            'individual_records', 'functional_data', 'phenotypes'
        ]

        for table in tables:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            stats[table] = cursor.fetchone()[0]

        # Count unique genes
        cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM variants")
        stats['unique_genes'] = cursor.fetchone()[0]

        # Count variants with penetrance data
        cursor.execute("""
            SELECT COUNT(DISTINCT variant_id)
            FROM penetrance_data
        """)
        stats['variants_with_penetrance'] = cursor.fetchone()[0]

        # Count affected vs unaffected individuals
        cursor.execute("""
            SELECT affected_status, COUNT(*)
            FROM individual_records
            GROUP BY affected_status
        """)
        for row in cursor.fetchall():
            stats[f'individuals_{row[0]}'] = row[1]

        return stats

    def list_genes(self) -> List[Dict[str, Any]]:
        """List all genes with variant counts."""
        cursor = self.conn.cursor()

        cursor.execute("""
            SELECT
                gene_symbol,
                COUNT(DISTINCT variant_id) as variant_count,
                COUNT(DISTINCT pmid) as paper_count
            FROM variants v
            LEFT JOIN variant_papers vp ON v.variant_id = vp.variant_id
            GROUP BY gene_symbol
            ORDER BY variant_count DESC
        """)

        return [dict(row) for row in cursor.fetchall()]

    def list_variants(self, gene_symbol: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        List variants, optionally filtered by gene.

        Args:
            gene_symbol: Filter by gene symbol (optional)

        Returns:
            List of variant dictionaries
        """
        cursor = self.conn.cursor()

        query = """
            SELECT
                v.variant_id,
                v.gene_symbol,
                v.protein_notation,
                v.cdna_notation,
                v.clinical_significance,
                COUNT(DISTINCT vp.pmid) as paper_count,
                COUNT(DISTINCT ir.record_id) as individual_count
            FROM variants v
            LEFT JOIN variant_papers vp ON v.variant_id = vp.variant_id
            LEFT JOIN individual_records ir ON v.variant_id = ir.variant_id
        """

        if gene_symbol:
            query += " WHERE v.gene_symbol = ?"
            cursor.execute(query + " GROUP BY v.variant_id ORDER BY paper_count DESC", (gene_symbol,))
        else:
            cursor.execute(query + " GROUP BY v.variant_id ORDER BY paper_count DESC")

        return [dict(row) for row in cursor.fetchall()]

    def get_variant_details(
        self,
        variant_id: Optional[int] = None,
        protein_notation: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a variant.

        Args:
            variant_id: Variant ID
            protein_notation: Protein notation (alternative to variant_id)

        Returns:
            Detailed variant dictionary
        """
        cursor = self.conn.cursor()

        # Get basic variant info
        if variant_id:
            cursor.execute("SELECT * FROM variants WHERE variant_id = ?", (variant_id,))
        elif protein_notation:
            cursor.execute("SELECT * FROM variants WHERE protein_notation = ?", (protein_notation,))
        else:
            return None

        variant_row = cursor.fetchone()
        if not variant_row:
            return None

        variant = dict(variant_row)
        vid = variant['variant_id']

        # Get associated papers
        cursor.execute("""
            SELECT vp.pmid, p.title, vp.source_location, vp.key_quotes
            FROM variant_papers vp
            LEFT JOIN papers p ON vp.pmid = p.pmid
            WHERE vp.variant_id = ?
        """, (vid,))
        variant['papers'] = [dict(row) for row in cursor.fetchall()]

        # Get penetrance data
        cursor.execute("""
            SELECT *
            FROM penetrance_data
            WHERE variant_id = ?
        """, (vid,))
        variant['penetrance_data'] = [dict(row) for row in cursor.fetchall()]

        # Get individual records
        cursor.execute("""
            SELECT *
            FROM individual_records
            WHERE variant_id = ?
        """, (vid,))
        variant['individual_records'] = [dict(row) for row in cursor.fetchall()]

        # Get functional data
        cursor.execute("""
            SELECT *
            FROM functional_data
            WHERE variant_id = ?
        """, (vid,))
        variant['functional_data'] = [dict(row) for row in cursor.fetchall()]

        # Get phenotypes
        cursor.execute("""
            SELECT *
            FROM phenotypes
            WHERE variant_id = ?
        """, (vid,))
        variant['phenotypes'] = [dict(row) for row in cursor.fetchall()]

        return variant

    def get_aggregated_penetrance(self, gene_symbol: str) -> List[Dict[str, Any]]:
        """
        Get aggregated penetrance statistics for a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            List of variant penetrance summaries
        """
        cursor = self.conn.cursor()

        cursor.execute("""
            SELECT
                v.protein_notation,
                v.cdna_notation,
                v.clinical_significance,
                SUM(pd.total_carriers_observed) as total_carriers,
                SUM(pd.affected_count) as total_affected,
                SUM(pd.unaffected_count) as total_unaffected,
                AVG(pd.penetrance_percentage) as avg_penetrance,
                COUNT(DISTINCT vp.pmid) as paper_count,
                COUNT(DISTINCT ir.record_id) as individual_count
            FROM variants v
            LEFT JOIN penetrance_data pd ON v.variant_id = pd.variant_id
            LEFT JOIN variant_papers vp ON v.variant_id = vp.variant_id
            LEFT JOIN individual_records ir ON v.variant_id = ir.variant_id
            WHERE v.gene_symbol = ?
            GROUP BY v.variant_id
            ORDER BY total_carriers DESC
        """, (gene_symbol,))

        return [dict(row) for row in cursor.fetchall()]

    def search_variants(
        self,
        gene: Optional[str] = None,
        clinical_significance: Optional[str] = None,
        min_carriers: int = 0
    ) -> List[Dict[str, Any]]:
        """
        Search variants with filters.

        Args:
            gene: Gene symbol filter
            clinical_significance: Clinical significance filter
            min_carriers: Minimum number of carriers

        Returns:
            List of matching variants
        """
        cursor = self.conn.cursor()

        query = """
            SELECT
                v.variant_id,
                v.gene_symbol,
                v.protein_notation,
                v.cdna_notation,
                v.clinical_significance,
                COALESCE(SUM(pd.total_carriers_observed), 0) as total_carriers,
                COUNT(DISTINCT vp.pmid) as paper_count
            FROM variants v
            LEFT JOIN penetrance_data pd ON v.variant_id = pd.variant_id
            LEFT JOIN variant_papers vp ON v.variant_id = vp.variant_id
            WHERE 1=1
        """

        params = []

        if gene:
            query += " AND v.gene_symbol = ?"
            params.append(gene)

        if clinical_significance:
            query += " AND v.clinical_significance LIKE ?"
            params.append(f"%{clinical_significance}%")

        query += " GROUP BY v.variant_id"

        if min_carriers > 0:
            query += " HAVING total_carriers >= ?"
            params.append(min_carriers)

        query += " ORDER BY total_carriers DESC"

        cursor.execute(query, params)

        return [dict(row) for row in cursor.fetchall()]


def main():
    """CLI entrypoint."""
    parser = argparse.ArgumentParser(
        description="Query Gene Variant SQLite Database"
    )

    parser.add_argument(
        "--db",
        type=str,
        default="variants.db",
        help="Path to SQLite database (default: variants.db)"
    )

    # Query modes
    parser.add_argument("--stats", action="store_true", help="Show database statistics")
    parser.add_argument("--list-genes", action="store_true", help="List all genes")
    parser.add_argument("--list-variants", action="store_true", help="List variants")
    parser.add_argument("--gene", type=str, help="Filter by gene symbol")
    parser.add_argument("--variant", type=str, help="Get variant details (protein notation)")
    parser.add_argument("--details", action="store_true", help="Show detailed information")
    parser.add_argument("--penetrance", action="store_true", help="Show penetrance summary")
    parser.add_argument("--search", action="store_true", help="Search variants")
    parser.add_argument("--clinical-significance", type=str, help="Filter by clinical significance")
    parser.add_argument("--min-carriers", type=int, default=0, help="Minimum carriers filter")

    # Output format
    parser.add_argument("--json", action="store_true", help="Output as JSON")

    args = parser.parse_args()

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found: {db_path}")
        return 1

    with VariantDatabaseQuery(str(db_path)) as db:

        # Database statistics
        if args.stats:
            stats = db.get_database_stats()
            if args.json:
                print(json.dumps(stats, indent=2))
            else:
                print("\n=== DATABASE STATISTICS ===")
                for key, value in stats.items():
                    print(f"{key.replace('_', ' ').title()}: {value}")

        # List genes
        elif args.list_genes:
            genes = db.list_genes()
            if args.json:
                print(json.dumps(genes, indent=2))
            else:
                print("\n=== GENES ===")
                print(tabulate(genes, headers="keys", tablefmt="grid"))

        # List variants
        elif args.list_variants:
            variants = db.list_variants(gene_symbol=args.gene)
            if args.json:
                print(json.dumps(variants, indent=2))
            else:
                print(f"\n=== VARIANTS{' FOR ' + args.gene if args.gene else ''} ===")
                print(tabulate(variants, headers="keys", tablefmt="grid"))

        # Variant details
        elif args.variant:
            variant = db.get_variant_details(protein_notation=args.variant)
            if args.json:
                print(json.dumps(variant, indent=2, default=str))
            else:
                if variant:
                    print(f"\n=== VARIANT DETAILS: {args.variant} ===")
                    print(f"\nGene: {variant['gene_symbol']}")
                    print(f"Protein: {variant['protein_notation']}")
                    print(f"cDNA: {variant['cdna_notation']}")
                    print(f"Clinical Significance: {variant['clinical_significance']}")
                    print(f"\nFound in {len(variant['papers'])} papers")
                    print(f"Individual records: {len(variant['individual_records'])}")

                    if args.details:
                        print("\n--- Papers ---")
                        for paper in variant['papers']:
                            print(f"  PMID {paper['pmid']}: {paper['title']}")

                        print("\n--- Penetrance Data ---")
                        for pd in variant['penetrance_data']:
                            print(f"  PMID {pd['pmid']}:")
                            print(f"    Total carriers: {pd['total_carriers_observed']}")
                            print(f"    Affected: {pd['affected_count']}")
                            print(f"    Unaffected: {pd['unaffected_count']}")
                            print(f"    Penetrance: {pd['penetrance_percentage']}%")

                        if variant['individual_records']:
                            print("\n--- Individual Records ---")
                            print(tabulate(
                                variant['individual_records'],
                                headers="keys",
                                tablefmt="grid"
                            ))
                else:
                    print(f"Variant not found: {args.variant}")

        # Penetrance summary
        elif args.penetrance and args.gene:
            penetrance = db.get_aggregated_penetrance(args.gene)
            if args.json:
                print(json.dumps(penetrance, indent=2, default=str))
            else:
                print(f"\n=== PENETRANCE SUMMARY FOR {args.gene} ===")
                print(tabulate(penetrance, headers="keys", tablefmt="grid"))

        # Search variants
        elif args.search:
            results = db.search_variants(
                gene=args.gene,
                clinical_significance=args.clinical_significance,
                min_carriers=args.min_carriers
            )
            if args.json:
                print(json.dumps(results, indent=2))
            else:
                print("\n=== SEARCH RESULTS ===")
                print(tabulate(results, headers="keys", tablefmt="grid"))

        else:
            parser.print_help()

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
