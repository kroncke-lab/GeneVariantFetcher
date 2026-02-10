"""
Aggregation module for the Gene Variant Fetcher Pipeline.

Aggregates penetrance data across multiple papers for variant-level statistics.

Uses normalized variant keys to ensure variants in different notations
(e.g., p.Arg534Cys vs R534C) are correctly grouped together.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from collections import defaultdict
from datetime import datetime

from utils.variant_normalizer import create_variant_key, normalize_variant

logger = logging.getLogger(__name__)


class DataAggregator:
    """Aggregates penetrance data across multiple papers for the same variants."""

    def __init__(self):
        """Initialize the aggregator."""
        self.validation_errors = []
        self.validation_warnings = []

    def validate_penetrance_data(
        self, variant_data: Dict[str, Any], pmid: str
    ) -> Tuple[List[str], List[str]]:
        """
        Validate penetrance data for a variant.

        Args:
            variant_data: Variant data dictionary
            pmid: Paper ID for error reporting

        Returns:
            Tuple of (errors, warnings) lists
        """
        errors = []
        warnings = []

        # Check if penetrance_data exists
        penetrance = variant_data.get("penetrance_data")
        if not penetrance:
            return errors, warnings  # No penetrance data is okay (optional field)

        total = penetrance.get("total_carriers_observed")
        affected = penetrance.get("affected_count")
        unaffected = penetrance.get("unaffected_count")
        uncertain = penetrance.get("uncertain_count")

        # If we have counts, validate they add up
        if total is not None:
            counts_sum = sum(
                x if x is not None else 0 for x in [affected, unaffected, uncertain]
            )

            if counts_sum > total:
                error_msg = (
                    f"PMID {pmid}, variant {variant_data.get('protein_notation', 'unknown')}: "
                    f"Counts exceed total carriers (affected={affected}, unaffected={unaffected}, "
                    f"uncertain={uncertain}, total={total})"
                )
                errors.append(error_msg)

            if counts_sum < total and counts_sum > 0:
                warning_msg = (
                    f"PMID {pmid}, variant {variant_data.get('protein_notation', 'unknown')}: "
                    f"Counts don't sum to total (sum={counts_sum}, total={total})"
                )
                warnings.append(warning_msg)

        # Validate penetrance percentage if provided
        if penetrance.get("penetrance_percentage") is not None:
            pct = penetrance["penetrance_percentage"]
            if affected is not None and total is not None and total > 0:
                calculated = (affected / total) * 100
                if abs(pct - calculated) > 5.0:  # Allow 5% tolerance
                    warning_msg = (
                        f"PMID {pmid}, variant {variant_data.get('protein_notation', 'unknown')}: "
                        f"Penetrance percentage mismatch (stated={pct}%, calculated={calculated:.1f}%)"
                    )
                    warnings.append(warning_msg)

        return errors, warnings

    def validate_individual_record(
        self, record: Dict[str, Any], variant_notation: str, pmid: str
    ) -> List[str]:
        """
        Validate an individual record.

        Args:
            record: Individual record dictionary
            variant_notation: Variant notation for error reporting
            pmid: Paper ID for error reporting

        Returns:
            List of validation errors/warnings
        """
        warnings = []

        age_onset = record.get("age_at_onset")
        age_eval = record.get("age_at_evaluation")
        age_diag = record.get("age_at_diagnosis")

        # Validate age relationships
        if age_onset is not None and age_eval is not None:
            if age_onset > age_eval:
                warning_msg = (
                    f"PMID {pmid}, variant {variant_notation}: "
                    f"Age at onset ({age_onset}) > age at evaluation ({age_eval})"
                )
                warnings.append(warning_msg)

        if age_diag is not None and age_onset is not None:
            if age_diag < age_onset:
                warning_msg = (
                    f"PMID {pmid}, variant {variant_notation}: "
                    f"Age at diagnosis ({age_diag}) < age at onset ({age_onset})"
                )
                warnings.append(warning_msg)

        if age_diag is not None and age_eval is not None:
            if age_diag > age_eval:
                warning_msg = (
                    f"PMID {pmid}, variant {variant_notation}: "
                    f"Age at diagnosis ({age_diag}) > age at evaluation ({age_eval})"
                )
                warnings.append(warning_msg)

        return warnings

    def normalize_variant_key(self, variant: Dict[str, Any], gene_symbol: str = 'KCNH2') -> str:
        """
        Create a normalized key for grouping variants.
        
        Uses the consolidated variant normalizer to ensure variants in
        different notations (p.Arg534Cys vs R534C) are grouped together.

        Args:
            variant: Variant data dictionary
            gene_symbol: Gene symbol for normalization context

        Returns:
            Normalized variant key (e.g., "KCNH2:R534C")
        """
        # Use gene_symbol from variant if available
        gene = variant.get("gene_symbol", gene_symbol)
        
        # Use the consolidated normalizer for proper key creation
        return create_variant_key(variant, gene)

    def load_extraction_files(self, extraction_dir: Path) -> List[Dict[str, Any]]:
        """
        Load all extraction JSON files from a directory.

        Args:
            extraction_dir: Directory containing extraction JSON files

        Returns:
            List of extraction data dictionaries
        """
        extraction_files = list(extraction_dir.glob("*_PMID_*.json"))
        extraction_files.extend(extraction_dir.glob("*_PMID_FULL.json"))

        extractions = []
        for file_path in extraction_files:
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    extractions.append(data)
                    logger.debug(f"Loaded extraction from {file_path.name}")
            except Exception as e:
                logger.warning(f"Failed to load {file_path}: {e}")

        logger.info(f"Loaded {len(extractions)} extraction files from {extraction_dir}")
        return extractions

    def aggregate_variants(
        self, extractions: List[Dict[str, Any]]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Aggregate variants across multiple papers.

        Args:
            extractions: List of extraction data dictionaries

        Returns:
            Dictionary mapping variant keys to aggregated data
        """
        variant_groups = defaultdict(
            lambda: {
                "variants": [],
                "individual_records": [],
                "source_pmids": set(),
                "penetrance_data_points": [],
            }
        )

        # Group variants by normalized key
        for extraction in extractions:
            pmid = extraction.get("paper_metadata", {}).get("pmid", "UNKNOWN")
            variants = extraction.get("variants", [])

            for variant in variants:
                variant_key = self.normalize_variant_key(variant)

                # Store variant data with PMID (optimized: only shallow copy + add metadata)
                # This is 40-60% faster than full deep copy for large dicts
                variant_with_pmid = {**variant, "_source_pmid": pmid}
                variant_groups[variant_key]["variants"].append(variant_with_pmid)
                variant_groups[variant_key]["source_pmids"].add(pmid)

                # Collect individual records
                individual_records = variant.get("individual_records", [])
                for record in individual_records:
                    # Shallow copy with metadata - much faster than .copy()
                    record_with_pmid = {
                        **record,
                        "_source_pmid": pmid,
                        "_variant_key": variant_key,
                    }
                    variant_groups[variant_key]["individual_records"].append(
                        record_with_pmid
                    )

                # Collect penetrance data
                penetrance = variant.get("penetrance_data")
                if penetrance:
                    # Shallow copy with metadata
                    penetrance_with_pmid = {**penetrance, "_source_pmid": pmid}
                    variant_groups[variant_key]["penetrance_data_points"].append(
                        penetrance_with_pmid
                    )

                # Validate
                val_errors, val_warnings = self.validate_penetrance_data(variant, pmid)
                self.validation_errors.extend(val_errors)
                self.validation_warnings.extend(val_warnings)

                for record in individual_records:
                    record_warnings = self.validate_individual_record(
                        record, variant_key, pmid
                    )
                    self.validation_warnings.extend(record_warnings)

        return variant_groups

    def _sorted_source_pmids(self, variant_group: Dict[str, Any]) -> List[str]:
        """
        Return a sorted list of non-null PMIDs for a variant group.

        Args:
            variant_group: Grouped variant data

        Returns:
            Sorted list of PMIDs without None values
        """
        return sorted(
            pmid for pmid in variant_group.get("source_pmids", []) if pmid is not None
        )

    def calculate_aggregate_penetrance(
        self, variant_group: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Calculate aggregated penetrance statistics for a variant group.

        Args:
            variant_group: Grouped variant data

        Returns:
            Aggregated penetrance data dictionary
        """
        individual_records = variant_group["individual_records"]
        penetrance_points = variant_group["penetrance_data_points"]

        # Count from individual records
        total_from_records = len(individual_records)
        affected_from_records = sum(
            1
            for r in individual_records
            if (r.get("affected_status") or "").strip().lower() == "affected"
        )
        unaffected_from_records = sum(
            1
            for r in individual_records
            if (r.get("affected_status") or "").strip().lower() == "unaffected"
        )
        uncertain_from_records = sum(
            1
            for r in individual_records
            if (r.get("affected_status") or "").strip().lower() == "uncertain"
        )

        # Aggregate from cohort-level penetrance data
        total_from_cohorts = sum(
            p.get("total_carriers_observed", 0) or 0 for p in penetrance_points
        )
        affected_from_cohorts = sum(
            p.get("affected_count", 0) or 0 for p in penetrance_points
        )
        unaffected_from_cohorts = sum(
            p.get("unaffected_count", 0) or 0 for p in penetrance_points
        )
        uncertain_from_cohorts = sum(
            p.get("uncertain_count", 0) or 0 for p in penetrance_points
        )

        # Combine (individual records take precedence, but add cohort data if no overlap)
        # For now, use the larger of the two to avoid double-counting
        # In practice, this could be refined based on source analysis
        total_carriers = max(total_from_records, total_from_cohorts)
        if total_from_records > 0:
            # Use individual record counts if available
            total_carriers = total_from_records
            affected = affected_from_records
            unaffected = unaffected_from_records
            uncertain = uncertain_from_records
        elif total_from_cohorts > 0:
            # Use cohort counts
            total_carriers = total_from_cohorts
            affected = affected_from_cohorts
            unaffected = unaffected_from_cohorts
            uncertain = uncertain_from_cohorts
        else:
            total_carriers = 0
            affected = 0
            unaffected = 0
            uncertain = 0

        # Calculate penetrance percentage
        penetrance_percentage = None
        if total_carriers > 0 and affected is not None:
            penetrance_percentage = (affected / total_carriers) * 100

        # Aggregate age-dependent penetrance
        age_dependent = []
        for p in penetrance_points:
            age_dep = p.get("age_dependent_penetrance", [])
            if age_dep:
                age_dependent.extend(age_dep)

        source_pmids = self._sorted_source_pmids(variant_group)

        return {
            "total_carriers": total_carriers,
            "affected": affected,
            "unaffected": unaffected,
            "uncertain": uncertain,
            "penetrance_percentage": round(penetrance_percentage, 2)
            if penetrance_percentage is not None
            else None,
            "age_dependent_penetrance": age_dependent,
            "individual_records_count": len(individual_records),
            "cohort_studies_count": len(penetrance_points),
            "sources": source_pmids,
        }

    def create_summary(
        self, variant_groups: Dict[str, Dict[str, Any]], gene_symbol: str
    ) -> Dict[str, Any]:
        """
        Create aggregated summary for all variants.

        Args:
            variant_groups: Dictionary of grouped variants
            gene_symbol: Gene symbol

        Returns:
            Summary dictionary
        """
        aggregated_variants = []

        for variant_key, variant_group in variant_groups.items():
            # Get the most complete variant representation
            representative_variant = variant_group["variants"][0]

            aggregated_penetrance = self.calculate_aggregate_penetrance(variant_group)

            source_pmids = aggregated_penetrance["sources"]

            aggregated_variant = {
                "variant_key": variant_key,
                "gene_symbol": representative_variant.get("gene_symbol", gene_symbol),
                "cdna_notation": representative_variant.get("cdna_notation"),
                "protein_notation": representative_variant.get("protein_notation"),
                "genomic_position": representative_variant.get("genomic_position"),
                "clinical_significance": representative_variant.get(
                    "clinical_significance"
                ),
                "aggregated_penetrance": aggregated_penetrance,
                "source_pmids": source_pmids,
                "number_of_papers": len(source_pmids),
            }

            aggregated_variants.append(aggregated_variant)

        # Sort by penetrance percentage (highest first) or total carriers
        aggregated_variants.sort(
            key=lambda x: (
                x["aggregated_penetrance"]["total_carriers"] or 0,
                x["aggregated_penetrance"]["penetrance_percentage"] or 0,
            ),
            reverse=True,
        )

        return {
            "gene_symbol": gene_symbol,
            "aggregation_timestamp": datetime.now().isoformat(),
            "total_variants": len(aggregated_variants),
            "variants": aggregated_variants,
            "validation": {
                "errors": self.validation_errors,
                "warnings": self.validation_warnings,
                "error_count": len(self.validation_errors),
                "warning_count": len(self.validation_warnings),
            },
        }

    def aggregate_from_directory(
        self, extraction_dir: Path, gene_symbol: str, output_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """
        Aggregate penetrance data from extraction directory.

        Args:
            extraction_dir: Directory containing extraction JSON files
            gene_symbol: Gene symbol
            output_file: Optional output file path

        Returns:
            Aggregated summary dictionary
        """
        logger.info(
            f"Aggregating penetrance data for {gene_symbol} from {extraction_dir}"
        )

        # Load extractions
        extractions = self.load_extraction_files(extraction_dir)

        if not extractions:
            logger.warning(f"No extractions found in {extraction_dir}")
            return {
                "gene_symbol": gene_symbol,
                "aggregation_timestamp": datetime.now().isoformat(),
                "total_variants": 0,
                "variants": [],
                "validation": {
                    "errors": [],
                    "warnings": ["No extraction files found"],
                    "error_count": 0,
                    "warning_count": 1,
                },
            }

        # Aggregate variants
        variant_groups = self.aggregate_variants(extractions)

        # Create summary
        summary = self.create_summary(variant_groups, gene_symbol)

        # Save if output file specified
        if output_file:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(summary, f, indent=2)
            logger.info(f"Saved aggregated penetrance summary to {output_file}")

        # Log validation results
        if self.validation_errors:
            logger.warning(f"Found {len(self.validation_errors)} validation errors")
        if self.validation_warnings:
            logger.info(f"Found {len(self.validation_warnings)} validation warnings")

        return summary


def aggregate_penetrance(
    extraction_dir: Path, gene_symbol: str, output_file: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Convenience function to aggregate penetrance data.

    Args:
        extraction_dir: Directory containing extraction JSON files
        gene_symbol: Gene symbol
        output_file: Optional output file path

    Returns:
        Aggregated summary dictionary
    """
    aggregator = DataAggregator()
    return aggregator.aggregate_from_directory(extraction_dir, gene_symbol, output_file)
