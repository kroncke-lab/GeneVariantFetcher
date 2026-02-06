"""
Run Manifest Generator for GeneVariantFetcher

Creates and manages run_manifest.json files that provide complete metadata
about the execution of a gene variant extraction workflow.
"""

import json
import logging
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


class RunManifest:
    """Manages the run manifest for workflow execution tracking."""
    
    def __init__(self, output_dir: Path, gene_symbol: str):
        """
        Initialize a new run manifest.
        
        Args:
            output_dir: The output directory for this run
            gene_symbol: The gene being processed
        """
        self.output_dir = Path(output_dir)
        self.gene_symbol = gene_symbol
        self.run_id = str(uuid.uuid4())
        self.manifest_file = self.output_dir / "run_manifest.json"
        
        # Initialize basic manifest structure
        self.data = {
            "run_id": self.run_id,
            "gene_symbol": gene_symbol,
            "start_timestamp": datetime.now().isoformat(),
            "end_timestamp": None,
            "duration_seconds": None,
            "status": "running",
            "config": {},
            "statistics": {
                "pmids_discovered": 0,
                "pmids_filtered_out": 0,
                "pmids_passed_filters": 0,
                "papers_downloaded": 0,
                "papers_download_failed": 0,
                "papers_extracted": 0,
                "papers_from_fulltext": 0,
                "papers_from_abstract_only": 0,
                "papers_extraction_failed": 0,
                "total_variants_found": 0,
                "variants_with_penetrance_data": 0,
                "total_carriers_observed": 0,
                "total_affected_carriers": 0,
                "success_rate": "0%",
            },
            "output_locations": {},
            "errors": [],
            "warnings": [],
            "cli_arguments": {},
            "environment": {
                "python_version": "",
                "gvf_version": "",
                "hostname": "",
                "working_directory": "",
            }
        }
        
    def set_config(self, **kwargs) -> None:
        """Set configuration parameters from CLI arguments."""
        self.data["config"].update(kwargs)
        self.data["cli_arguments"] = kwargs.copy()
        
    def update_statistics(self, **kwargs) -> None:
        """Update statistics with final workflow metrics."""
        self.data["statistics"].update(kwargs)
        
    def update_output_locations(self, **kwargs) -> None:
        """Update output file locations."""
        self.data["output_locations"].update(kwargs)
        
    def add_error(self, error: str) -> None:
        """Add an error to the manifest."""
        self.data["errors"].append({
            "timestamp": datetime.now().isoformat(),
            "message": str(error)
        })
        
    def add_warning(self, warning: str) -> None:
        """Add a warning to the manifest."""
        self.data["warnings"].append({
            "timestamp": datetime.now().isoformat(),
            "message": str(warning)
        })
        
    def update_status(self, status: str) -> None:
        """Update the overall run status."""
        self.data["status"] = status
        
    def save(self) -> Path:
        """Save the manifest to disk and return the file path."""
        self.data["updated_at"] = datetime.now().isoformat()
        
        with open(self.manifest_file, 'w') as f:
            json.dump(self.data, f, indent=2)
            
        logger.info(f"Saved run manifest: {self.manifest_file}")
        return self.manifest_file
        
    def finalize(self, success: bool = True) -> Path:
        """
        Finalize the manifest with completion status and timings.
        
        Args:
            success: Whether the workflow completed successfully
            
        Returns:
            Path to the manifest file
        """
        self.data["end_timestamp"] = datetime.now().isoformat()
        start_time = datetime.fromisoformat(self.data["start_timestamp"])
        end_time = datetime.fromisoformat(self.data["end_timestamp"])
        duration = int((end_time - start_time).total_seconds())
        
        self.data["duration_seconds"] = duration
        self.data["status"] = "completed" if success else "failed"
        
        return self.save()


class RunManifestManager:
    """Static helper methods for run manifest management."""
    
    @staticmethod
    def create_for_workflow(gene_symbol: str, output_dir: str) -> 'RunManifest':
        """
        Create a new run manifest for a gene workflow.
        
        Args:
            gene_symbol: Gene symbol being processed
            output_dir: Output directory for this run
            
        Returns:
            New RunManifest instance
        """
        output_path = Path(output_dir)
        manifest = RunManifest(output_path, gene_symbol)
        
        # Set environment information
        import platform
        import socket
        import sys
        
        try:
            import pkg_resources
            gvf_version = pkg_resources.get_distribution("genevariantfetcher").version
        except:
            gvf_version = "development"
            
        manifest.data["environment"].update({
            "python_version": sys.version,
            "gvf_version": gvf_version,
            "hostname": socket.gethostname(),
            "working_directory": str(Path.cwd())
        })
        
        return manifest
    
    @staticmethod
    def load_from_directory(output_dir: Path) -> Optional[Dict[str, Any]]:
        """Load an existing manifest from a directory."""
        manifest_file = output_dir / "run_manifest.json"
        if manifest_file.exists():
            with open(manifest_file, 'r') as f:
                return json.load(f)
        return None
        
    @staticmethod
    def list_runs_in_directory(base_dir: Path) -> list[Dict[str, Any]]:
        """
        List all runs in a base directory (recursive).
        
        Args:
            base_dir: Base directory to search in
            
        Returns:
            List of run manifest data dicts
        """
        runs = []
        for gene_dir in base_dir.iterdir():
            if gene_dir.is_dir():
                for run_dir in gene_dir.iterdir():
                    if run_dir.is_dir():
                        manifest = RunManifestManager.load_from_directory(run_dir)
                        if manifest:
                            runs.append(manifest)
                            
        # Sort by start time, most recent first
        runs.sort(key=lambda x: x.get('start_timestamp', ''), reverse=True)
        return runs


def generate_run_manifest_cli():
    """CLI entry point for manifest inspection."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="Inspect GeneVariantFetcher run manifests"
    )
    parser.add_argument(
        "action",
        choices=["list", "show", "cleanup"],
        help="Action to perform"
    )
    parser.add_argument(
        "--directory", "-d",
        type=Path,
        default=Path("."),
        help="Directory to inspect (default: current)"
    )
    parser.add_argument(
        "--run-id", "-r",
        help="Specific run ID to show details for"
    )
    parser.add_argument(
        "--days", "-t",
        type=int,
        default=30,
        help="Age threshold in days for cleanup (default: 30)"
    )
    
    args = parser.parse_args()
    
    if args.action == "list":
        runs = RunManifestManager.list_runs_in_directory(args.directory)
        print(f"Found {len(runs)} runs:")
        for run in runs:
            print(f"  {run['gene_symbol']}: {run['run_id']}")
            print(f"    Status: {run['status']}")
            print(f"    Started: {run['start_timestamp']}")
            print(f"    Duration: {run.get('duration_seconds', 'running')}s")
            print()
            
    elif args.action == "show":
        if args.run_id:
            runs = RunManifestManager.list_runs_in_directory(args.directory)
            run = next((r for r in runs if r["run_id"] == args.run_id), None)
            if run:
                print(json.dumps(run, indent=2))
            else:
                print(f"Run {args.run_id} not found")
                sys.exit(1)
        else:
            print("--run-id required for show action")
            sys.exit(1)
            
    elif args.action == "cleanup":
        # This would be implemented when old manifest cleanup is needed
        print("Manual cleanup not implemented - use workflow cleanup commands")


if __name__ == "__main__":
    generate_run_manifest_cli()