import json
import os
from typing import Dict, List


def get_pmid_status(output_dir: str, pmid: str) -> Dict:
    """
    Returns the status of a PMID from its JSON file.
    """
    status_dir = os.path.join(output_dir, "pmid_status")
    status_file = os.path.join(status_dir, f"{pmid}.json")
    try:
        with open(status_file, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def get_failed_pmids(output_dir: str) -> List[str]:
    """
    Returns a list of failed PMIDs (includes 'failed' and 'paywalled' statuses).
    """
    status_dir = os.path.join(output_dir, "pmid_status")
    failed_pmids = []
    failed_statuses = {"failed", "paywalled"}
    for filename in os.listdir(status_dir):
        if filename.endswith(".json"):
            pmid = filename[:-5]  # Remove ".json"
            status = get_pmid_status(output_dir, pmid)
            if status and status.get("status") in failed_statuses:
                failed_pmids.append(pmid)
    return failed_pmids


def get_stats_summary(output_dir: str) -> Dict[str, int]:
    """
    Returns a summary of the status counts.
    """
    status_dir = os.path.join(output_dir, "pmid_status")
    status_counts = {
        "downloaded": 0,
        "extracted": 0,
        "failed": 0,
        "paywalled": 0,
        "no_variants": 0,
    }
    for filename in os.listdir(status_dir):
        if filename.endswith(".json"):
            pmid = filename[:-5]  # Remove ".json"
            status = get_pmid_status(output_dir, pmid)
            if status:
                status_counts[status["status"]] += 1
    return status_counts
