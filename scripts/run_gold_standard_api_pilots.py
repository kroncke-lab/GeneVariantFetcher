#!/usr/bin/env python3
"""Run API smoke checks over the deterministic gold-standard pilot PMIDs.

Default mode is a dry plan: it verifies the pilot package and reports which
checks would run. Use ``--live`` to exercise project API client classes. Live
mode requires the GVF Python dependencies to be installed in the active
environment.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


BASE_DIR = Path(__file__).resolve().parents[1]
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))
DEFAULT_PILOT_DIR = BASE_DIR / "gene_variant_fetcher_gold_standard" / "pilots"
DEFAULT_OUT = BASE_DIR / "comparison_results" / "gold_standard_api_pilots"
CHECK_COLUMNS = (
    "gene",
    "pmid",
    "pilot_case",
    "pubmed_metadata",
    "europepmc_metadata",
    "pmcid",
    "doi",
    "unpaywall_status",
    "publisher_route",
    "publisher_client_configured",
    "error",
)


def load_local_env() -> None:
    """Load .env when python-dotenv is installed; dry runs do not require it."""
    try:
        from dotenv import load_dotenv
    except Exception:
        return
    load_dotenv(BASE_DIR / ".env")


@dataclass
class ApiPilotResult:
    gene: str
    pmid: str
    pilot_case: str
    pubmed_metadata: str = "not_run"
    europepmc_metadata: str = "not_run"
    pmcid: str = ""
    doi: str = ""
    unpaywall_status: str = "not_run"
    publisher_route: str = ""
    publisher_client_configured: str = ""
    error: str = ""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[ApiPilotResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHECK_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def filter_rows(
    rows: list[dict[str, str]], genes: set[str] | None, limit_per_gene: int | None
) -> list[dict[str, str]]:
    filtered = [
        row for row in rows if genes is None or row.get("gene", "").upper() in genes
    ]
    if limit_per_gene is None:
        return filtered

    kept: list[dict[str, str]] = []
    counts: dict[str, int] = {}
    for row in filtered:
        gene = row.get("gene", "").upper()
        if counts.get(gene, 0) >= limit_per_gene:
            continue
        kept.append(row)
        counts[gene] = counts.get(gene, 0) + 1
    return kept


def parse_gene_filter(value: str | None) -> set[str] | None:
    if not value:
        return None
    genes = {item.strip().upper() for item in value.split(",") if item.strip()}
    return genes or None


def import_project_clients():
    """Import the project API classes lazily so dry-run has no dependencies."""
    from harvesting.elsevier_api import ElsevierAPIClient
    from harvesting.pmc_api import PMCAPIClient
    from harvesting.springer_api import SpringerAPIClient
    from harvesting.unpaywall_api import UnpaywallClient
    from harvesting.wiley_api import WileyAPIClient
    from utils import pubmed_utils

    return {
        "PMCAPIClient": PMCAPIClient,
        "UnpaywallClient": UnpaywallClient,
        "ElsevierAPIClient": ElsevierAPIClient,
        "SpringerAPIClient": SpringerAPIClient,
        "WileyAPIClient": WileyAPIClient,
        "pubmed_utils": pubmed_utils,
    }


def publisher_route(doi: str, clients: dict) -> tuple[str, str]:
    if not doi:
        return "", ""

    configured_clients = {
        "elsevier": clients["ElsevierAPIClient"](
            api_key=os.getenv("ELSEVIER_API_KEY"),
            insttoken=os.getenv("ELSEVIER_INSTTOKEN"),
        ),
        "springer": clients["SpringerAPIClient"](api_key=os.getenv("SPRINGER_API_KEY")),
        "wiley": clients["WileyAPIClient"](api_key=os.getenv("WILEY_API_KEY")),
    }
    checks = {
        "elsevier": configured_clients["elsevier"].is_elsevier_doi,
        "springer": configured_clients["springer"].is_springer_doi,
        "wiley": configured_clients["wiley"].is_wiley_doi,
    }
    matches = [name for name, check in checks.items() if check(doi)]
    if len(matches) > 1:
        configured = ",".join(
            f"{name}={configured_clients[name].is_available}" for name in matches
        )
        return f"ambiguous:{'|'.join(matches)}", configured
    if matches:
        name = matches[0]
        return name, str(configured_clients[name].is_available)
    return "generic_or_unknown", ""


def run_live_check(
    row: dict[str, str], clients: dict, email: str | None
) -> ApiPilotResult:
    pmid = row["pmid"]
    result = ApiPilotResult(
        gene=row["gene"],
        pmid=pmid,
        pilot_case=row.get("pilot_case", ""),
    )

    try:
        pubmed_metadata = clients["pubmed_utils"].fetch_paper_metadata(pmid)
        result.pubmed_metadata = "ok" if pubmed_metadata else "missing"
    except Exception as exc:  # noqa: BLE001
        result.pubmed_metadata = "error"
        result.error = f"pubmed:{exc}"

    try:
        # Europe PMC coverage via the production discovery path (gene search):
        # does Europe PMC surface this PMID for its gene? (pmcid/doi are filled
        # from the PMC API block below.)
        epmc_pmids = clients["pubmed_utils"].query_europepmc(
            row["gene"], max_results=1000
        )
        result.europepmc_metadata = "ok" if pmid in epmc_pmids else "missing"
    except Exception as exc:  # noqa: BLE001
        result.europepmc_metadata = "error"
        result.error = append_error(result.error, f"europepmc:{exc}")

    try:
        pmc = clients["PMCAPIClient"]()
        pmcid = pmc.pmid_to_pmcid(pmid)
        if pmcid and not result.pmcid:
            result.pmcid = pmcid
        doi = pmc.get_doi_from_pmid(pmid)
        if doi and not result.doi:
            result.doi = doi
    except Exception as exc:  # noqa: BLE001
        result.error = append_error(result.error, f"pmc:{exc}")

    result.publisher_route, result.publisher_client_configured = publisher_route(
        result.doi, clients
    )

    if email and result.doi:
        try:
            unpaywall = clients["UnpaywallClient"](email=email)
            oa_result, error = unpaywall.find_open_access(result.doi)
            if error:
                result.unpaywall_status = f"error:{error}"
            elif oa_result:
                result.unpaywall_status = (
                    f"oa={oa_result.get('is_oa')}:{oa_result.get('oa_status')}"
                )
            else:
                result.unpaywall_status = "missing"
        except Exception as exc:  # noqa: BLE001
            result.unpaywall_status = "error"
            result.error = append_error(result.error, f"unpaywall:{exc}")
    elif result.doi:
        result.unpaywall_status = "skipped:no_email"
    else:
        result.unpaywall_status = "skipped:no_doi"

    return result


def append_error(existing: str, new: str) -> str:
    return f"{existing}; {new}" if existing else new


def run_dry_check(row: dict[str, str]) -> ApiPilotResult:
    return ApiPilotResult(
        gene=row["gene"],
        pmid=row["pmid"],
        pilot_case=row.get("pilot_case", ""),
        pubmed_metadata="planned",
        europepmc_metadata="planned",
        unpaywall_status="planned_if_doi_and_email",
        publisher_route="planned_after_doi",
        publisher_client_configured="planned",
    )


def main() -> int:
    load_local_env()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-dir",
        type=Path,
        default=DEFAULT_PILOT_DIR,
        help="Pilot package built by scripts/build_gold_standard_pilots.py.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT,
        help="Directory for API pilot results.",
    )
    parser.add_argument("--genes", help="Comma-separated gene filter.")
    parser.add_argument(
        "--limit-per-gene",
        type=int,
        default=5,
        help="Maximum pilot PMIDs per gene to check. Use 0 for all.",
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="Exercise project API client classes instead of writing a dry plan.",
    )
    parser.add_argument(
        "--email",
        default=os.getenv("UNPAYWALL_EMAIL") or os.getenv("NCBI_EMAIL"),
        help="Email for Unpaywall. Defaults to UNPAYWALL_EMAIL or NCBI_EMAIL.",
    )
    args = parser.parse_args()

    pilot_pmids = args.pilot_dir.expanduser() / "pilot_pmids.csv"
    if not pilot_pmids.exists():
        parser.error(
            f"Missing {pilot_pmids}. Run scripts/build_gold_standard_pilots.py first."
        )

    limit = None if args.limit_per_gene == 0 else args.limit_per_gene
    rows = filter_rows(read_csv(pilot_pmids), parse_gene_filter(args.genes), limit)
    if not rows:
        parser.error("No pilot PMIDs selected.")

    clients = None
    import_error = ""
    if args.live:
        try:
            clients = import_project_clients()
        except Exception as exc:  # noqa: BLE001
            import_error = str(exc)

    if args.live and clients is None:
        results = [
            ApiPilotResult(
                gene=row["gene"],
                pmid=row["pmid"],
                pilot_case=row.get("pilot_case", ""),
                error=f"project-client imports failed: {import_error}",
            )
            for row in rows
        ]
        mode = "live_import_failed"
        exit_code = 2
    elif args.live:
        results = [run_live_check(row, clients, args.email) for row in rows]
        mode = "live_project_clients"
        exit_code = 0
    else:
        results = [run_dry_check(row) for row in rows]
        mode = "dry_plan"
        exit_code = 0

    out_dir = args.out_dir.expanduser()
    write_csv(out_dir / "api_pilot_results.csv", results)
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": mode,
        "pilot_dir": str(args.pilot_dir.expanduser()),
        "selected_pmid_count": len(rows),
        "genes": sorted({row["gene"] for row in rows}),
        "live": args.live,
        "email_configured": bool(args.email),
        "import_error": import_error,
        "results_csv": str(out_dir / "api_pilot_results.csv"),
    }
    write_json(out_dir / "summary.json", summary)

    print(f"Mode: {mode}")
    print(f"PMIDs checked: {len(rows)}")
    print(f"Results: {out_dir / 'api_pilot_results.csv'}")
    if import_error:
        print(f"Import error: {import_error}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
