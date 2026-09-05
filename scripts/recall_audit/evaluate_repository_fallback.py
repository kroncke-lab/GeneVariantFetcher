#!/usr/bin/env python3
"""Live source-only evaluation on an explicit PMID panel; no gold values or LLMs.

All downloads are isolated. Use build_source_corpus.py --dry-run (its default),
then --apply to sync acceptable source upgrades. Never rewrite scored runs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.repository_pdf_fallback import (
    RepositoryPDFRecovery,
    write_repository_source,
)  # noqa: E402


def evaluate(panel_path, out, email, workers=2):
    config = json.loads(panel_path.read_text())
    papers = config["papers"]
    if len({str(p["pmid"]) for p in papers}) != len(papers):
        raise ValueError("panel must contain unique PMIDs")
    if out.exists() and any(out.iterdir()):
        raise ValueError("use a fresh output directory; prior evidence is immutable")
    out.mkdir(parents=True, exist_ok=True)
    (out / "panel.json").write_text(json.dumps(config, indent=2) + "\n")

    def run(paper):
        started = time.monotonic()
        gene, pmid = paper["gene"], str(paper["pmid"])
        root = out / "sources" / gene
        result = RepositoryPDFRecovery(email=email).recover(pmid=pmid, output_dir=root)
        if result.success:
            write_repository_source(result, root, gene=gene)
        row = {
            "gene": gene,
            "pmid": pmid,
            "role": paper["role"],
            "success": result.success,
            "doi": result.doi,
            "title": result.title,
            "source_url": result.source_url,
            "pdf_sha256": result.pdf_sha256,
            "pages": result.page_count,
            "characters": len(result.markdown),
            "provider": result.candidate.get("provider", ""),
            "identity_reason": result.identity_reason,
            "error": result.error,
            "seconds": round(time.monotonic() - started, 2),
            "source_status": "body_only" if result.success else "not_recovered",
        }
        print(json.dumps(row), flush=True)
        return row

    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(run, papers))
    with (out / "results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "panel_sha256": hashlib.sha256(panel_path.read_bytes()).hexdigest(),
        "purpose": "source-only live network evaluation; no extraction score or phenotype completeness claim",
        "attempted": len(rows),
        "recovered": sum(r["success"] for r in rows),
        "additional_attempted": sum(
            r["role"] != "requested reproduction" for r in rows
        ),
        "additional_recovered": sum(
            r["success"] and r["role"] != "requested reproduction" for r in rows
        ),
        "runtime_sha256": {
            name: hashlib.sha256((REPO / name).read_bytes()).hexdigest()
            for name in (
                "harvesting/repository_pdf_fallback.py",
                "harvesting/unpaywall_api.py",
                "harvesting/orchestrator.py",
                "harvesting/free_text_flow.py",
                "scripts/fetch_paywalled.py",
            )
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--email", required=True)
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=2)
    args = parser.parse_args()
    print(
        json.dumps(
            evaluate(
                args.panel.resolve(), args.out_dir.resolve(), args.email, args.workers
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
