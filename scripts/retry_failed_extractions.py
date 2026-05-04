"""Retry failed extractions sequentially with rate-limit-aware backoff.

Default Tier 3 model is selected via TIER3_MODELS env. When the upstream
deployment is broken (e.g. azure_ai/grok-4-20-reasoning returns
InternalServerError on every call), pass --models to override.
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.extraction import ExpertExtractor
from utils.models import Paper


_RATE_LIMIT_MARKERS = ("RateLimitError", "RateLimitReached")
_SERVICE_DOWN_MARKERS = (
    "InternalServerError",
    "Model service is unavailable",
    "ServiceUnavailable",
)


def _is_rate_limit(text: str) -> bool:
    return any(m in text for m in _RATE_LIMIT_MARKERS)


def _is_service_down(text: str) -> bool:
    return any(m in text for m in _SERVICE_DOWN_MARKERS)


def main(
    run_dir: Path,
    gene: str = "KCNH2",
    base_delay: float = 5.0,
    rate_limit_cooldown: float = 90.0,
    service_down_cooldown: float = 180.0,
    max_consecutive_pauses: int = 5,
    models: list[str] | None = None,
) -> None:
    failures_csv = run_dir / "pmid_status" / "extraction_failures.csv"
    extractions_dir = run_dir / "extractions"
    harvest_dir = run_dir / "pmc_fulltext"
    abstract_dir = run_dir / "abstract_json"

    failed_pmids: list[str] = []
    with open(failures_csv) as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if row:
                failed_pmids.append(row[0])

    already_done = {
        p.name.replace(f"{gene}_PMID_", "").replace(".json", "")
        for p in extractions_dir.glob(f"{gene}_PMID_*.json")
    }
    pending = [p for p in failed_pmids if p not in already_done]
    print(
        f"failed in csv: {len(failed_pmids)}, already extracted: "
        f"{len(already_done & set(failed_pmids))}, pending retry: {len(pending)}",
        flush=True,
    )
    if models:
        print(f"using Tier 3 models: {models}", flush=True)

    extractor = ExpertExtractor(
        models=models, tier_threshold=1, fulltext_dir=str(harvest_dir)
    )

    success = 0
    fail = 0
    new_failures: list[tuple[str, str]] = []
    consecutive_pauses = 0
    i = 0
    queue = list(pending)
    while queue:
        i += 1
        pmid = queue.pop(0)

        ds = harvest_dir / f"{pmid}_DATA_ZONES.md"
        cl = harvest_dir / f"{pmid}_CLEANED.md"
        fc = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
        ab = abstract_dir / f"{pmid}.json"

        md_path = (
            ds if ds.exists() else cl if cl.exists() else fc if fc.exists() else None
        )

        try:
            if md_path is not None:
                content = md_path.read_text(encoding="utf-8")
                paper = Paper(pmid=pmid, full_text=content, gene_symbol=gene)
                result = extractor.extract(paper)
                source = "fulltext"
            elif ab.exists():
                with open(ab) as f:
                    record = json.load(f)
                abstract_text = record.get("abstract")
                if not abstract_text:
                    new_failures.append((pmid, "no abstract"))
                    fail += 1
                    continue
                metadata = record.get("metadata", {})
                paper = Paper(
                    pmid=pmid,
                    title=metadata.get("title", f"Paper {pmid}"),
                    abstract=abstract_text,
                    gene_symbol=gene,
                )
                result = extractor.extract(paper)
                if result.success and result.extracted_data is not None:
                    result.extracted_data.setdefault("extraction_metadata", {})[
                        "abstract_only"
                    ] = True
                source = "abstract"
            else:
                new_failures.append((pmid, "no source files"))
                fail += 1
                continue

            if result.success:
                out = extractions_dir / f"{gene}_PMID_{pmid}.json"
                with open(out, "w") as f:
                    json.dump(result.extracted_data, f, indent=2)
                success += 1
                consecutive_pauses = 0
                print(f"[{i}/{i + len(queue)}] {pmid} OK ({source})", flush=True)
            else:
                err = result.error or "unknown"
                if _is_rate_limit(err):
                    if consecutive_pauses >= max_consecutive_pauses:
                        new_failures.append((pmid, err))
                        fail += 1
                        print(
                            f"[{i}/{i + len(queue)}] {pmid} RL-ABANDONED after "
                            f"{max_consecutive_pauses} consecutive cooldowns",
                            flush=True,
                        )
                    else:
                        consecutive_pauses += 1
                        print(
                            f"[{i}/{i + len(queue)}] {pmid} RATE-LIMITED, "
                            f"sleeping {rate_limit_cooldown:.0f}s "
                            f"(pause {consecutive_pauses}/{max_consecutive_pauses}) "
                            f"and re-queueing",
                            flush=True,
                        )
                        queue.insert(0, pmid)
                        i -= 1
                        time.sleep(rate_limit_cooldown)
                        continue
                elif _is_service_down(err):
                    if consecutive_pauses >= max_consecutive_pauses:
                        new_failures.append((pmid, err))
                        fail += 1
                        print(
                            f"[{i}/{i + len(queue)}] {pmid} SVC-ABANDONED",
                            flush=True,
                        )
                    else:
                        consecutive_pauses += 1
                        print(
                            f"[{i}/{i + len(queue)}] {pmid} SERVICE-DOWN, "
                            f"sleeping {service_down_cooldown:.0f}s "
                            f"(pause {consecutive_pauses}/{max_consecutive_pauses}) "
                            f"and re-queueing",
                            flush=True,
                        )
                        queue.insert(0, pmid)
                        i -= 1
                        time.sleep(service_down_cooldown)
                        continue
                else:
                    new_failures.append((pmid, err))
                    fail += 1
                    consecutive_pauses = 0
                    print(
                        f"[{i}/{i + len(queue)}] {pmid} FAIL: {err[:200]}",
                        flush=True,
                    )
        except Exception as e:
            err = str(e)
            new_failures.append((pmid, err))
            fail += 1
            consecutive_pauses = 0
            print(f"[{i}/{i + len(queue)}] {pmid} EXC: {err[:200]}", flush=True)

        time.sleep(base_delay)

    print(f"\nDone. Success: {success}, Fail: {fail}", flush=True)
    with open(
        run_dir / "pmid_status" / "extraction_retry_failures.csv", "w", newline=""
    ) as f:
        w = csv.writer(f)
        w.writerow(["PMID", "Error"])
        w.writerows(new_failures)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--gene", default="KCNH2")
    parser.add_argument("--base-delay", type=float, default=5.0)
    parser.add_argument("--rate-limit-cooldown", type=float, default=90.0)
    parser.add_argument("--service-down-cooldown", type=float, default=180.0)
    parser.add_argument("--max-consecutive-pauses", type=int, default=5)
    parser.add_argument(
        "--models",
        nargs="+",
        help="Override Tier 3 model list (e.g. azure_ai/Kimi-K2.6-1)",
    )
    args = parser.parse_args()
    main(
        args.run_dir,
        gene=args.gene,
        base_delay=args.base_delay,
        rate_limit_cooldown=args.rate_limit_cooldown,
        service_down_cooldown=args.service_down_cooldown,
        max_consecutive_pauses=args.max_consecutive_pauses,
        models=args.models,
    )
