"""Retry failed extractions sequentially to avoid Azure rate-limit thrash."""

import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.extraction import ExpertExtractor
from utils.models import Paper


def main(run_dir: Path, gene: str = "KCNH2", min_delay: float = 2.0) -> None:
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
        f"failed in csv: {len(failed_pmids)}, already extracted: {len(already_done & set(failed_pmids))}, pending retry: {len(pending)}"
    )

    extractor = ExpertExtractor(tier_threshold=1, fulltext_dir=str(harvest_dir))

    success = 0
    fail = 0
    new_failures: list[tuple[str, str]] = []

    for i, pmid in enumerate(pending, 1):
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
                print(f"[{i}/{len(pending)}] {pmid} OK ({source})")
            else:
                new_failures.append((pmid, result.error or "unknown"))
                fail += 1
                print(f"[{i}/{len(pending)}] {pmid} FAIL: {result.error}")
        except Exception as e:
            new_failures.append((pmid, str(e)))
            fail += 1
            print(f"[{i}/{len(pending)}] {pmid} EXC: {e}")

        time.sleep(min_delay)

    print(f"\nDone. Success: {success}, Fail: {fail}")
    with open(
        run_dir / "pmid_status" / "extraction_retry_failures.csv", "w", newline=""
    ) as f:
        w = csv.writer(f)
        w.writerow(["PMID", "Error"])
        w.writerows(new_failures)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd())
