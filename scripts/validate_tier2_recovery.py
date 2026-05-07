"""Validate the Tier 2 (InternFilter) fix on a labeled set of KCNH2 abstracts.

Compares two cohorts pulled from results/KCNH2/20260506_102238/:
  - GOLD: PMIDs in comparison_results/missing_in_sqlite.csv that were dropped
    at Tier 2 in that run. These are false negatives we need to recover.
  - CONTROL: A random sample of PMIDs that Tier 2 correctly rejected with a
    real reason (not the silent 'No reason provided' fall-through). These
    should still be FAIL after the fix; recovering too many indicates
    excessive false positives.

Reads the abstract from results/.../abstract_json/<pmid>.json and runs the
current InternFilter (which uses .env's TIER2_MODEL — Kimi-K2.6-1).
"""

from __future__ import annotations

import csv
import json
import os
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dotenv import load_dotenv

load_dotenv(REPO / ".env")

from pipeline.filters import InternFilter
from utils.models import Paper

RUN_DIR = REPO / "results/KCNH2/20260506_102238"
ABS_DIR = RUN_DIR / "abstract_json"
FILTER_JSONL = RUN_DIR / "pmid_status/filter_progress.jsonl"
GOLD_CSV = REPO / "comparison_results/missing_in_sqlite.csv"


def load_gold_dropped() -> List[str]:
    gold = set()
    with open(GOLD_CSV) as f:
        for row in csv.DictReader(f):
            if row.get("missing_in_sqlite") == "True":
                p = (row.get("pmid") or "").strip()
                if p:
                    gold.add(p)

    dropped = []
    with open(FILTER_JSONL) as f:
        for line in f:
            rec = json.loads(line)
            if (
                rec.get("pmid") in gold
                and rec.get("stage") == "tier2"
                and rec.get("final_decision") == "FAIL"
            ):
                dropped.append(rec["pmid"])
    return sorted(dropped, key=int)


def load_control_set(seed: int = 7, n: int = 30) -> List[str]:
    """Pick n PMIDs that Tier 2 rejected with a real reason (not 'No reason
    provided'). These represent papers Kimi confidently classified as off-topic
    — the fix should not flip them to PASS."""
    pool = []
    with open(FILTER_JSONL) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("stage") != "tier2" or rec.get("final_decision") != "FAIL":
                continue
            t2 = rec.get("tier2") or {}
            reason = (t2.get("reason") or "").strip()
            if not reason or reason == "No reason provided":
                continue
            if t2.get("confidence", 0) < 0.7:
                continue
            pool.append(rec["pmid"])
    rng = random.Random(seed)
    rng.shuffle(pool)
    return pool[:n]


def load_paper(pmid: str) -> Optional[Paper]:
    fp = ABS_DIR / f"{pmid}.json"
    if not fp.exists():
        return None
    data = json.loads(fp.read_text())
    title = (data.get("title") or "").strip()
    abstract = (data.get("abstract") or data.get("abstract_text") or "").strip()
    if not abstract:
        return None
    if not title:
        # The abstract dump prefixes the title within the abstract for many
        # records; use the first non-citation line as a stand-in title.
        for line in abstract.splitlines():
            line = line.strip()
            if line and not line.startswith(("1.", "PMID", "DOI", "©")):
                title = line[:200]
                break
    return Paper(pmid=pmid, title=title or pmid, abstract=abstract)


def run_filter(pmids: List[str], label: str) -> Tuple[int, int, List[Dict]]:
    flt = InternFilter()
    rows = []
    passed = 0
    failed = 0
    for i, pmid in enumerate(pmids, 1):
        paper = load_paper(pmid)
        if paper is None:
            print(f"  [{label}] {pmid}: SKIP (no abstract)")
            continue
        try:
            res = flt.filter(paper)
        except Exception as exc:  # noqa: BLE001
            print(f"  [{label}] {pmid}: ERROR {exc}")
            continue
        decision = (
            res.decision.value if hasattr(res.decision, "value") else str(res.decision)
        )
        is_pass = decision.lower() == "pass"
        if is_pass:
            passed += 1
        else:
            failed += 1
        rows.append(
            {
                "pmid": pmid,
                "decision": decision,
                "confidence": res.confidence,
                "reason": (res.reason or "")[:160],
                "fail_open": (res.metadata or {}).get("fail_open", False),
            }
        )
        marker = "PASS" if is_pass else "FAIL"
        print(
            f"  [{label}] {i:3d}/{len(pmids)} {pmid}: {marker} conf={res.confidence:.2f} "
            f"{'(fail-open)' if (res.metadata or {}).get('fail_open') else ''} "
            f"{(res.reason or '')[:90]}"
        )
    return passed, failed, rows


def main() -> int:
    gold_pmids = load_gold_dropped()
    control_pmids = load_control_set()
    print(f"GOLD (dropped, should recover): {len(gold_pmids)} PMIDs")
    print(f"CONTROL (correctly rejected, should stay FAIL): {len(control_pmids)} PMIDs")
    print()

    print("=== Running on GOLD set ===")
    g_pass, g_fail, g_rows = run_filter(gold_pmids, "gold")
    print(f"\nGOLD recovery: {g_pass}/{g_pass + g_fail} PASS  ({g_fail} still FAIL)")

    print("\n=== Running on CONTROL set ===")
    c_pass, c_fail, c_rows = run_filter(control_pmids, "ctrl")
    print(
        f"\nCONTROL leakage: {c_pass}/{c_pass + c_fail} flipped to PASS "
        f"({c_fail} still FAIL — desired)"
    )

    out = REPO / "results/tier2_validation.json"
    out.write_text(
        json.dumps(
            {
                "gold": {"pass": g_pass, "fail": g_fail, "rows": g_rows},
                "control": {"pass": c_pass, "fail": c_fail, "rows": c_rows},
            },
            indent=2,
        )
    )
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
