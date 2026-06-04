#!/usr/bin/env python3
"""Backfill papers.{first_author, journal, publication_date, doi, pmc_id} into a
scored DB from on-disk caches — NO LLM, NO network.

Sources (all local):
  * abstract_json/{pmid}.json (in run dirs)        -> authors[0], journal, year
  * corpus/<GENE>/<PMID>/{pmid}_artifacts.json     -> doi, pmcid

Only fills columns that are currently NULL/empty (never overwrites a real value).
Run after migrate / refresh, or any time, to enrich paper-level provenance for
the dashboard and the eventual variant website. Idempotent.

Usage:
  python scripts/backfill_paper_metadata.py --db results/KCNH2/<ts>/KCNH2.db
  python scripts/backfill_paper_metadata.py --db a.db --db b.db --corpus corpus
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def build_abstract_map(roots: list[Path]) -> dict[str, dict]:
    """pmid -> {first_author, journal, publication_date} from abstract_json caches."""
    out: dict[str, dict] = {}
    for root in roots:
        if not root.exists():
            continue
        for f in root.rglob("abstract_json/*.json"):
            pmid = f.stem
            if not pmid.isdigit() or pmid in out:
                continue
            try:
                m = json.loads(f.read_text()).get("metadata", {})
            except (json.JSONDecodeError, OSError):
                continue
            authors = m.get("authors") or []
            out[pmid] = {
                "first_author": authors[0] if authors else None,
                "journal": m.get("journal"),
                "publication_date": str(m["year"]) if m.get("year") else None,
            }
    return out


def build_artifact_map(corpus: Path) -> dict[str, dict]:
    """pmid -> {doi, pmc_id} from corpus artifacts.json (picks first non-null seen)."""
    out: dict[str, dict] = {}
    if not corpus.exists():
        return out
    for f in corpus.rglob("*_artifacts.json"):
        pmid = f.parent.name
        if not pmid.isdigit():
            continue
        try:
            d = json.loads(f.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        cur = out.setdefault(pmid, {"doi": None, "pmc_id": None})
        cur["doi"] = cur["doi"] or d.get("doi")
        cur["pmc_id"] = cur["pmc_id"] or d.get("pmcid")
    return out


COLS = ("first_author", "journal", "publication_date", "doi", "pmc_id")


def backfill_db(db_path: Path, abstracts: dict, artifacts: dict) -> dict:
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    have = {r[1] for r in con.execute("PRAGMA table_info(papers)")}
    for col in COLS:
        if col not in have:
            con.execute(f"ALTER TABLE papers ADD COLUMN {col} TEXT")
    filled = {c: 0 for c in COLS}
    for row in con.execute("SELECT * FROM papers").fetchall():
        pmid = str(row["pmid"])
        src = {
            **abstracts.get(pmid, {}),
            **{k: v for k, v in artifacts.get(pmid, {}).items()},
        }
        sets, vals = [], []
        for col in COLS:
            cur = row[col] if col in row.keys() else None
            new = src.get(col)
            if new and not (cur and str(cur).strip()):
                sets.append(f"{col} = ?")
                vals.append(new)
                filled[col] += 1
        if sets:
            vals.append(pmid)
            con.execute(f"UPDATE papers SET {', '.join(sets)} WHERE pmid = ?", vals)
    con.commit()
    con.close()
    return filled


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--db", action="append", required=True, help="DB path (repeatable)."
    )
    ap.add_argument(
        "--corpus",
        default=str(REPO / "corpus"),
        help="Corpus dir (artifacts.json source).",
    )
    ap.add_argument(
        "--roots",
        nargs="*",
        default=["results", "validation_runs"],
        help="Roots to scan for abstract_json/.",
    )
    args = ap.parse_args()

    abstracts = build_abstract_map([REPO / r for r in args.roots])
    artifacts = build_artifact_map(Path(args.corpus).expanduser())
    print(
        f"sources: {len(abstracts)} abstract_json metadata, {len(artifacts)} artifacts.json"
    )
    for db in args.db:
        p = Path(db).expanduser()
        if not p.exists():
            print(f"  SKIP (missing): {p}")
            continue
        filled = backfill_db(p, abstracts, artifacts)
        print(f"  {p}: " + ", ".join(f"{c}+{n}" for c, n in filled.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
