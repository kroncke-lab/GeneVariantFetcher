#!/usr/bin/env python3
"""Backfill papers.{first_author, journal, publication_date, doi, pmc_id} into a
scored DB. Local caches first (NO network), with an optional PubMed fallback for
the columns the local caches don't carry (mainly doi / pmc_id).

Local sources (no network, no LLM):
  * abstract_json/{pmid}.json (in run dirs)        -> authors[0], journal, year
  * corpus/<GENE>/<PMID>/{pmid}_artifacts.json     -> doi, pmcid

Optional network fallback (``--fetch-missing``):
  * PubMed ESummary for bibliographic columns still missing after the local pass.
  * PubMed EFetch for first-author affiliation / author-country fields.

Only fills columns that are currently NULL/empty (never overwrites a real value),
so it is idempotent and safe to re-run after migrate / refresh. Local values always
win over fetched ones.

Usage:
  python scripts/backfill_paper_metadata.py --db results/KCNH2/<ts>/KCNH2.db
  python scripts/backfill_paper_metadata.py --db a.db --fetch-missing --email you@x.org
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

REPO = Path(__file__).resolve().parents[1]
# Allow ``from utils...`` when run as a standalone script from scripts/.
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# ESummary-derivable columns (also sourced from local abstract/artifact caches).
COLS = ("first_author", "journal", "publication_date", "doi", "pmc_id")
# Author-affiliation columns — a weaker paper-level location signal (EFetch only,
# never in local caches). Downstream consumers should keep this distinct from a
# stated cohort origin.
AUTHOR_COLS = ("author_affiliation", "author_country")
ALL_COLS = COLS + AUTHOR_COLS


def _clean(value: Any) -> Optional[str]:
    if value in (None, ""):
        return None
    text = str(value).strip()
    return text or None


def _layer(*maps: Dict[str, Any]) -> Dict[str, Any]:
    """Merge column maps; later maps win, but only with a non-empty value.

    Lets us pass ``_layer(fetched, abstracts, artifacts)`` so the local caches
    override PubMed where present, and PubMed only fills the gaps.
    """
    out: Dict[str, Any] = {}
    for source in maps:
        for key, value in (source or {}).items():
            cleaned = _clean(value)
            if cleaned is not None:
                out[key] = cleaned
    return out


def esummary_to_columns(md: Dict[str, Any]) -> Dict[str, Optional[str]]:
    """Map one PubMed ESummary document to our ``papers`` columns."""
    article_ids = md.get("ArticleIds") or {}
    if not isinstance(article_ids, dict):
        article_ids = {}
    authors = md.get("AuthorList") or []
    pmc = _clean(article_ids.get("pmc"))
    if pmc and not pmc.upper().startswith("PMC"):
        pmc = f"PMC{pmc}"
    return {
        "first_author": _clean(authors[0]) if authors else None,
        "journal": _clean(md.get("FullJournalName")) or _clean(md.get("Source")),
        "publication_date": _clean(md.get("PubDate"))
        or _clean(md.get("EPubDate"))
        or _clean(md.get("SortPubDate")),
        "doi": _clean(md.get("DOI")) or _clean(article_ids.get("doi")),
        "pmc_id": pmc,
    }


def build_abstract_map(roots: List[Path]) -> Dict[str, dict]:
    """pmid -> {first_author, journal, publication_date} from abstract_json caches."""
    out: Dict[str, dict] = {}
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


def build_artifact_map(corpus: Path) -> Dict[str, dict]:
    """pmid -> {doi, pmc_id} from corpus artifacts.json (picks first non-null seen)."""
    out: Dict[str, dict] = {}
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


def db_pmids(db_path: Path) -> List[str]:
    """Numeric PMIDs present in a DB's papers table."""
    con = sqlite3.connect(str(db_path))
    try:
        rows = con.execute("SELECT pmid FROM papers").fetchall()
    except sqlite3.OperationalError:
        return []
    finally:
        con.close()
    return [str(r[0]) for r in rows if r[0] is not None and str(r[0]).strip().isdigit()]


def pmids_missing_after_local(
    pmids: Iterable[str],
    abstracts: Dict[str, dict],
    artifacts: Dict[str, dict],
) -> List[str]:
    """PMIDs still missing at least one column once local caches are applied."""
    # author_* columns are never in local caches, so any PMID is "missing" them —
    # that's intentional: it makes fetch_missing pull affiliations for every paper.
    missing: List[str] = []
    for pmid in pmids:
        have = _layer(abstracts.get(pmid, {}), artifacts.get(pmid, {}))
        if any(col not in have for col in ALL_COLS):
            missing.append(pmid)
    return missing


def build_fetched_map(
    pmids: List[str], email: Optional[str], max_fetch: Optional[int] = None
) -> Dict[str, dict]:
    """pmid -> columns from PubMed ESummary + author affiliation (network)."""
    if not pmids:
        return {}
    if max_fetch is not None:
        pmids = pmids[:max_fetch]
    from utils.geo_ancestry import country_from_affiliation
    from utils.pubmed_utils import batch_fetch_affiliations, batch_fetch_metadata

    pmids = list(pmids)
    out = {
        pmid: esummary_to_columns(doc)
        for pmid, doc in batch_fetch_metadata(pmids, email=email).items()
    }
    for pmid, affiliation in batch_fetch_affiliations(pmids, email=email).items():
        entry = out.setdefault(pmid, {})
        entry["author_affiliation"] = affiliation
        entry["author_country"] = country_from_affiliation(affiliation)
    return out


def backfill_db(
    db_path: Path,
    abstracts: Dict[str, dict],
    artifacts: Dict[str, dict],
    fetched: Optional[Dict[str, dict]] = None,
) -> Dict[str, int]:
    fetched = fetched or {}
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    have = {r[1] for r in con.execute("PRAGMA table_info(papers)")}
    for col in ALL_COLS:
        if col not in have:
            con.execute(f"ALTER TABLE papers ADD COLUMN {col} TEXT")
    filled = {c: 0 for c in ALL_COLS}
    for row in con.execute("SELECT * FROM papers").fetchall():
        pmid = str(row["pmid"])
        # Local caches win over PubMed; PubMed only fills what's still missing.
        src = _layer(
            fetched.get(pmid, {}),
            abstracts.get(pmid, {}),
            artifacts.get(pmid, {}),
        )
        sets, vals = [], []
        for col in ALL_COLS:
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


def run_backfill(
    dbs: List[str],
    *,
    corpus: Optional[Path] = None,
    roots: Optional[List[Path]] = None,
    fetch_missing: bool = False,
    email: Optional[str] = None,
    max_fetch: Optional[int] = None,
) -> Dict[str, Dict[str, int]]:
    """Backfill one or more DBs. Reusable entry point (also called from gvf-run)."""
    corpus = corpus or (REPO / "corpus")
    roots = roots or [REPO / "results", REPO / "validation_runs"]
    abstracts = build_abstract_map(roots)
    artifacts = build_artifact_map(corpus)

    fetched: Dict[str, dict] = {}
    if fetch_missing:
        candidates: set[str] = set()
        for db in dbs:
            p = Path(db).expanduser()
            if p.exists():
                candidates |= set(
                    pmids_missing_after_local(db_pmids(p), abstracts, artifacts)
                )
        fetched = build_fetched_map(sorted(candidates), email, max_fetch)

    results: Dict[str, Dict[str, int]] = {}
    for db in dbs:
        p = Path(db).expanduser()
        if not p.exists():
            continue
        results[db] = backfill_db(p, abstracts, artifacts, fetched)
    return results


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
    ap.add_argument(
        "--fetch-missing",
        action="store_true",
        help="Fall back to PubMed ESummary for columns local caches lack (doi/pmc_id).",
    )
    ap.add_argument(
        "--email", help="NCBI contact email (else NCBI_EMAIL / Entrez default)."
    )
    ap.add_argument(
        "--max-fetch", type=int, default=None, help="Cap PMIDs fetched from PubMed."
    )
    args = ap.parse_args()

    results = run_backfill(
        args.db,
        corpus=Path(args.corpus).expanduser(),
        roots=[REPO / r for r in args.roots],
        fetch_missing=args.fetch_missing,
        email=args.email,
        max_fetch=args.max_fetch,
    )
    for db, filled in results.items():
        print(f"  {db}: " + ", ".join(f"{c}+{n}" for c, n in filled.items()))
    if not results:
        print("  (no DBs updated)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
