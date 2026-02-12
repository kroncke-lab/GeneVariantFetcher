"""Paywall audit utilities.

Goal: distinguish truly paywalled papers from papers that are free/open but blocked
by captcha/anti-bot or missed by our download logic.

This is intentionally lightweight and safe to run after a harvest. It does NOT
attempt downloads; it only queries metadata sources and classifies outcomes.
"""

from __future__ import annotations

import csv
import json
import os
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests


@dataclass
class AuditResult:
    pmid: str
    doi: Optional[str]
    original_reason: Optional[str]
    classification: str
    unpaywall: Dict[str, Any]
    europepmc: Dict[str, Any]
    doi_landing: Dict[str, Any]


CAPTCHA_KEYWORDS = [
    "captcha",
    "recaptcha",
    "hcaptcha",
    "cloudflare",
    "please verify",
    "robot",
    "challenge",
    "access denied",
]

PAYWALL_KEYWORDS = [
    "purchase",
    "subscribe",
    "buy this article",
    "rent this article",
    "institutional access",
    "sign in to access",
    "get access",
]


def _session(ncbi_email: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": f"GVF-PaywallAudit/1.0 (mailto:{ncbi_email})"})
    return s


def load_paywalled_pmids(paywalled_csv: Path) -> Dict[str, str]:
    pmids: Dict[str, str] = {}
    if not paywalled_csv.exists():
        return pmids
    with open(paywalled_csv, newline="") as f:
        r = csv.reader(f)
        _ = next(r, None)  # header
        for row in r:
            if not row:
                continue
            pmid = str(row[0]).strip()
            reason = str(row[1]).strip() if len(row) > 1 else ""
            if pmid:
                pmids[pmid] = reason
    return pmids


def pubmed_get_doi(sess: requests.Session, pmid: str) -> Optional[str]:
    try:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=pubmed&id={pmid}&retmode=json"
        )
        r = sess.get(url, timeout=15)
        r.raise_for_status()
        data = r.json()
        item = data.get("result", {}).get(str(pmid), {})
        for artid in item.get("articleids", []) or []:
            if artid.get("idtype") == "doi":
                return artid.get("value")
    except Exception:
        return None
    return None


def unpaywall_check(sess: requests.Session, doi: Optional[str], email: str) -> Dict[str, Any]:
    if not doi:
        return {}
    try:
        url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
        r = sess.get(url, timeout=15)
        if r.status_code != 200:
            return {"error": f"HTTP {r.status_code}"}
        data = r.json()
        best = data.get("best_oa_location") or {}
        return {
            "is_oa": data.get("is_oa"),
            "oa_status": data.get("oa_status"),
            "journal_is_oa": data.get("journal_is_oa"),
            "best_oa_url": best.get("url"),
            "best_oa_host": best.get("host_type"),
        }
    except Exception as e:
        return {"error": str(e)}


def europepmc_check(sess: requests.Session, pmid: str) -> Dict[str, Any]:
    try:
        url = (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
            f"?query=EXT_ID:{pmid}&format=json&resultType=core"
        )
        r = sess.get(url, timeout=15)
        if r.status_code != 200:
            return {"error": f"HTTP {r.status_code}"}
        data = r.json()
        results = data.get("resultList", {}).get("result", []) or []
        if not results:
            return {}
        paper = results[0]
        urls = [
            {"url": u.get("url"), "availabilityCode": u.get("availabilityCode")}
            for u in (paper.get("fullTextUrlList", {}) or {}).get("fullTextUrl", [])
        ]
        return {
            "isOpenAccess": paper.get("isOpenAccess"),
            "inEPMC": paper.get("inEPMC"),
            "inPMC": paper.get("inPMC"),
            "hasPDF": paper.get("hasPDF"),
            "fullTextUrls": urls,
        }
    except Exception as e:
        return {"error": str(e)}


def doi_landing_check(sess: requests.Session, doi: Optional[str]) -> Dict[str, Any]:
    if not doi:
        return {}
    try:
        r = sess.get(f"https://doi.org/{doi}", timeout=20, allow_redirects=True)
        text = (r.text or "")[:6000].lower()
        has_captcha = any(k in text for k in CAPTCHA_KEYWORDS)
        has_paywall = any(k in text for k in PAYWALL_KEYWORDS)
        return {
            "status_code": r.status_code,
            "final_url": r.url,
            "has_captcha": has_captcha,
            "has_paywall_language": has_paywall,
            "content_length": len(r.text or ""),
        }
    except Exception as e:
        return {"error": str(e)}


def classify(unpaywall: Dict[str, Any], europepmc: Dict[str, Any], landing: Dict[str, Any]) -> str:
    if unpaywall.get("is_oa"):
        return "FREE_BUT_BLOCKED" if landing.get("has_captcha") else "FREE_BUT_MISSED"

    if europepmc.get("isOpenAccess") == "Y" or europepmc.get("inPMC") == "Y":
        return "FREE_BUT_BLOCKED" if landing.get("has_captcha") else "FREE_BUT_MISSED"

    if landing.get("has_captcha") and not landing.get("has_paywall_language"):
        return "FREE_BUT_BLOCKED"

    if landing.get("has_paywall_language"):
        return "TRULY_PAYWALLED"

    return "LIKELY_PAYWALLED"


def run_paywall_audit(
    paywalled_csv: Path,
    out_json: Path,
    out_md: Path,
    ncbi_email: Optional[str] = None,
    sleep_s: float = 0.35,
    limit: Optional[int] = None,
) -> Dict[str, Any]:
    email = (
        ncbi_email
        or os.environ.get("NCBI_EMAIL")
        or os.environ.get("UNPAYWALL_EMAIL")
        or "gvf@example.com"
    )
    sess = _session(email)

    pmids = load_paywalled_pmids(paywalled_csv)
    items = list(pmids.items())
    if limit:
        items = items[:limit]

    results: List[AuditResult] = []
    counts = Counter()

    for idx, (pmid, reason) in enumerate(items, 1):
        doi = pubmed_get_doi(sess, pmid)
        time.sleep(sleep_s)

        unp = unpaywall_check(sess, doi, email)
        time.sleep(sleep_s)

        epmc = europepmc_check(sess, pmid)
        time.sleep(sleep_s)

        landing = doi_landing_check(sess, doi)
        time.sleep(sleep_s)

        cls = classify(unp, epmc, landing)
        counts[cls] += 1
        results.append(
            AuditResult(
                pmid=str(pmid),
                doi=doi,
                original_reason=reason,
                classification=cls,
                unpaywall=unp,
                europepmc=epmc,
                doi_landing=landing,
            )
        )

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump([r.__dict__ for r in results], f, indent=2)

    # minimal markdown report
    free_missed = [r for r in results if r.classification == "FREE_BUT_MISSED"]
    free_blocked = [r for r in results if r.classification == "FREE_BUT_BLOCKED"]

    md = [
        f"# Paywall Audit Report ({len(results)} papers)",
        "",
        "## Summary",
        "",
        "| Category | Count |",
        "|---|---:|",
    ]
    for k in ["TRULY_PAYWALLED", "LIKELY_PAYWALLED", "FREE_BUT_BLOCKED", "FREE_BUT_MISSED"]:
        md.append(f"| {k} | {counts.get(k, 0)} |")

    md += [
        "",
        f"## FREE_BUT_MISSED ({len(free_missed)})",
        "| PMID | DOI | Unpaywall URL |",
        "|---|---|---|",
    ]
    for r in free_missed:
        md.append(
            f"| {r.pmid} | {r.doi or ''} | {r.unpaywall.get('best_oa_url','') if r.unpaywall else ''} |"
        )

    md += [
        "",
        f"## FREE_BUT_BLOCKED ({len(free_blocked)})",
        "| PMID | DOI | Final URL |",
        "|---|---|---|",
    ]
    for r in free_blocked:
        md.append(f"| {r.pmid} | {r.doi or ''} | {r.doi_landing.get('final_url','') if r.doi_landing else ''} |")

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")

    return {
        "n": len(results),
        "counts": dict(counts),
        "out_json": str(out_json),
        "out_md": str(out_md),
    }
