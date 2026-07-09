#!/usr/bin/env python3
"""Acquisition expected-value (EV) prototype — abstract-only yield prediction.

Goal (see docs/VARIANT_BROWSER_INTEGRATION.md discussion / TASKS acquisition lever):
for a paper we have NOT downloaded, predict from its *abstract alone* whether it is
worth acquiring, i.e. how many phenotyped carriers and distinct genetic variants a
full extraction would likely yield. Rank the un-downloaded tail by that expected
value so manual acquisition effort goes to the highest-payoff papers first.

This module is the deterministic v1 + its honest evaluation. It does NOT fetch full
text and it does NOT call an LLM. Every feature is computed from the PubMed abstract
+ title + metadata, so the eval faithfully simulates the deployment case (abstract is
all we have). The score is a FIXED, pre-registered formula (no fitting), so a win on
held-out gold cannot be overfitting.

Evaluation target = the *true full-paper* gold counts (which include supplement
tables the abstract never shows). We measure, on the cardiac ion-channel gold
standards (KCNH2, KCNQ1, SCN5A, RYR2):

  * Spearman rank correlation of the score (and each raw signal) vs true carriers
    and true unique variants;
  * a yield-capture curve: if we acquire the top-K% of papers by score, what
    fraction of all gold carriers / variants do we capture, vs a random baseline
    and an oracle (sort-by-truth) ceiling.

"Not just noise" = the score's capture materially beats random and Spearman is
positive and significant. The regex signal set mirrors pipeline/extraction_priority.py
(kept inline so this evaluator is standalone and runnable without the pipeline .env).

Usage:
  python scripts/acquisition_ev/predict_yield.py --email you@example.com
  python scripts/acquisition_ev/predict_yield.py --genes KCNH2,SCN5A --refresh
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
GOLD_DIR = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
CACHE = HERE / ".cache" / "abstracts.json"
RESULTS = HERE / "eval_output"

CARDIAC_GENES = ["KCNH2", "KCNQ1", "SCN5A", "RYR2"]

# Aliases so gene-mention counting is realistic (mirrors config alias intent).
GENE_ALIASES = {
    "KCNH2": ["KCNH2", "hERG", "HERG", "LQT2", "Kv11.1"],
    "KCNQ1": ["KCNQ1", "KvLQT1", "KCNA9", "LQT1", "Kv7.1"],
    "SCN5A": ["SCN5A", "Nav1.5", "NaV1.5", "hH1"],
    "RYR2": ["RYR2", "RyR2", "ryanodine receptor 2", "cardiac ryanodine"],
}

# --- Regex signal set (mirrors pipeline/extraction_priority.py) ---------------
VARIANT_PATTERNS = (
    re.compile(r"\bc\.\d+[A-Za-z0-9_>*+\-?]+", re.IGNORECASE),
    re.compile(r"\bp\.[A-Za-z]{1,3}\d+[A-Za-z*]{1,3}\b", re.IGNORECASE),
    re.compile(r"\brs\d{4,}\b", re.IGNORECASE),
    re.compile(r"\bIVS\d+[+\-]\d+[A-Z]>[A-Z]\b", re.IGNORECASE),
    re.compile(r"\b\d{2,6}(?:del|ins|dup)[A-Za-z0-9]*\b", re.IGNORECASE),
    # single-letter protein change e.g. R190Q, G628S (common in abstracts)
    re.compile(r"\b[ACDEFGHIKLMNPQRSTVWY]\d{2,4}[ACDEFGHIKLMNPQRSTVWYX*]\b"),
)
CARRIER_COUNT_RE = re.compile(
    r"\b(?:carrier|carriers|proband|probands|patient|patients|individual|individuals|"
    r"family|families|relative|relatives|subject|subjects|case|cases|mutation[- ]?positive)\b",
    re.IGNORECASE,
)
ORIGINAL_DATA_RE = re.compile(
    r"\b(?:mutation|variant|genotyp\w*|phenotyp\w*|sequenc\w*|screen\w*|cohort|"
    r"pedigree|proband|clinical|affected|unaffected|penetrance|carrier)\b",
    re.IGNORECASE,
)
REVIEW_TITLE_RE = re.compile(
    r"\b(?:review|meta[- ]analysis|editorial|comment|guideline|consensus|"
    r"perspective|erratum|correction|protocol)\b",
    re.IGNORECASE,
)
# cohort-size: an integer adjacent to a people/family noun, or n = NNN.
COHORT_NEAR_RE = re.compile(
    r"\b(\d{1,6})\s+(?:unrelated\s+|consecutive\s+|index\s+|additional\s+)?"
    r"(?:patients|probands|individuals|subjects|cases|carriers|relatives|families|"
    r"members|participants|children|men|women|LQTS|CPVT|Brugada)\b",
    re.IGNORECASE,
)
N_EQUALS_RE = re.compile(r"\b[nN]\s*=\s*(\d{1,6})\b")
COHORT_CAP = 50000  # ignore absurd numbers (years, base pairs, etc.)


def alias_regex(gene: str) -> re.Pattern:
    parts = [re.escape(a) for a in GENE_ALIASES.get(gene, [gene])]
    return re.compile(
        r"(?<![A-Za-z0-9])(?:" + "|".join(parts) + r")(?![A-Za-z0-9])", re.IGNORECASE
    )


# --- Gold truth ---------------------------------------------------------------
def load_gold(genes: list[str]) -> dict[str, dict[str, dict]]:
    """gene -> pmid -> {variants, carriers, affected}. True full-paper counts."""
    out: dict[str, dict[str, dict]] = {}
    for g in genes:
        path = GOLD_DIR / f"{g}_recall_input.csv"
        if not path.exists():
            print(f"  WARN: no gold file for {g} at {path}; skipping")
            continue
        per = defaultdict(lambda: {"vset": set(), "carriers": 0, "affected": 0})
        with open(path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                pmid = (row.get("pmid") or "").strip()
                var = (row.get("variant") or "").strip()
                if not pmid:
                    continue
                if var:
                    per[pmid]["vset"].add(var)
                for k in ("carriers", "affected"):
                    try:
                        per[pmid][k] += int(float(row.get(k) or 0))
                    except (TypeError, ValueError):
                        pass
        out[g] = {
            pmid: {
                "variants": len(d["vset"]),
                "carriers": d["carriers"],
                "affected": d["affected"],
            }
            for pmid, d in per.items()
        }
    return out


# --- Abstract acquisition (PubMed efetch, cached) -----------------------------
def _clean(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", text or "")).strip()


def fetch_abstracts(
    pmids: list[str], email: str, refresh: bool = False
) -> dict[str, dict]:
    cache: dict[str, dict] = {}
    if CACHE.exists() and not refresh:
        cache = json.loads(CACHE.read_text(encoding="utf-8"))
    todo = [p for p in pmids if p not in cache]
    if todo:
        print(f"  fetching {len(todo)} abstracts from PubMed ({len(cache)} cached)...")
    for i in range(0, len(todo), 200):
        batch = todo[i : i + 200]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "abstract",
            "retmode": "xml",
            "email": email,
        }
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            + urllib.parse.urlencode(params)
        )
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as response:
                    data = response.read()
                root = ET.fromstring(data)
                break
            except Exception as exc:  # noqa: BLE001
                if attempt == 2:
                    print(f"  WARN batch {i}: {exc}")
                    root = None
                time.sleep(1.5 * (attempt + 1))
        if root is None:
            continue
        for art in root.findall(".//PubmedArticle"):
            pmid = (art.findtext(".//PMID") or "").strip()
            if not pmid:
                continue
            title = _clean(
                " ".join(t.text or "" for t in art.findall(".//ArticleTitle"))
            )
            abst = _clean(
                " ".join(
                    _clean(ET.tostring(a, "unicode"))
                    for a in art.findall(".//Abstract/AbstractText")
                )
            )
            journal = _clean(art.findtext(".//Journal/Title") or "")
            year = (
                art.findtext(".//JournalIssue/PubDate/Year")
                or art.findtext(".//PubDate/Year")
                or ""
            ).strip()
            ptypes = [(_clean(p.text or "")) for p in art.findall(".//PublicationType")]
            cache[pmid] = {
                "title": title,
                "abstract": abst,
                "journal": journal,
                "year": year,
                "pubtypes": ptypes,
            }
        time.sleep(0.4)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache), encoding="utf-8")
    return cache


# --- Features + score (FIXED formula, no fitting) -----------------------------
def _count(patterns, text: str) -> int:
    if isinstance(patterns, (list, tuple)):
        return sum(len(p.findall(text)) for p in patterns)
    return len(patterns.findall(text))


def compute_features(gene: str, rec: dict) -> dict:
    title = rec.get("title", "") or ""
    abstract = rec.get("abstract", "") or ""
    text = f"{title}. {abstract}"
    ptypes = " ".join(rec.get("pubtypes", []))

    cohort_vals = [int(m) for m in COHORT_NEAR_RE.findall(text)] + [
        int(m) for m in N_EQUALS_RE.findall(text)
    ]
    cohort_vals = [v for v in cohort_vals if 0 < v <= COHORT_CAP]
    cohort = max(cohort_vals) if cohort_vals else 0

    review = bool(REVIEW_TITLE_RE.search(title)) or ("Review" in ptypes)
    return {
        "has_abstract": int(bool(abstract)),
        "gene_hits": _count(alias_regex(gene), text),
        "variant_mentions": _count(VARIANT_PATTERNS, text),
        "carrier_mentions": _count(CARRIER_COUNT_RE, text),
        "original_data": _count(ORIGINAL_DATA_RE, text),
        "cohort_size": cohort,
        "review": int(review),
        "abstract_len": len(abstract),
    }


def score(feat: dict) -> dict:
    """Two exposed components + combined EV. Fixed, interpretable, unfit."""
    # p_relevant: does it contain phenotyped individuals WITH variant identities?
    z = (
        0.55 * feat["original_data"]
        + 0.80 * feat["carrier_mentions"]
        + 0.60 * min(feat["variant_mentions"], 12)
        + 0.50 * min(feat["gene_hits"], 6)
        - 2.5 * feat["review"]
        - 3.0
    )
    p_relevant = 1.0 / (1.0 + math.exp(-0.35 * z)) if feat["has_abstract"] else 0.05

    # Volume estimates (floors — abstracts under-report supplement content).
    est_carriers = float(
        feat["cohort_size"] if feat["cohort_size"] else feat["carrier_mentions"]
    )
    est_variants = float(feat["variant_mentions"])

    # Combined ranking EV (log-additive; robust to outliers).
    gate = 0.3 if feat["review"] else 1.0
    ev = gate * (
        math.log1p(feat["variant_mentions"])
        + math.log1p(feat["carrier_mentions"])
        + 1.2 * math.log1p(feat["cohort_size"])
        + 0.5 * math.log1p(feat["original_data"])
    )
    ev *= p_relevant
    return {
        "p_relevant": round(p_relevant, 4),
        "est_variants": est_variants,
        "est_carriers": est_carriers,
        "ev_score": round(ev, 4),
    }


# --- Evaluation ---------------------------------------------------------------
def _ranks(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> tuple[float, float]:
    n = len(xs)
    if n < 3:
        return 0.0, 1.0
    rx, ry = _ranks(xs), _ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    vy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if vx == 0 or vy == 0:
        return 0.0, 1.0
    rho = cov / (vx * vy)
    # two-sided p via t approximation
    denom = 1 - rho * rho
    if denom <= 0.0:  # rho at +/-1, or float underflow near it: p ~ 0
        return rho, 0.0
    t = rho * math.sqrt((n - 2) / denom)
    p = 2 * (1 - NormalDist().cdf(abs(t)))
    return rho, p


def capture_at(items: list[dict], score_key: str, truth_key: str, fracs) -> dict:
    """Fraction of total truth captured in the top-frac papers by score_key."""
    total = sum(it[truth_key] for it in items) or 1
    ordered = sorted(items, key=lambda it: it[score_key], reverse=True)
    n = len(ordered)
    out = {}
    for f in fracs:
        k = max(1, round(f * n))
        out[f] = sum(it[truth_key] for it in ordered[:k]) / total
    return out


def evaluate(items: list[dict], label: str) -> dict:
    fracs = [0.1, 0.2, 0.3, 0.5]
    res = {"label": label, "n": len(items)}
    for truth in ("carriers", "variants"):
        rho_ev, p_ev = spearman(
            [it["ev_score"] for it in items], [it[truth] for it in items]
        )
        cap_ev = capture_at(items, "ev_score", truth, fracs)
        cap_or = capture_at(items, truth, truth, fracs)  # oracle
        cap_len = capture_at(items, "abstract_len", truth, fracs)  # naive baseline
        res[truth] = {
            "spearman_ev": round(rho_ev, 3),
            "spearman_p": f"{p_ev:.1e}",
            "capture_ev": {f"{int(f * 100)}%": round(cap_ev[f], 3) for f in fracs},
            "capture_random": {f"{int(f * 100)}%": f for f in fracs},
            "capture_oracle": {f"{int(f * 100)}%": round(cap_or[f], 3) for f in fracs},
            "capture_len_baseline": {
                f"{int(f * 100)}%": round(cap_len[f], 3) for f in fracs
            },
            "lift_at_20pct": round(cap_ev[0.2] / 0.2, 2),
        }
    # raw-signal diagnostics vs carriers
    res["signal_spearman_vs_carriers"] = {
        s: round(
            spearman([it["feat"][s] for it in items], [it["carriers"] for it in items])[
                0
            ],
            3,
        )
        for s in (
            "gene_hits",
            "variant_mentions",
            "carrier_mentions",
            "original_data",
            "cohort_size",
            "review",
        )
    }
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--email",
        default="brett.kroncke@gmail.com",
        help="NCBI contact email for efetch",
    )
    ap.add_argument("--genes", default=",".join(CARDIAC_GENES))
    ap.add_argument("--refresh", action="store_true", help="ignore abstract cache")
    args = ap.parse_args()
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]

    print(f"Loading gold truth for {genes} ...")
    gold = load_gold(genes)
    all_pmids = sorted({p for g in genes for p in gold[g]})
    print(f"  {len(all_pmids)} gold PMIDs")

    abstracts = fetch_abstracts(all_pmids, args.email, refresh=args.refresh)

    items: list[dict] = []
    no_abs = 0
    for g in genes:
        for pmid, truth in gold[g].items():
            rec = abstracts.get(pmid, {})
            if not rec.get("abstract"):
                no_abs += 1
            feat = compute_features(g, rec)
            sc = score(feat)
            items.append(
                {
                    "gene": g,
                    "pmid": pmid,
                    "carriers": truth["carriers"],
                    "variants": truth["variants"],
                    "abstract_len": feat["abstract_len"],
                    "feat": feat,
                    **sc,
                }
            )
    print(f"  scored {len(items)} papers ({no_abs} had no abstract available)")

    RESULTS.mkdir(parents=True, exist_ok=True)
    # per-paper CSV
    with open(RESULTS / "per_paper_scores.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(
            [
                "gene",
                "pmid",
                "true_carriers",
                "true_variants",
                "p_relevant",
                "est_carriers",
                "est_variants",
                "ev_score",
                "cohort_size",
                "variant_mentions",
                "carrier_mentions",
                "review",
            ]
        )
        for it in items:
            w.writerow(
                [
                    it["gene"],
                    it["pmid"],
                    it["carriers"],
                    it["variants"],
                    it["p_relevant"],
                    it["est_carriers"],
                    it["est_variants"],
                    it["ev_score"],
                    it["feat"]["cohort_size"],
                    it["feat"]["variant_mentions"],
                    it["feat"]["carrier_mentions"],
                    it["feat"]["review"],
                ]
            )

    metrics = {"pooled": evaluate(items, "POOLED (cardiac)")}
    for g in genes:
        metrics[g] = evaluate([it for it in items if it["gene"] == g], g)
    (RESULTS / "metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )

    # console summary
    print("\n" + "=" * 74)
    print("ACQUISITION-EV EVALUATION vs true cardiac gold counts")
    print("=" * 74)
    for key in ["pooled"] + genes:
        m = metrics[key]
        print(f"\n[{m['label']}]  n={m['n']}")
        for truth in ("carriers", "variants"):
            t = m[truth]
            print(
                f"  {truth:9} Spearman(EV) = {t['spearman_ev']:+.3f} (p={t['spearman_p']})  lift@20%={t['lift_at_20pct']}x"
            )
            print(
                f"            capture@20%: EV={t['capture_ev']['20%']:.2f}  random=0.20  oracle={t['capture_oracle']['20%']:.2f}  len-baseline={t['capture_len_baseline']['20%']:.2f}"
            )
        print(f"  raw-signal Spearman vs carriers: {m['signal_spearman_vs_carriers']}")

    p = metrics["pooled"]
    verdict_ok = (
        p["carriers"]["spearman_ev"] >= 0.25
        and p["variants"]["spearman_ev"] >= 0.25
        and p["carriers"]["lift_at_20pct"] >= 1.5
        and p["variants"]["lift_at_20pct"] >= 1.5
    )
    print("\n" + "=" * 74)
    print(
        f"VERDICT: {'SIGNAL — beats noise, worth pushing' if verdict_ok else 'WEAK — looks like noise, do NOT push as-is'}"
    )
    print("=" * 74)
    print(f"Artifacts: {RESULTS / 'metrics.json'} , {RESULTS / 'per_paper_scores.csv'}")
    return 0 if verdict_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
