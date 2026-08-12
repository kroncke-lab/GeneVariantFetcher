#!/usr/bin/env python3
"""Score a rollout tier twice: on all gold, and on gold the source could support.

A tier gate is supposed to measure the *reading protocol*. When a paper's body
never landed — expired EZproxy session, revoked publisher access, a supplement
behind a login — no prompt, model, or router could have produced the row, so
counting it as a miss measures the fetcher and reports it as an extraction
failure. This module separates the two and makes the excluded set explicit
rather than silent.

Why per gold ROW and not per paper
----------------------------------
Whole-paper exclusion is the obvious design and it is wrong here. SCN5A PMID
27566755 advertises one supplement that is not on disk, so the per-paper
predicate in ``source_stratified_metrics.py`` calls it SOURCE-INCOMPLETE — yet
the locked Tier 1 baseline extracted 51/51 of its gold rows with zero count
error, because the table it needed was folded into the body anyway. Dropping
that paper would discard a perfectly handled one. Conversely a paper can hold a
usable body and still have most of its gold rows sitting in an unfetched
supplement. Reachability is a property of the row, not the article.

The three classes
-----------------
``present_in_body``
    Some normalized form of the gold variant appears in the source on disk. The
    protocol had its chance. **Penalized.**
``body_absent``
    The paper has no usable body at all: nothing on disk, an abstract-only
    fallback, or a stub too short to be an article, and no other on-disk
    representation exposes the row. **Excluded and noted.**
``variant_absent_from_body``
    A usable body is on disk, no form of the variant appears in any searchable
    file, the notation is a plain substitution (the class where the probe is
    reliable), and there are no figure images that could be hiding it.
    **Excluded and noted.**
``notation_inconclusive``
    Absent from the searchable text, but the notation is an indel/frameshift/
    structural form the probe cannot reliably find (the source may write
    ``c.4389_4396delCCTCTTTA`` where gold writes ``P.N1380del``).
    **Penalized** — deliberately. A weak probe must never manufacture credit.
``absent_but_figures_present``
    Absent from the searchable text, but the paper's figure images are on disk
    and the pipeline has a figure-reading lane. The source was acquired; we
    simply cannot string-search a PNG. **Penalized** — this is a probe blind
    spot, not an acquisition failure. On Tier 1 this class holds the rows the
    figure layer really did read, which a text-only probe would have excluded.

The probe's false-exclusion rate is reported with every run and belongs in the
reading of any excluded stratum: it is the fraction of gold rows a paper layer
demonstrably extracted from the body that plain string search still fails to
find, restricted to substitutions. On the locked Tier 1 baseline it is a few
percent, driven by papers that write only the cDNA form of a protein variant.

Two integrity properties hold by construction:

* Classification reads the gold row and the on-disk source only. It never sees
  a prediction, a score, or a match, so it cannot be tuned to flatter a result.
* No scoring logic is reimplemented. The gold loader, the count resolver (which
  honours the adjudicated ``gold_v2_*`` columns via
  ``utils.gold_standard.authoritative_gold_count``), and the greedy one-to-one
  matcher are imported from ``run_eval`` itself, so a stratum cannot drift from
  the harness. When the run's population is exactly the tier's gold-bearing
  papers, the ALL-GOLD stratum is checked against the published ``overall``
  block: a matched-pair disagreement means this module drifted and the numbers
  are withheld (``pairing_ok: false``), while a value-only disagreement is
  reported as answer-key drift, since re-scoring locked predictions against an
  adjudicated key is the intended use rather than an error. A run that does not
  cover the tier gets ``applicable: false`` and its strata are labelled a
  partial subset instead of being presented as the tier result.

  python scripts/recall_audit/tier_source_reachability.py \
      --tier benchmarks/evaluation_tiers/tier1_gold_50.tsv \
      --run-dir benchmarks/codex_paper_eval/runs/20260726_fixed48_production \
      --corpus corpus
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.run_eval import (  # noqa: E402
    COUNT_FIELDS,
    DEFAULT_GOLD,
    extract_pdf_text,
    load_gold,
    matches,
    read_paper_manifest,
)
from scripts.recall_audit.source_stratified_metrics import (  # noqa: E402
    ABSTRACT_ONLY_MARKER,
    MIN_FULL_TEXT_CHARS,
)
from utils.env_utils import local_data_discovery_disabled  # noqa: E402
from utils.variant_normalizer import AA_MAP, AA_MAP_REVERSE, VariantNormalizer  # noqa: E402

# Layers that read the paper. A row matched only by clinvar/pubtator linkage was
# attributed from a database, so the body text is not evidence about it either
# way — it is excluded from the probe's false-exclusion calibration.
PAPER_LAYERS = frozenset(
    {"regex_table", "llm_table", "llm_text", "regex_text", "mixed", "figure"}
)

# A plain substitution: optional p., residue, position, residue. The only
# notation class where absence from the body text is trustworthy evidence that
# the row is unreachable — an indel or frameshift row is routinely written in a
# form string search will not find (gold ``P.N1380del`` vs body
# ``c.4389_4396delCCTCTTTA``).
#
# Both residues must be real amino-acid codes. A bare length check is not
# enough: ``P.F1617DEL`` has a three-character tail and would be accepted as a
# substitution, which is exactly the false exclusion this guards against.
_SUBSTITUTION_SHAPE = re.compile(
    r"^[Pp]?\.?(?P<ref>[A-Za-z]{1,3})(?P<pos>\d+)(?P<alt>[A-Za-z]{1,3}|\*)$"
)
_AA_CODES = frozenset(
    {a.upper() for a in AA_MAP} | {a.upper() for a in AA_MAP_REVERSE} | {"X", "*"}
)


def is_substitution(variant: str) -> bool:
    """True for a residue-to-residue substitution in any standard notation."""
    m = _SUBSTITUTION_SHAPE.match((variant or "").strip())
    if not m:
        return False
    return m.group("ref").upper() in _AA_CODES and m.group("alt").upper() in _AA_CODES


PENALIZED = frozenset(
    {"present_in_body", "notation_inconclusive", "absent_but_figures_present"}
)
EXCLUDED = frozenset({"body_absent", "variant_absent_from_body"})


def _squash(text: str) -> str:
    """Whitespace-free upper case, so ``p. Arg176 Trp`` matches ``p.Arg176Trp``."""
    return re.sub(r"\s+", "", text or "").upper()


# Text-bearing supplement formats. A row can live in a converted supplement that
# was never folded into FULL_CONTEXT, and searching only FULL_CONTEXT declared 46
# of Tier 1's 246 candidate exclusions unreachable when their variant is sitting
# in a sibling file on the same disk.
SEARCHABLE_SUFFIXES = frozenset(
    {".md", ".txt", ".csv", ".tsv", ".json", ".xml", ".html", ".htm"}
)
# One pathological supplement must not stall a sweep over hundreds of papers.
MAX_SEARCHABLE_BYTES = 64 * 1024 * 1024


@dataclass
class PaperSource:
    """Every searchable byte we hold for one paper, plus what we cannot search."""

    text: str = ""
    unusable_reason: str = ""
    body_chars: int = 0
    supplement_files: int = 0
    pdf_files: int = 0
    figure_files: int = 0


class BodyIndex:
    """All on-disk searchable text for a paper, read once.

    Deliberately wider than the article body. "Could the source we hold have
    supported this row" is a question about everything on disk for that paper --
    ``_FULL_CONTEXT.md``, the possibly-richer ``_CLEANED.md``, and any converted
    supplement -- not about whichever single file the extractor happened to pick.
    """

    def __init__(self, corpus: Path) -> None:
        self.corpus = corpus
        self._cache: dict[tuple[str, str], PaperSource] = {}

    def get(self, gene: str, pmid: str) -> PaperSource:
        key = (gene, pmid)
        if key in self._cache:
            return self._cache[key]
        self._cache[key] = self._load(gene, pmid)
        return self._cache[key]

    def _load(self, gene: str, pmid: str) -> PaperSource:
        paper = self.corpus / gene / pmid
        out = PaperSource()

        image_suffixes = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
        out.figure_files = sum(
            1
            for f in paper.rglob("*")
            if f.is_file() and f.suffix.lower() in image_suffixes
        )

        # A run can select CLEANED when FULL_CONTEXT is absent or weaker, and a
        # PDF/supplement/figure may contain a row even when the body is only a
        # stub. Inventory every representation before deciding that the source
        # is unusable; otherwise the diagnostic can credit a miss that the
        # benchmark reader actually had a chance to recover.
        body_candidates = [
            paper / f"{pmid}_FULL_CONTEXT.md",
            paper / f"{pmid}_CLEANED.md",
        ]
        parts: list[str] = []
        body_states: list[str] = []
        for candidate in body_candidates:
            if not candidate.is_file():
                continue
            try:
                raw = candidate.read_text(errors="replace")
            except OSError as exc:
                body_states.append(f"unreadable: {exc.strerror or exc}")
                continue
            if raw.lstrip().startswith(ABSTRACT_ONLY_MARKER):
                body_states.append("abstract-only fallback")
                continue
            if len(raw) < MIN_FULL_TEXT_CHARS:
                body_states.append("source too short to be a body")
            else:
                out.body_chars = max(out.body_chars, len(raw))
                parts.append(raw)

        budget = MAX_SEARCHABLE_BYTES - sum(len(part) for part in parts)
        supplements = paper / f"{pmid}_supplements"
        if supplements.is_dir():
            for f in sorted(supplements.rglob("*")):
                if budget <= 0:
                    break
                if not f.is_file() or f.suffix.lower() not in SEARCHABLE_SUFFIXES:
                    continue
                try:
                    if f.stat().st_size > budget:
                        continue
                    parts.append(f.read_text(errors="replace"))
                except OSError:
                    continue
                out.supplement_files += 1
                budget -= f.stat().st_size
        # Match the benchmark reader's PDF representation: pdftotext output is
        # source the model can actually receive, even when the Markdown body
        # does not contain the variant. Omitting this route could incorrectly
        # exclude a hard row and make the diagnostic stratum look better.
        pdfs = sorted(f for f in paper.rglob("*.pdf") if f.is_file())
        out.pdf_files = len(pdfs)
        if pdfs:
            pdf_text = extract_pdf_text(
                [str(path.resolve()) for path in pdfs], MAX_SEARCHABLE_BYTES
            )
            if pdf_text:
                parts.append(pdf_text)
        out.text = _squash("\n".join(parts))
        if not out.body_chars:
            out.unusable_reason = (
                body_states[0]
                if body_states
                else (
                    "no usable article body on disk"
                    if paper.exists()
                    else "no source on disk"
                )
            )
        return out


class FormIndex:
    """Normalized notation variants of a gold variant string."""

    def __init__(self) -> None:
        self._cache: dict[tuple[str, str], frozenset[str]] = {}

    def get(self, gene: str, variant: str) -> frozenset[str]:
        key = (gene, variant)
        if key in self._cache:
            return self._cache[key]
        forms = {variant.strip()}
        try:
            for value in VariantNormalizer(gene).get_all_forms(variant).values():
                if value:
                    forms.add(str(value))
        except Exception:  # noqa: BLE001 - a normalizer refusal is not fatal here
            pass
        for form in list(forms):
            forms.add(form.replace("p.", ""))
            if re.match(r"^[A-Za-z]{3}\d+", form):
                forms.add("p." + form)
        # Two characters cannot identify a variant and would match everywhere.
        result = frozenset(f for f in (_squash(x) for x in forms) if len(f) >= 3)
        self._cache[key] = result
        return result


@dataclass
class GoldRow:
    gene: str
    pmid: str
    variant: str
    values: dict
    reachability: str = "present_in_body"
    reason: str = ""


def classify_row(row: GoldRow, bodies: BodyIndex, forms: FormIndex) -> tuple[str, str]:
    """Classify one gold row's source reachability. Blind to any prediction."""
    src = bodies.get(row.gene, row.pmid)
    if any(form in src.text for form in forms.get(row.gene, row.variant)):
        return "present_in_body", "variant string present in on-disk source"
    if src.figure_files:
        # We hold figure images for this paper and cannot string-search them, and
        # the pipeline has a figure-reading lane. "Absent from the text" is
        # therefore not evidence the source could not support the row, so it is
        # not an acquisition failure and must stay penalized.
        return (
            "absent_but_figures_present",
            f"absent from searchable text, but {src.figure_files} unsearchable "
            "figure image(s) are on disk",
        )
    if src.unusable_reason:
        return "body_absent", src.unusable_reason
    if not is_substitution(row.variant):
        return (
            "notation_inconclusive",
            "absent from searchable text but notation class is not reliably searchable",
        )
    return (
        "variant_absent_from_body",
        "substitution absent from every searchable file on disk, and no figures",
    )


@dataclass
class Stratum:
    """Recall and count error for one population of gold rows."""

    label: str
    gold_rows: int = 0
    matched_rows: int = 0
    gold_pmids: set = field(default_factory=set)
    matched_pmids: set = field(default_factory=set)
    predicted: dict = field(default_factory=lambda: dict.fromkeys(COUNT_FIELDS, 0))
    # Gold rows that state a value for this field at all (non-empty cell).
    gold_valued: dict = field(default_factory=lambda: dict.fromkeys(COUNT_FIELDS, 0))
    abs_err: dict = field(default_factory=lambda: {f: 0.0 for f in COUNT_FIELDS})
    sq_err: dict = field(default_factory=lambda: {f: 0.0 for f in COUNT_FIELDS})
    n_err: dict = field(default_factory=lambda: dict.fromkeys(COUNT_FIELDS, 0))
    # Proper coverage-aware loss. A missing prediction is charged as predicting
    # zero, with a one-unit floor so abstaining on an explicit gold zero is not
    # free. MAE/RMSE above remain conditional on both sides supplying a value.
    gold_row_abs_err: dict = field(
        default_factory=lambda: {f: 0.0 for f in COUNT_FIELDS}
    )
    missing_predictions: dict = field(
        default_factory=lambda: dict.fromkeys(COUNT_FIELDS, 0)
    )

    def to_dict(self) -> dict:
        counts = {}
        for f in COUNT_FIELDS:
            n = self.n_err[f]
            counts[f] = {
                "gold_asserted": self.gold_rows,
                "gold_valued": self.gold_valued[f],
                "predicted": self.predicted[f],
                # Mirrors the harness: denominator is every gold row in the
                # stratum, whether or not its cell carries a value.
                "recall": (self.predicted[f] / self.gold_rows)
                if self.gold_rows
                else None,
                # The stricter, more meaningful denominator.
                "recall_vs_valued_gold": (self.predicted[f] / self.gold_valued[f])
                if self.gold_valued[f]
                else None,
                "n_compared": n,
                "missing_predictions": self.missing_predictions[f],
                "mae": (self.abs_err[f] / n) if n else None,
                "rmse": math.sqrt(self.sq_err[f] / n) if n else None,
                # Monotone in both goals: total error spread over every gold row
                # in the stratum, so leaving a count null is never free.
                "error_per_gold_row": (self.gold_row_abs_err[f] / self.gold_rows)
                if self.gold_rows
                else None,
            }
        return {
            "label": self.label,
            "gold_rows": self.gold_rows,
            "matched_rows": self.matched_rows,
            "row_recall": (self.matched_rows / self.gold_rows)
            if self.gold_rows
            else None,
            "gold_pmids": len(self.gold_pmids),
            "matched_pmids": len(self.matched_pmids),
            "counts": counts,
        }


def _num(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _predictions_by_paper(predictions: Path) -> dict:
    """``{(gene, pmid): [variant_record, ...]}`` in locked prediction order."""
    payload = json.loads(predictions.read_text(encoding="utf-8"))
    out: dict = {}
    for paper in payload.get("papers", []):
        key = ((paper.get("gene") or "").upper(), str(paper.get("pmid") or ""))
        out[key] = list(paper.get("variants", []))
    return out


def pair_paper(gene: str, predicted: list[dict], gold: list[dict]) -> tuple:
    """Reproduce ``run_eval.score_one``'s greedy one-to-one pairing exactly.

    Prediction-driven: each predicted row claims the first not-yet-claimed gold
    row it matches. Order matters, so the locked prediction order is preserved
    and the gold order comes from the harness's own loader. Returns
    ``(pairs, missed_gold_indices)``.
    """
    used: set[int] = set()
    pairs: list[tuple[dict, dict]] = []
    for p in predicted:
        hit = next(
            (
                i
                for i, g in enumerate(gold)
                if i not in used and matches(p["variant"], g["variant"], gene)
            ),
            None,
        )
        if hit is None:
            continue
        used.add(hit)
        pairs.append((p, gold[hit]))
    missed = [i for i in range(len(gold)) if i not in used]
    return pairs, missed


def score_tier(
    tier: list[tuple[str, str]],
    run_dir: Path,
    corpus: Path,
    gold_dir: Path = DEFAULT_GOLD,
) -> dict:
    """Classify every gold row in the tier and re-stratify the locked metrics."""
    bodies = BodyIndex(corpus)
    forms = FormIndex()
    predictions = _predictions_by_paper(run_dir / "predictions.json")
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))

    strata = {
        "all_gold": Stratum("ALL GOLD"),
        "source_reachable": Stratum("SOURCE-REACHABLE"),
        "excluded_unreachable": Stratum("EXCLUDED — SOURCE UNREACHABLE"),
    }
    by_class: dict = {}
    per_gene: dict = {}
    noted: list[dict] = []
    calibration = {"paper_matched_substitutions": 0, "probe_missed": 0}
    scored_papers = 0
    unscored: list[dict] = []

    for gene, pmid in tier:
        gold_rows = load_gold(gold_dir, gene, pmid)
        if not gold_rows:
            continue
        key = (gene, pmid)
        # Reachability is a property of the gold row and the corpus, so it is
        # computed for the whole tier whether or not a run has touched the paper
        # yet. That makes the exclusion list available *before* paying for an
        # extraction. Only the metric strata need run coverage.
        covered = key in predictions
        if not covered:
            unscored.append({"gene": gene, "pmid": pmid, "gold_rows": len(gold_rows)})
        else:
            scored_papers += 1
        gene_bucket = per_gene.setdefault(
            gene, {"gold_rows": 0, "reachable": 0, "excluded": 0, "classes": {}}
        )

        paired_by_index = {}
        if covered:
            pairs, _ = pair_paper(gene, predictions[key], gold_rows)
            for p, g in pairs:
                paired_by_index[id(g)] = p

        for gold_row in gold_rows:
            variant = (gold_row.get("variant") or "").strip()
            if not variant:
                continue
            row = GoldRow(gene=gene, pmid=pmid, variant=variant, values=gold_row)
            row.reachability, row.reason = classify_row(row, bodies, forms)
            by_class[row.reachability] = by_class.get(row.reachability, 0) + 1
            gene_bucket["gold_rows"] += 1
            gene_bucket["classes"][row.reachability] = (
                gene_bucket["classes"].get(row.reachability, 0) + 1
            )

            penalized = row.reachability in PENALIZED
            gene_bucket["reachable" if penalized else "excluded"] += 1

            pred = paired_by_index.get(id(gold_row))
            matched = pred is not None

            if not penalized:
                noted.append(
                    {
                        "gene": gene,
                        "pmid": pmid,
                        "variant": variant,
                        "class": row.reachability,
                        "reason": row.reason,
                        "gold": {f: gold_row.get(f) for f in COUNT_FIELDS},
                        "matched_anyway": matched,
                        "matched_by_layer": (pred or {}).get("source_layer"),
                        "run_covered": covered,
                    }
                )

            if not covered:
                continue

            # Calibration: only a row a PAPER layer pulled out of this body is
            # evidence about the probe. A clinvar/pubtator match says nothing
            # about what the body contains.
            if (
                matched
                and (pred.get("source_layer") or "unknown") in PAPER_LAYERS
                and is_substitution(variant)
            ):
                calibration["paper_matched_substitutions"] += 1
                if row.reachability == "variant_absent_from_body":
                    calibration["probe_missed"] += 1

            targets = [strata["all_gold"]]
            targets.append(
                strata["source_reachable"]
                if penalized
                else strata["excluded_unreachable"]
            )
            for stratum in targets:
                stratum.gold_rows += 1
                stratum.gold_pmids.add(f"{gene}:{pmid}")
                if matched:
                    stratum.matched_rows += 1
                    stratum.matched_pmids.add(f"{gene}:{pmid}")
                for f in COUNT_FIELDS:
                    if gold_row.get(f) is not None:
                        stratum.gold_valued[f] += 1

            for f in COUNT_FIELDS:
                gold_value = gold_row.get(f)
                if gold_value is None:
                    continue
                got = pred.get(f) if matched else None
                # The harness counts an observation only when BOTH sides carry a
                # value; a predicted 0 is a value, a null is not.
                if got is None:
                    penalty = max(1, abs(gold_value))
                    for stratum in targets:
                        stratum.gold_row_abs_err[f] += penalty
                        stratum.missing_predictions[f] += 1
                    continue
                for stratum in targets:
                    stratum.predicted[f] += 1
                diff = got - gold_value
                for stratum in targets:
                    stratum.abs_err[f] += abs(diff)
                    stratum.gold_row_abs_err[f] += abs(diff)
                    stratum.sq_err[f] += diff * diff
                    stratum.n_err[f] += 1

    n = calibration["paper_matched_substitutions"]
    calibration["false_exclusion_rate"] = calibration["probe_missed"] / n if n else None

    gold_papers = scored_papers + len(unscored)
    outside = sorted(f"{g}:{p}" for (g, p) in predictions if (g, p) not in set(tier))
    coverage = {
        "tier_attempts": len(tier),
        "tier_attempts_with_gold": gold_papers,
        "scored_papers": scored_papers,
        "papers_in_tier_without_run_coverage": unscored,
        "run_papers_outside_tier": outside,
        # The published `overall` block describes the run's own population. It
        # is only a valid reference when that population is exactly the tier's
        # gold-covered papers; otherwise a mismatch says nothing about this
        # module's fidelity and the strata are a partial subset, not the tier.
        "run_matches_tier": not unscored and not outside,
        # Distinct from the above: the strata cover every gold-bearing paper the
        # tier asks for. A run that additionally holds papers outside the tier
        # (the BRCA2 arm's 8-paper run against the active 2-paper cohort) still
        # produces complete tier strata; only the published overall block, which
        # describes the wider population, stops being a valid reference.
        "tier_fully_covered": not unscored,
    }

    payload = {
        "coverage": coverage,
        "reachability_classes": by_class,
        "per_gene": per_gene,
        "probe_calibration": calibration,
        "strata": {k: v.to_dict() for k, v in strata.items()},
        "noted_exclusions": noted,
    }
    if coverage["run_matches_tier"]:
        payload["reproduction"] = _check_reproduction(
            payload["strata"]["all_gold"], report
        )
    else:
        if coverage["tier_fully_covered"]:
            reason = (
                f"This run also holds {len(outside)} paper(s) outside the tier, "
                "so its published overall block describes a wider population "
                "and cannot verify this module's pairing. The strata below are "
                "complete for the tier."
            )
        else:
            reason = (
                f"This run covers {scored_papers} of the tier's {gold_papers} "
                "gold-bearing papers, so the strata below are a partial subset "
                "and must not be quoted as the tier result. The reachability "
                "classification is still complete."
            )
        payload["reproduction"] = {
            "applicable": False,
            "pairing_ok": None,
            "reason": reason,
            "pairing_checks": [],
            "answer_key_drift": [],
        }
    return payload


def _check_reproduction(all_gold: dict, report: dict) -> dict:
    """Check the ALL-GOLD stratum against the run's published overall block.

    Two very different things can disagree here and they must not share a
    verdict:

    *Pairing* — how many gold rows matched, and how many matched pairs carry a
    value on both sides. That is pure re-implementation fidelity. Any mismatch
    means this module's matcher drifted from the harness and nothing below can
    be quoted, so it is a hard gate.

    *Values* — ``gold_asserted`` and the error metrics. These legitimately move
    when the answer key is adjudicated after a run was scored: a new
    ``gold_v2_*`` row changes a gold value, and an
    ``excluded_duplicate_current_cohort`` status removes a row outright. The
    locked predictions did not change, so this is answer-key drift, not a
    scoring defect, and re-scoring against today's key is the *point* rather
    than an error. It is reported, not failed.
    """
    published = (report.get("overall") or {}).get("count") or {}
    pairing_checks, value_drift = [], []
    pairing_ok = True

    published_tp = (report.get("overall") or {}).get("tp")
    if published_tp is not None:
        match = published_tp == all_gold["matched_rows"]
        pairing_ok = pairing_ok and match
        pairing_checks.append(
            {
                "metric": "matched_rows",
                "published": published_tp,
                "recomputed": all_gold["matched_rows"],
                "match": match,
            }
        )

    for f in COUNT_FIELDS:
        want = published.get(f) or {}
        got = all_gold["counts"][f]
        a, b = want.get("predicted"), got.get("predicted")
        match = a == b
        pairing_ok = pairing_ok and match
        pairing_checks.append(
            {
                "metric": f"{f}.predicted",
                "published": a,
                "recomputed": b,
                "match": match,
            }
        )
        for metric in ("gold_asserted", "mae"):
            a, b = want.get(metric), got.get(metric)
            same = (a is None and b is None) or (
                isinstance(a, (int, float))
                and isinstance(b, (int, float))
                and abs(a - b) < 1e-9
            )
            if not same:
                value_drift.append(
                    {"metric": f"{f}.{metric}", "as_scored": a, "today": b}
                )

    return {
        "pairing_ok": pairing_ok,
        "pairing_checks": pairing_checks,
        "answer_key_drift": value_drift,
    }


def _pct(value) -> str:
    return "n/a" if value is None else f"{value * 100:.1f}%"


def _fmt(value) -> str:
    return "n/a" if value is None else f"{value:.3f}"


def render(tier_name: str, payload: dict) -> str:
    lines = [f"# Source-reachability stratified metrics — {tier_name}\n"]
    rep = payload["reproduction"]
    partial = rep.get("applicable") is False
    if partial:
        lines.append(f"> **Partial run coverage.** {rep['reason']}\n")
    elif not rep["pairing_ok"]:
        lines.append(
            "> **WITHHELD.** This module's pairing does not reproduce the run's\n"
            "> published matched-pair counts, so it has drifted from the harness\n"
            "> matcher and nothing below can be quoted. Failing checks:\n"
        )
        for c in rep["pairing_checks"]:
            if not c["match"]:
                lines.append(
                    f"> - {c['metric']}: published {c['published']} "
                    f"vs recomputed {c['recomputed']}"
                )
        lines.append("")
    if rep["answer_key_drift"]:
        lines.append(
            "> **Answer-key drift.** Pairing reproduces exactly, but the gold\n"
            "> standard has been adjudicated since this run was scored. The\n"
            "> predictions are unchanged; the numbers below are the locked\n"
            "> predictions re-scored against *today's* key, which is the\n"
            "> comparable figure. As-scored vs today:\n"
        )
        for d in rep["answer_key_drift"]:
            lines.append(f"> - {d['metric']}: {d['as_scored']} → **{d['today']}**")
        lines.append("")

    cal = payload["probe_calibration"]
    lines.append(
        "ALL GOLD remains the primary turnkey acceptance number. SOURCE-REACHABLE\n"
        "is a secondary reader diagnostic: gold rows whose variant is\n"
        "present in the source on disk, plus rows whose notation the probe cannot\n"
        "reliably search (kept penalized on purpose). EXCLUDED rows had no usable\n"
        "body, or a usable body that demonstrably does not contain them — nothing\n"
        "the reader could have done. The calibration sample is conditional on\n"
        "rows a paper layer matched and does not bound false exclusions among\n"
        "missed rows, so it cannot promote this diagnostic to a gate. Observed\n"
        "probe false-exclusion rate on matched rows: "
        f"{_pct(cal['false_exclusion_rate'])} "
        f"({cal['probe_missed']}/{cal['paper_matched_substitutions']}).\n"
    )

    lines.append("\n## Reachability of gold rows\n")
    lines.append("| class | gold rows | disposition |")
    lines.append("|---|---:|---|")
    for name, count in sorted(
        payload["reachability_classes"].items(), key=lambda kv: -kv[1]
    ):
        lines.append(
            f"| {name} | {count} | "
            f"{'penalized' if name in PENALIZED else 'excluded and noted'} |"
        )

    cov = payload["coverage"]
    heading = "Metrics by stratum"
    if partial and not cov["tier_fully_covered"]:
        heading += (
            f" — PARTIAL SUBSET ({cov['scored_papers']} of "
            f"{cov['tier_attempts_with_gold']} gold-bearing tier papers)"
        )
    elif partial:
        heading += " — complete for the tier, pairing unverified"
    lines.append(f"\n## {heading}\n")
    lines.append(
        "| stratum | rows | row recall | carriers recall / MAE "
        "| affected recall / MAE | unaffected recall / MAE |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|")
    for key in ("all_gold", "source_reachable", "excluded_unreachable"):
        s = payload["strata"][key]
        c = s["counts"]
        lines.append(
            f"| {s['label']} | {s['matched_rows']}/{s['gold_rows']} "
            f"| {_pct(s['row_recall'])} "
            f"| {_pct(c['carriers']['recall'])} / {_fmt(c['carriers']['mae'])} "
            f"| {_pct(c['affected']['recall'])} / {_fmt(c['affected']['mae'])} "
            f"| {_pct(c['unaffected']['recall'])} / {_fmt(c['unaffected']['mae'])} |"
        )

    lines.append(
        "\n`error_per_gold_row` is the form that is monotone in both goals — "
        "total absolute error divided by every gold row in the stratum, so a "
        "null count is never cheaper than a wrong one:\n"
    )
    lines.append("| stratum | carriers | affected | unaffected |")
    lines.append("|---|---:|---:|---:|")
    for key in ("all_gold", "source_reachable"):
        s = payload["strata"][key]
        c = s["counts"]
        lines.append(
            f"| {s['label']} | {_fmt(c['carriers']['error_per_gold_row'])} "
            f"| {_fmt(c['affected']['error_per_gold_row'])} "
            f"| {_fmt(c['unaffected']['error_per_gold_row'])} |"
        )

    noted = payload["noted_exclusions"]
    lines.append(f"\n## Noted exclusions ({len(noted)} gold rows)\n")
    if not noted:
        lines.append("None — every gold row in this tier is source-reachable.\n")
    else:
        lines.append("| gene | pmid | variant | class | matched anyway |")
        lines.append("|---|---|---|---|---|")
        for row in noted:
            lines.append(
                f"| {row['gene']} | {row['pmid']} | {row['variant']} "
                f"| {row['class']} | {'yes' if row['matched_anyway'] else 'no'} |"
            )
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--tier", type=Path, required=True, help="gene<TAB>pmid manifest")
    ap.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="codex_paper_eval run dir holding predictions.json and report.json",
    )
    ap.add_argument(
        "--corpus",
        type=Path,
        default=None,
        help="Source corpus root. Defaults to <repo>/corpus unless implicit "
        "local-data discovery is disabled, in which case it is required.",
    )
    ap.add_argument("--gold-dir", type=Path, default=DEFAULT_GOLD)
    ap.add_argument(
        "--genes",
        default="",
        help="Comma-separated subset of the tier. Tier 1 deliberately holds two "
        "arms with different gold provenance and different production runs, so "
        "score the cardiac four and BRCA2 separately rather than pooling them.",
    )
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args()

    corpus = args.corpus
    if corpus is None:
        if local_data_discovery_disabled():
            print(
                "--corpus is required when GVF_DISABLE_LOCAL_DATA is set; "
                "refusing to guess <repo>/corpus."
            )
            return 2
        corpus = REPO / "corpus"
    if not corpus.is_dir():
        print(
            f"Corpus not reachable at {corpus}. Mount the external volume "
            "(see CLAUDE.md 'Operating Shape') or pass --corpus."
        )
        return 2

    tier = read_paper_manifest(args.tier)
    wanted = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
    if wanted:
        tier = [(gene, pmid) for gene, pmid in tier if gene in wanted]
        if not tier:
            print(f"No {args.tier} entries for genes: {sorted(wanted)}")
            return 2
    payload = score_tier(tier, args.run_dir, corpus, gold_dir=args.gold_dir)
    payload["inputs"] = {
        "tier": str(args.tier),
        "genes": sorted(wanted) or "all",
        "run_dir": str(args.run_dir),
        "corpus": str(corpus),
        "gold_dir": str(args.gold_dir),
    }

    name = args.tier.stem + (f" ({','.join(sorted(wanted))})" if wanted else "")
    md = render(name, payload)
    print(md)
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"Wrote {args.out_json}")
    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(md, encoding="utf-8")
        print(f"Wrote {args.out_md}")
    return 0 if payload["reproduction"]["pairing_ok"] is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
