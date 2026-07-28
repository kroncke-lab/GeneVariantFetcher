"""Targeted carrier-count recovery for variants already found in a paper.

The count stack in this repo is entirely *subtractive*: :mod:`carrier_guard`
NULLs implausible counts, :mod:`count_outlier_guard` clears statistical
outliers, and :mod:`count_classifier` refuses counts whose provenance is not
per-variant. Nothing fills a count that was never emitted — and omission, not
error, is the dominant count deficit.

Measured on the locked 48-paper cardiac eval set (2026-07-26,
``benchmarks/codex_paper_eval/runs/20260726_fixed48_production``): of the 789
gold variants production correctly identified, it left ``total_carriers_observed``
NULL on 461 of them (58.4%). A single-model protocol reading the same sources
left only 11.3% null. When production *does* commit a count it is accurate
(85.5% precision on count-bearing predictions), so the gap is not an accuracy
problem that a reviewer could catch — it is a blank that has to be filled.

This module fills those blanks and nothing else:

* It only asks about variants **already stored** for the paper. It is a fill,
  not an extraction — a returned variant we did not ask about is discarded.
* It only writes into NULL slots. An existing count is never overwritten, so
  the pass cannot regress a count the extractor already got right.
* Every accepted count must be grounded in a short verbatim quote from the
  source that locally binds the requested value to the requested variant.
  Detached table blocks and numbers that only occur inside HGVS notation are
  dropped, not stored with a caveat.
* Writes are logged to ``count_recovery_log`` with the quote and model, so the
  pass is auditable and fully reversible.

Gold-free by construction: it compares the database against the paper, never
against a curated answer, so it works on new gene-disease pairs.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import sqlite3
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence

logger = logging.getLogger(__name__)

COUNT_RECOVERY_VERSION = "count-recovery-v2"

#: Fields this pass may fill, mapped to their ``penetrance_data`` column.
RECOVERABLE_FIELDS = {
    "carriers": "total_carriers_observed",
    "affected": "affected_count",
    "unaffected": "unaffected_count",
}

#: Upper bound on a plausible per-variant count, matching
#: :func:`pipeline.carrier_guard.apply_carrier_guard`'s default threshold so
#: this pass cannot introduce a row the guard would immediately quarantine.
PLAUSIBLE_COUNT_CEILING = 100_000
GROUNDING_QUOTE_MAX_CHARS = 2_000

#: Spelled-out numbers accepted as grounding for a numeric count. Without this
#: a quote reading "the variant was seen once" cannot support carriers=1, which
#: is a documented failure mode of strict digit-only grounding.
_WORD_NUMBERS = {
    "no one": 0,
    "neither": 0,
    "no": 0,
    "none": 0,
    "zero": 0,
    "a single": 1,
    "one": 1,
    "once": 1,
    "two": 2,
    "twice": 2,
    "both": 2,
    "three": 3,
    "thrice": 3,
    "four": 4,
    "five": 5,
    "six": 6,
    "seven": 7,
    "eight": 8,
    "nine": 9,
    "ten": 10,
    "eleven": 11,
    "twelve": 12,
}

_WS_RE = re.compile(r"\s+")
_DIGIT_COUNT_RE = re.compile(r"(?<![\w.])\d[\d,]*(?![\w.])")
_COMPOUND_NUMBER_TOKENS = frozenset(
    {
        *(_WORD_NUMBERS.keys()),
        "thirteen",
        "fourteen",
        "fifteen",
        "sixteen",
        "seventeen",
        "eighteen",
        "nineteen",
        "twenty",
        "thirty",
        "forty",
        "fifty",
        "sixty",
        "seventy",
        "eighty",
        "ninety",
        "hundred",
        "thousand",
        "million",
    }
)
_FIELD_CONTEXT_RE = {
    "carriers": re.compile(
        r"\b(?:carrier|carriers|carried|carries|harbou?rs?|harbou?red|"
        r"patients?|individuals?|subjects?|probands?|participants?|siblings?|cases?|"
        r"identified|observed|seen|found)\b",
        re.IGNORECASE,
    ),
    "affected": re.compile(
        r"\b(?:affected|symptomatic|phenotype[- ]positive|events?|cases?|disease)\b",
        re.IGNORECASE,
    ),
    "unaffected": re.compile(
        r"\b(?:unaffected|asymptomatic|phenotype[- ]negative|controls?)\b",
        re.IGNORECASE,
    ),
}
#: Population-scale / denominator vocabulary. Strong enough to disqualify a
#: number anywhere in its ±35-character neighbourhood, and applied to BOTH the
#: prose and the table branch: an explicitly labelled ``allele count`` /
#: ``gnomAD AC`` / ``families screened`` column is precisely the PMID 33013630
#: annotation-table failure class, and the table branch previously applied no
#: role check and no negative-context check at all.
_NON_COUNT_CONTEXT_RE = re.compile(
    r"\b(?:"
    r"allele (?:number|numbers|count|counts|frequency|frequencies)|"
    r"gnomad|exac|topmed|1000 ?genomes|esp6500|1000g|"
    r"minor allele|frequency|frequencies|percentile|prevalence|"
    r"sample size|cohort size|study (?:size|total|population|cohort)|"
    r"total (?:cohort|studied|screened|population)|"
    r"famil(?:y|ies)|pedigrees?|kindreds?|"
    r"screened|screening|enrolled|tested|genotyped|sequenced|referred|recruited"
    r")\b",
    re.IGNORECASE,
)
#: Measurement / annotation vocabulary. These words legitimately appear near a
#: real carrier count in prose ("in exon 5, 3 carriers"), so they only disqualify
#: a number when they sit *immediately* beside it or label its column.
_MEASURE_UNIT_RE = re.compile(
    r"\b(?:"
    r"age|ages|onset|years?|yrs?|months?|days?|follow[- ]?up|duration|"
    r"qtc|qt|jtc|msec|ms|bpm|beats|heart rate|mmhg|kg|cm|bmi|"
    r"exon|intron|codon|position|nucleotide|residue|chromosome|chr|"
    r"p[- ]?value|odds ratio|hazard ratio|confidence interval|"
    r"score|scores|cadd|cadd_?phred|sift|polyphen|revel|grantham"
    r")\b",
    re.IGNORECASE,
)
#: How far a measurement/annotation word reaches in prose. Deliberately tight.
_MEASURE_UNIT_WINDOW = 14
#: Table-column labels that are not count columns. Only safe to match inside a
#: header/label cell: ``\bAN\b`` would otherwise hit the English article "an",
#: and a prose sentence may legitimately cite "Table 2" beside a real count (the
#: cited number is itself a candidate integer, so the ≥2-candidate rule already
#: makes that quote fail closed).
_NON_COUNT_TABLE_LABEL_RE = re.compile(
    r"(?:^|[\s|(\[/])"
    r"(?:ac|an|af|maf|hom|het|alleles?|n_?alleles|"
    r"table|figure|fig|ref|refs|reference|references|page|supplementary)"
    r"(?:$|[\s|)\]/:.,])",
    re.IGNORECASE,
)
#: Unambiguous denominator markers: the number that follows is a cohort size.
_STRONG_DENOMINATOR_MARKER_RE = re.compile(
    r"(?:\bamong\b|\bamongst\b|\bout of\b|/)\s*(?:the\s+)?$", re.IGNORECASE
)
#: Ambiguous markers. "identified in 44 carriers" is not a denominator, but
#: "12 of 812" and "12 in 812" are — so these only count when another integer
#: (or an ``n=`` head) immediately precedes the marker.
_WEAK_DENOMINATOR_MARKER_RE = re.compile(
    r"(?:\bof\b|\bin\b|\bfrom\b|\bper\b|\bwithin\b|\bacross\b)\s*(?:the\s+)?$",
    re.IGNORECASE,
)
_N_EQUALS_TAIL_RE = re.compile(r"\bn\s*=\s*\d[\d,]*\s*$", re.IGNORECASE)
#: "<n> cases of <phenotype>" / "<n> patients with <disease>" — a disease-cohort
#: total. The integer describes the cohort, not carriers of the target variant.
_COHORT_HEAD_RE = re.compile(
    r"^\s*(?:cases?|patients?|individuals?|subjects?|probands?|participants?|"
    r"controls?|famil(?:y|ies)|family members?|relatives?)\s+"
    r"(?:of|with|who|having|diagnosed|referred|affected|screened|from)\b",
    re.IGNORECASE,
)
#: Explicit per-variant role labels. When a quote states more than one candidate
#: integer, only a value carrying one of these *between the neighbouring
#: candidates* can be trusted — otherwise the pass fails closed.
_EXPLICIT_ROLE_RE = {
    "carriers": re.compile(
        r"\b(?:carrier|carriers|carried|carries|harbou?rs?|harbou?red|"
        r"variant[- ]positive|genotype[- ]positive|mutation[- ]positive|"
        r"heterozygotes?|homozygotes?)\b",
        re.IGNORECASE,
    ),
    "affected": re.compile(
        r"\b(?:affected|symptomatic|phenotype[- ]positive)\b", re.IGNORECASE
    ),
    "unaffected": re.compile(
        r"\b(?:unaffected|asymptomatic|phenotype[- ]negative)\b", re.IGNORECASE
    ),
}

#: The structured count role a recovered value must carry to be written. The
#: model is asked to declare one; a non-per-variant declaration is a veto even
#: when the deterministic gate would otherwise pass, and the accepted role is
#: always the one the deterministic evidence proves.
PER_VARIANT_ROLES = {
    "carriers": "per_variant_carriers",
    "affected": "per_variant_affected",
    "unaffected": "per_variant_unaffected",
}
#: Roles the model may declare that disqualify the value outright.
NON_PER_VARIANT_ROLES = frozenset(
    {
        "cohort_total",
        "study_total",
        "denominator",
        "allele_count",
        "allele_number",
        "allele_frequency",
        "population_frequency",
        "family_count",
        "measurement",
        "other",
        "unknown",
    }
)
_GENE_PREFIXED_VARIANT_RE = re.compile(
    r"\b(?P<gene>[A-Z][A-Z0-9]{1,15})\s*[-:/]\s*"
    r"(?P<variant>(?:[pc]\.)?[A-Za-z0-9_()>*]+)"
)
_VARIANT_LIKE_CELL_RE = re.compile(
    r"(?:[pc]\.\(?)?(?:[A-Za-z]{1,3}\d+[A-Za-z*]{1,3}|"
    r"\d+[ACGT]>[ACGT]|\d+(?:_\d+)?(?:del|dup|ins))",
    re.IGNORECASE,
)
_TABLE_FIELD_LABEL_RE = {
    "carriers": re.compile(
        r"\b(?:carriers?|variant[- ]positive|genotype[- ]positive)\b",
        re.IGNORECASE,
    ),
    "affected": re.compile(
        r"\b(?:affected|symptomatic|phenotype[- ]positive)\b",
        re.IGNORECASE,
    ),
    "unaffected": re.compile(
        r"\b(?:unaffected|asymptomatic|phenotype[- ]negative)\b",
        re.IGNORECASE,
    ),
}


class CountRecoveryError(RuntimeError):
    """The model's response cannot support an honest count fill."""


#: Layers whose rows were actually read out of the paper. Rows that exist only
#: via ClinVar / PubTator linkage name a variant the paper may never mention, so
#: asking a model to find its count *in that paper* is meaningless -- measured on
#: KCNH2 10973849, 583 of its "gaps" were linkage rows whose notation does not
#: occur anywhere in the 51,783-char source, and the model correctly returned
#: nothing for every one of them.
PAPER_DERIVED_LAYERS = frozenset(
    {"llm_table", "llm_text", "regex_table", "regex_text", "figure", "mixed"}
)


@dataclass
class VariantGap:
    """One stored variant that is missing at least one count for one paper.

    ``variant_id`` is the row a recovered count is written to. ``duplicate_ids``
    are other ``variants`` rows in the same paper that denote the same variant
    under a different notation; they are deliberately NOT written to, because a
    second count row for one variant in one paper inflates carrier aggregates.
    """

    variant_id: int
    notation: str
    missing: list[str]  # subset of RECOVERABLE_FIELDS keys
    duplicate_ids: list[int] = field(default_factory=list)
    parts: tuple[str, ...] = ()  # component notations, for identity matching
    paper_derived: bool = False


@dataclass
class PaperGap:
    gene: str
    pmid: str
    variants: list[VariantGap]

    @property
    def missing_total(self) -> int:
        return sum(len(v.missing) for v in self.variants)


@dataclass
class RecoveredCount:
    """One accepted fill, with the structured role and locator that justify it.

    ``count_role`` and ``evidence_locator`` are not decoration: the trust gate
    reads role from ``variant_papers.count_provenance``, so a recovered value
    written without them cannot be scored for role consistency and stays
    ``trusted`` by DDL default no matter how it was derived.
    """

    variant_id: int
    notation: str
    field: str
    value: int
    quote: str
    count_role: str = ""
    evidence_locator: dict = field(default_factory=dict)
    model_declared_role: Optional[str] = None


@dataclass
class PaperResult:
    gene: str
    pmid: str
    accepted: list[RecoveredCount] = field(default_factory=list)
    rejected: list[dict] = field(default_factory=list)
    error: Optional[str] = None


def _normalize(text: str) -> str:
    return _WS_RE.sub(" ", text or "").strip().lower()


def find_count_gaps(
    db: Path | str,
    gene: str,
    *,
    pmids: Optional[Iterable[str]] = None,
    fields: Sequence[str] = ("carriers",),
    paper_derived_only: bool = True,
) -> list[PaperGap]:
    """Return per-paper variants missing one or more of ``fields``.

    A variant counts as missing a field when it has no ``penetrance_data`` row
    for that paper at all, or has one whose column is NULL.

    ``paper_derived_only`` (default) restricts the result to variants at least
    one of whose rows came from reading the paper -- see
    :data:`PAPER_DERIVED_LAYERS`. Pass False only to inspect the full gap set;
    recovering counts for linkage-only variants asks the model to find numbers
    the paper does not contain.
    """
    unknown = set(fields) - set(RECOVERABLE_FIELDS)
    if unknown:
        raise ValueError(f"unrecoverable fields: {sorted(unknown)}")

    con = sqlite3.connect(f"file:{Path(db)}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            """
            SELECT vp.pmid AS pmid, vp.variant_id AS variant_id,
                   vp.source_layer AS source_layer,
                   v.protein_notation AS protein, v.cdna_notation AS cdna,
                   pd.total_carriers_observed AS carriers,
                   pd.affected_count AS affected,
                   pd.unaffected_count AS unaffected
            FROM variant_papers vp
            JOIN variants v ON v.variant_id = vp.variant_id
            LEFT JOIN penetrance_data pd
                   ON pd.variant_id = vp.variant_id AND pd.pmid = vp.pmid
            WHERE UPPER(COALESCE(v.gene_symbol, '')) = UPPER(?)
            """,
            (gene,),
        ).fetchall()
    finally:
        con.close()

    wanted = {str(p) for p in pmids} if pmids is not None else None
    by_pmid: dict[str, dict[int, VariantGap]] = {}
    for r in rows:
        pmid = str(r["pmid"])
        if wanted is not None and pmid not in wanted:
            continue
        parts = tuple(
            dict.fromkeys(
                p.strip() for p in (r["protein"], r["cdna"]) if p and p.strip()
            )
        )
        notation = " ".join(parts)
        if not notation:
            continue  # nothing to name the variant by; cannot ask about it
        layers = {t.strip() for t in (r["source_layer"] or "").split(",") if t.strip()}
        is_paper_derived = bool(layers & PAPER_DERIVED_LAYERS)
        missing = [f for f in fields if r[f] is None]
        # A variant can appear on several joined rows (multiple source layers,
        # multiple penetrance rows). A field is only a gap when it is NULL on
        # *every* row, so a row that carries the count must still be folded in
        # -- skipping it here would report a populated variant as missing.
        slot = by_pmid.setdefault(pmid, {})
        existing = slot.get(r["variant_id"])
        if existing is None:
            slot[r["variant_id"]] = VariantGap(
                r["variant_id"],
                notation,
                list(missing),
                parts=parts,
                paper_derived=is_paper_derived,
            )
        else:
            existing.missing = [f for f in existing.missing if f in missing]
            existing.paper_derived = existing.paper_derived or is_paper_derived
            if len(notation) > len(existing.notation):
                existing.notation, existing.parts = notation, parts

    out = []
    for pmid, slot in sorted(by_pmid.items()):
        out.append(
            PaperGap(
                gene,
                pmid,
                _collapse_duplicate_variants(
                    list(slot.values()),
                    gene,
                    fields,
                    paper_derived_only=paper_derived_only,
                ),
            )
        )
    return [p for p in out if p.variants]


def _gap_rank(gap: VariantGap) -> tuple[int, int]:
    """Prefer the richest notation as a cluster's representative row."""
    return (-len(gap.notation), gap.variant_id)


def _same_variant(a: VariantGap, b: VariantGap, gene: str) -> bool:
    """Compare on component notations, never on the combined display string.

    ``notation`` is a space-joined "p.Leu552Ser c.1655T>C" built for the prompt;
    ``variants_match`` cannot parse that (it returns False against "p.L552S"
    even though the protein halves clearly agree). Matching each component
    pairwise is what actually detects the duplicate.
    """
    from utils.variant_normalizer import variants_match

    for left in a.parts or (a.notation,):
        for right in b.parts or (b.notation,):
            try:
                if variants_match(left, right, gene):
                    return True
            except Exception:  # noqa: BLE001 - a normalizer failure must not merge
                continue
            # Protein-only and cDNA-only rows cannot be proven identical without
            # a transcript-aware consequence mapper. Conservatively collapse a
            # same-codon pair anyway: suppressing a possible recovery is safer
            # than writing two carrier counts for one biological variant.
            if _cdna_protein_positions_match(left, right):
                return True
    return False


def _collapse_duplicate_variants(
    entries: list[VariantGap],
    gene: str,
    fields: Sequence[str],
    *,
    paper_derived_only: bool,
) -> list[VariantGap]:
    """Merge entries whose notations denote the same variant.

    Production stores one variant under several ``variant_id``s when different
    layers write different notations — e.g. ``p.Leu552Ser c.1655T>C`` from
    llm_text alongside ``p.L552S`` from a PubTator linkage. Treating those as
    two gaps is worse than wasteful: if the extractor already counted one of
    them, filling the other creates a second carrier count for the same variant
    in the same paper, which inflates penetrance. So a field is only a gap when
    it is missing across the whole cluster, and only the representative row is
    ever written to.
    """
    clusters: list[list[VariantGap]] = []
    for entry in sorted(entries, key=_gap_rank):
        for cluster in clusters:
            if _same_variant(entry, cluster[0], gene):
                cluster.append(entry)
                break
        else:
            clusters.append([entry])

    collapsed: list[VariantGap] = []
    for cluster in clusters:
        eligible = (
            [m for m in cluster if m.paper_derived] if paper_derived_only else cluster
        )
        if not eligible:
            continue
        rep = min(eligible, key=_gap_rank)
        # Missing only where every member of the cluster is missing it.
        missing = [f for f in fields if all(f in m.missing for m in cluster)]
        if not missing:
            continue
        collapsed.append(
            VariantGap(
                rep.variant_id,
                rep.notation,
                missing,
                duplicate_ids=[m.variant_id for m in cluster if m is not rep],
                parts=rep.parts,
                paper_derived=rep.paper_derived,
            )
        )
    return collapsed


def build_prompt(
    gene: str,
    pmid: str,
    source_text: str,
    variants: Sequence[VariantGap],
    *,
    fields: Sequence[str] = ("carriers",),
    max_source_chars: int = 120_000,
) -> str:
    listing = "\n".join(
        f'  {i}. "{v.notation}"  (need: {", ".join(v.missing)})'
        for i, v in enumerate(variants, 1)
    )
    field_defs = {
        "carriers": (
            "carriers = number of HUMAN individuals reported to carry this "
            "variant in this paper"
        ),
        "affected": (
            "affected = of those carriers, how many are reported symptomatic / "
            "phenotype-positive"
        ),
        "unaffected": (
            "unaffected = of those carriers, how many are reported "
            "asymptomatic / phenotype-negative"
        ),
    }
    defs = "\n".join(f"  - {field_defs[f]}" for f in fields if f in field_defs)
    source = source_text[:max_source_chars]
    truncated = len(source_text) > max_source_chars

    return f"""You are completing missing counts in an existing {gene} curation record.

The variants below were ALREADY extracted from this paper (PMID {pmid}). Do not
look for new variants and do not correct the notations. For each listed variant,
report only the counts that are missing, or null.

Variants:
{listing}

Definitions:
{defs}

Hard rules:
1. Every non-null count MUST be supported by a verbatim quote copied exactly
   from the source below. The number must be stated in that quote, either as
   digits or as a spelled-out number ("once", "two patients").
2. The quote must be one contiguous sentence or the target variant's local table
   row/cell group. It must name that variant (the paper's equivalent notation is
   fine). Never quote a whole table or a detached block of variant names and
   numbers; if row alignment is lost, return null.
3. If the paper does not state the count for that specific variant, return null.
   Never infer, apportion, or compute a count. Never use a study-wide N, a
   cohort total, a screened denominator, a family count, an allele count, or a
   population/gnomAD frequency as a per-variant count.
4. Counts are per-variant and per-this-paper only.
5. Report only variants from the list above, using the notation exactly as given.
6. Declare the role of the number you report in "count_role". Use
   "per_variant_carriers", "per_variant_affected" or "per_variant_unaffected"
   ONLY when the paper states that number for this one variant. If the number is
   a cohort total, a denominator ("12 of 812"), a family count, an allele
   count/number/frequency, or a measurement (age, QTc, follow-up), say so with
   "cohort_total", "denominator", "family_count", "allele_count" or
   "measurement" and set the count to null.
7. In "evidence_locator", say where the number is: for a table give
   {{"kind": "table_row", "table": "<caption or id>", "column": "<header>",
   "row": "<the target variant's row>"}}; for prose give
   {{"kind": "sentence", "section": "<section name if known>"}}.

Return strict JSON, no prose and no markdown fences:

{{"variants": [{{"variant": "<notation exactly as listed>",
                "carriers": <int|null>, "affected": <int|null>,
                "unaffected": <int|null>,
                "count_role": "<per_variant_carriers|per_variant_affected|"
                              "per_variant_unaffected|cohort_total|denominator|"
                              "family_count|allele_count|measurement>",
                "evidence_locator": {{...}},
                "quote": "<verbatim span supporting the numbers, or null>"}}]}}

If no count can be grounded for any listed variant, return {{"variants": []}}.

--- SOURCE{" (truncated)" if truncated else ""} ---
{source}
"""


def parse_response(raw: str) -> list[dict]:
    """Extract the variants array, tolerating fences and prose preface."""
    if not raw or not raw.strip():
        raise CountRecoveryError("empty response")
    text = raw.strip()
    if text.startswith("```"):
        text = re.sub(r"^```[a-zA-Z0-9_]*\s*", "", text)
        text = re.sub(r"\s*```\s*$", "", text)
    decoder = json.JSONDecoder()
    for start in (i for i, ch in enumerate(text) if ch == "{"):
        try:
            obj, _ = decoder.raw_decode(text[start:])
        except ValueError:
            continue
        if isinstance(obj, dict) and "variants" in obj:
            got = obj["variants"]
            if not isinstance(got, list):
                raise CountRecoveryError("'variants' is not a list")
            return [v for v in got if isinstance(v, dict)]
    raise CountRecoveryError("no JSON object with a 'variants' key")


def _single_table_row(quote: str) -> bool:
    lines = [line.strip() for line in quote.splitlines() if line.strip()]
    return len(lines) == 1 and "|" in lines[0]


def _adjacent_word(text: str, offset: int, *, before: bool) -> Optional[str]:
    fragment = text[:offset] if before else text[offset:]
    pattern = r"([a-z]+)\s*$" if before else r"^\s*([a-z]+)"
    match = re.search(pattern, fragment)
    return match.group(1) if match else None


def _number_occurrences(quote: str, value: int) -> list[tuple[int, int]]:
    """Return spans that state ``value`` without matching HGVS positions."""
    norm = _normalize(quote)
    spans: list[tuple[int, int]] = []
    for match in _DIGIT_COUNT_RE.finditer(norm):
        try:
            if int(match.group(0).replace(",", "")) == value:
                spans.append(match.span())
        except ValueError:
            continue

    for phrase, number in sorted(_WORD_NUMBERS.items(), key=lambda item: -len(item[0])):
        if number != value:
            continue
        pattern = re.compile(rf"(?<!\w){re.escape(phrase)}(?!\w)")
        for match in pattern.finditer(norm):
            before = _adjacent_word(norm, match.start(), before=True)
            after = _adjacent_word(norm, match.end(), before=False)
            if before in _COMPOUND_NUMBER_TOKENS or after in _COMPOUND_NUMBER_TOKENS:
                continue
            spans.append(match.span())
    return spans


def _non_count_label(text: str) -> bool:
    """True when a header/label/segment names something other than a per-variant count.

    Used on the *table* branch, where the label is the whole role signal. This is
    the check whose absence let ``| p.Leu552Ser | allele count | 12 |``,
    ``| Variant | gnomAD AC |``, ``| Variant | Age at diagnosis |`` and
    ``| p.Leu552Ser | families screened | 3 |`` all be written as carrier counts.
    """
    if not text:
        return False
    return bool(
        _NON_COUNT_CONTEXT_RE.search(text)
        or _NON_COUNT_TABLE_LABEL_RE.search(text)
        or _MEASURE_UNIT_RE.search(text)
    )


def _is_denominator(
    norm: str,
    start: int,
    end: int,
    other_spans: Sequence[tuple[int, int]],
) -> bool:
    """True when ``norm[start:end]`` is a cohort denominator, not a count.

    ``12 of 812 subjects``, ``5/120 patients``, ``among 913 individuals``,
    ``44 cases of Long QT`` — every one of these was accepted as a per-variant
    carrier count by the prior gate.
    """
    before = norm[:start]
    if _STRONG_DENOMINATOR_MARKER_RE.search(before):
        return True
    weak = _WEAK_DENOMINATOR_MARKER_RE.search(before)
    if weak:
        prefix = before[: weak.start()]
        if _N_EQUALS_TAIL_RE.search(prefix):
            return True
        # "<int> of <int>" — a preceding candidate integer makes the marker a
        # ratio. Bare "identified in 44 carriers" has no such integer, so it is
        # not treated as a denominator.
        if any(
            len(prefix) - 12 <= span_end <= len(prefix) for _, span_end in other_spans
        ):
            return True
    return bool(_COHORT_HEAD_RE.match(norm[end:]))


def _local_segment(
    norm: str,
    start: int,
    end: int,
    other_spans: Sequence[tuple[int, int]],
) -> str:
    """The text around one number, bounded by the nearest other candidate numbers.

    Bounding at the neighbours is what stops ``n=7 among 913 individuals carried
    p.Leu552Ser`` from lending "carried" to the 7: the segment for 7 ends at 913.
    """
    left = max(
        (span_end for _, span_end in other_spans if span_end <= start), default=0
    )
    right = min(
        (span_start for span_start, _ in other_spans if span_start >= end),
        default=len(norm),
    )
    return norm[left:right]


def _has_explicit_role(
    norm: str,
    start: int,
    end: int,
    other_spans: Sequence[tuple[int, int]],
    field: str,
) -> bool:
    pattern = _EXPLICIT_ROLE_RE.get(field)
    if pattern is None:
        return False
    return bool(pattern.search(_local_segment(norm, start, end, other_spans)))


def quote_supports(quote: str, value: int, *, field: Optional[str] = None) -> bool:
    """True when ``value`` is stated as a count, not a notation position.

    ``field=None`` means "the caller already established the role" (a header, an
    inline label, or a single-number local table segment) and only asks whether
    the number is stated here. It is NOT a bypass: the negative-context checks
    still run, because the prior ``field is None -> return True`` short-circuit
    was how every table-branch value reached acceptance with no role check.

    With a ``field``, three things must hold for at least one occurrence:

    1. no population/denominator vocabulary within ±35 characters, and no
       measurement/annotation word immediately beside it;
    2. the number is not in a denominator grammar (``X of Y``, ``X/Y``,
       ``among Y``, ``N cases of <phenotype>``);
    3. if the quote states more than one candidate integer, the requested one
       carries an explicit role label in its own local segment — otherwise the
       pass fails closed, matching what the table branch already did for rows
       with several unlabeled integers.
    """
    if not quote:
        return False
    norm = _normalize(quote)
    occurrences = _number_occurrences(norm, value)
    if not occurrences:
        return False

    context_pattern = None
    if field is not None:
        context_pattern = _FIELD_CONTEXT_RE.get(field)
        if context_pattern is None:
            return False

    all_spans = sorted(_number_spans(norm))
    ambiguous = len(all_spans) > 1
    for start, end in occurrences:
        others = [
            span for span in all_spans if not (span[0] <= start and end <= span[1])
        ]
        if _NON_COUNT_CONTEXT_RE.search(
            norm[max(0, start - 35) : min(len(norm), end + 35)]
        ):
            continue
        if _MEASURE_UNIT_RE.search(
            norm[
                max(0, start - _MEASURE_UNIT_WINDOW) : min(
                    len(norm), end + _MEASURE_UNIT_WINDOW
                )
            ]
        ):
            continue
        if field is None:
            return True
        if _is_denominator(norm, start, end, others):
            continue
        if ambiguous and not _has_explicit_role(norm, start, end, others, field):
            continue
        window = norm[max(0, start - 90) : min(len(norm), end + 90)]
        if context_pattern.search(window):
            return True
    return False


def _cdna_protein_positions_match(left: str, right: str) -> bool:
    """Flag a possible cDNA/protein same-codon alias for deduplication only."""
    cdna, protein = (left, right) if "c." in left.lower() else (right, left)
    cdna_match = re.search(r"\bc\.\(?(\d+)", cdna, re.IGNORECASE)
    protein_match = re.search(
        r"(?:\bp\.\(?)?[A-Za-z]{1,3}(\d+)[A-Za-z*]{1,3}",
        protein,
    )
    if not cdna_match or not protein_match:
        return False
    return (int(cdna_match.group(1)) + 2) // 3 == int(protein_match.group(1))


def _notation_matches_target(notation: str, target: VariantGap, gene: str) -> bool:
    from utils.variant_normalizer import variants_match

    parts = target.parts or (target.notation,)
    for part in parts:
        try:
            if variants_match(notation, part, gene):
                return True
        except Exception:  # noqa: BLE001 - malformed scanner tokens fail closed
            continue
    return False


def _literal_variant_forms(notation: str, gene: str) -> set[str]:
    from utils.variant_normalizer import normalize_variant

    forms = {_normalize(notation)}
    try:
        forms.add(_normalize(normalize_variant(notation, gene)))
    except Exception:  # noqa: BLE001 - raw form remains available
        pass
    expanded = set(forms)
    for form in forms:
        if form.startswith(("p.", "c.")):
            expanded.add(form[2:])
        expanded.add(form.replace("(", "").replace(")", ""))
        expanded.add(form.replace("\\", "").replace("(", "").replace(")", ""))
        expanded.add(
            form.replace("\\", "").replace("(", "").replace(")", "").replace("_", "")
        )
    return {form for form in expanded if form}


def _text_mentions_target(text: str, target: VariantGap, gene: str) -> bool:
    if _notation_matches_target(text, target, gene):
        return True
    norm_text = _normalize(text)
    text_forms = {
        norm_text,
        norm_text.replace("\\", "").replace("(", "").replace(")", ""),
        norm_text.replace("\\", "").replace("(", "").replace(")", "").replace("_", ""),
    }
    for part in target.parts or (target.notation,):
        for form in _literal_variant_forms(part, gene):
            for text_form in text_forms:
                if re.search(
                    rf"(?<![a-z0-9]){re.escape(form)}(?![a-z0-9])",
                    text_form,
                ):
                    return True
                if re.search(r"(?:ins|del)$", form) and form in text_form:
                    return True
    for candidate in re.findall(
        r"(?<![A-Za-z0-9.])(?:c\.)?\d+[ACGT]>[ACGT](?![A-Za-z0-9])",
        text,
        re.IGNORECASE,
    ):
        if _notation_matches_target(
            candidate if candidate.lower().startswith("c.") else f"c.{candidate}",
            target,
            gene,
        ):
            return True
    return False


def _variant_position(notation: str) -> Optional[int]:
    cdna_match = re.search(r"\bc\.\(?(\d+)", notation, re.IGNORECASE)
    if cdna_match:
        return (int(cdna_match.group(1)) + 2) // 3
    protein_match = re.search(
        r"(?:\bp\.\(?)?[A-Za-z]{1,3}[\s._-]*(\d+)",
        notation,
        re.IGNORECASE,
    )
    return int(protein_match.group(1)) if protein_match else None


def _explicitly_binds_target_to_other_gene(
    quote: str,
    target: VariantGap,
    gene: str,
) -> bool:
    """Reject a matching notation explicitly prefixed by another gene symbol."""
    expected = gene.upper()
    return any(
        match.group("gene").upper() != expected
        and _notation_matches_target(match.group("variant"), target, gene)
        for match in _GENE_PREFIXED_VARIANT_RE.finditer(quote)
    )


def _table_cells(line: str) -> list[str]:
    cells = [cell.strip() for cell in line.split("|")]
    if cells and not cells[0]:
        cells = cells[1:]
    if cells and not cells[-1]:
        cells = cells[:-1]
    return cells


def _number_spans(norm: str) -> set[tuple[int, int]]:
    """Every standalone candidate-number span in already-normalized text.

    Excludes HGVS position digits and decimal fragments (``0.21``), and drops
    word numbers that are part of a compound (``one hundred``).
    """
    spans = {match.span() for match in _DIGIT_COUNT_RE.finditer(norm)}
    for phrase in sorted(_WORD_NUMBERS, key=len, reverse=True):
        pattern = re.compile(rf"(?<!\w){re.escape(phrase)}(?!\w)")
        for match in pattern.finditer(norm):
            before = _adjacent_word(norm, match.start(), before=True)
            after = _adjacent_word(norm, match.end(), before=False)
            if before in _COMPOUND_NUMBER_TOKENS or after in _COMPOUND_NUMBER_TOKENS:
                continue
            if any(
                start <= match.start() and match.end() <= end for start, end in spans
            ):
                continue
            spans.add(match.span())
    return spans


def _table_number_spans(text: str) -> set[tuple[int, int]]:
    """Return standalone table-number spans, excluding HGVS and decimals."""
    return _number_spans(_normalize(text))


def _table_field_label(cell: str, field: str) -> bool:
    """True when a cell explicitly names the requested per-variant role.

    A cell that also names a denominator, an allele count, or a measurement
    ("Carriers screened", "Age at diagnosis", "gnomAD AC") is refused even when
    the role word is present — the label has to be unambiguous to carry the role
    on its own.
    """
    pattern = _TABLE_FIELD_LABEL_RE.get(field)
    if not (pattern and pattern.search(cell)):
        return False
    return not _non_count_label(cell)


def _labeled_table_value(
    cells: Sequence[str], field: str, value: int
) -> Optional[dict[str, Any]]:
    """Accept a value in, or immediately after, an explicit field label."""
    for index, cell in enumerate(cells):
        if not _table_field_label(cell, field):
            continue
        if quote_supports(cell, value):
            return {"kind": "table_inline_label", "label": cell, "cell": cell}
        if index + 1 < len(cells) and quote_supports(cells[index + 1], value):
            return {
                "kind": "table_inline_label",
                "label": cell,
                "cell": cells[index + 1],
            }
    return None


def _value_column(cells: Sequence[str], value: int) -> Optional[int]:
    for index, cell in enumerate(cells):
        if _number_occurrences(_normalize(cell), value):
            return index
    return None


def _table_row_evidence(
    quote: str,
    target: VariantGap,
    gene: str,
    field: str,
    value: int,
) -> Optional[dict[str, Any]]:
    """Bind a table count only when it is in the target's local cell group.

    Returns a structured evidence locator, or None. Every acceptance path now
    carries a role check: a column header, an inline label, or — for the
    single-number segment — the absence of any non-count vocabulary in the
    segment *and* in the number's own column header.
    """
    lines = [line.strip() for line in quote.splitlines() if "|" in line]
    target_line_indexes: list[int] = []
    for line_index, line in enumerate(lines):
        if _text_mentions_target(line, target, gene):
            target_line_indexes.append(line_index)
        cells = _table_cells(line)
        headers: list[str] = []
        if line_index > 0:
            candidate_headers = _table_cells(lines[line_index - 1])
            if len(candidate_headers) == len(cells):
                headers = candidate_headers
        target_cells = [
            index
            for index, cell in enumerate(cells)
            if _text_mentions_target(cell, target, gene)
        ]
        if not target_cells:
            continue

        # Conventional header row: map the requested role by column name.
        for column, header in enumerate(headers):
            if _table_field_label(header, field) and quote_supports(
                cells[column], value
            ):
                return {
                    "kind": "table_header_column",
                    "header": header,
                    "column": column,
                    "row": line,
                }

        # Inline "carriers | 73" or "carriers: 73" labels are also explicit.
        inline = _labeled_table_value(cells, field, value)
        if inline is not None:
            return {**inline, "row": line}

        # Repeating horizontal records are safe when the target's segment,
        # ending at the next different variant cell, contains exactly one
        # standalone number. Rows with several unlabeled integers are ambiguous
        # (carriers/affected/families/etc.) and fail closed.
        start = min(target_cells)
        end = len(cells)
        for index in range(max(target_cells) + 1, len(cells)):
            cell = cells[index]
            if _VARIANT_LIKE_CELL_RE.search(cell) and not _text_mentions_target(
                cell, target, gene
            ):
                end = index
                break
        segment_cells = list(cells[start:end])
        local_text = " | ".join(segment_cells)
        if (
            len(_table_number_spans(local_text)) == 1
            and not _non_count_label(local_text)
            and quote_supports(local_text, value)
        ):
            column = _value_column(cells, value)
            column_header = headers[column] if headers and column is not None else None
            if column_header is None or not _non_count_label(column_header):
                return {
                    "kind": "table_local_segment",
                    "segment": local_text,
                    "column": column,
                    "header": column_header,
                    "row": line,
                }

    # Some markdown renderers split one logical record into a variant heading
    # followed by label/value rows. Permit at most the next two rows and require
    # an explicit field label in the same row as its value. This retains local
    # groups such as "Genotype-positive subjects | 9" without letting a detached
    # "Study total" or "families | 3" supply the carrier count.
    for line_index in target_line_indexes:
        local_lines = lines[line_index : line_index + 3]
        if not quote_binds_variant("\n".join(local_lines), target, gene):
            continue
        for local_line in local_lines[1:]:
            inline = _labeled_table_value(_table_cells(local_line), field, value)
            if inline is not None:
                return {
                    **inline,
                    "kind": "table_multiline_label",
                    "row": local_line,
                    "heading_row": lines[line_index],
                }
    return None


def _table_row_binds_count(
    quote: str,
    target: VariantGap,
    gene: str,
    field: str,
    value: int,
) -> bool:
    """Boolean view of :func:`_table_row_evidence` (kept for callers/tests)."""
    return _table_row_evidence(quote, target, gene, field, value) is not None


def quote_binds_variant(
    quote: str,
    target: VariantGap,
    gene: str,
    *,
    value: Optional[int] = None,
) -> bool:
    """Require quote variant tokens to denote the target or a safe collective."""
    from utils.variant_scanner import scan_document_for_variants

    if _explicitly_binds_target_to_other_gene(quote, target, gene):
        return False
    exact_target = _text_mentions_target(quote, target, gene)
    try:
        scanned = scan_document_for_variants(quote, gene, source="count_recovery_quote")
    except Exception:  # noqa: BLE001 - exact notation is the only safe fallback
        return exact_target

    if not scanned.variants:
        return exact_target
    matches = [
        any(
            _notation_matches_target(token, target, gene)
            for token in (variant.raw_text, variant.normalized)
            if token
        )
        for variant in scanned.variants
    ]
    if any(matches) and all(matches):
        return True
    if exact_target:
        target_positions = {
            position
            for part in target.parts or (target.notation,)
            if (position := _variant_position(part)) is not None
        }
        scanned_positions = [
            _variant_position(variant.raw_text) or _variant_position(variant.normalized)
            for variant in scanned.variants
        ]
        if (
            target_positions
            and scanned_positions
            and all(position in target_positions for position in scanned_positions)
        ):
            return True
    if (
        value == 0
        and any(matches)
        and re.search(
            r"\b(?:none|no)\b.*\b(?:these|the)\s+variants?\b|"
            r"\b(?:these|the)\s+variants?\b.*\b(?:none|no)\b",
            _normalize(quote),
        )
    ):
        return True
    return False


def _quote_has_table_row(quote: str) -> bool:
    return any("|" in line for line in quote.splitlines())


def _count_evidence(
    quote: str,
    target: VariantGap,
    gene: str,
    field: str,
    value: int,
) -> Optional[dict[str, Any]]:
    """Resolve the structured role + locator that justify writing ``value``.

    Returns ``None`` when nothing in the quote proves the number is a
    per-variant count of the requested role. The returned dict is what the
    writer persists into ``count_provenance`` / ``trust_sources`` so the trust
    gate can reason about role, not just magnitude.
    """
    if _quote_has_table_row(quote):
        locator = _table_row_evidence(quote, target, gene, field, value)
    elif quote_binds_variant(quote, target, gene, value=value) and quote_supports(
        quote, value, field=field
    ):
        locator = {"kind": "prose_sentence", "sentence": quote.strip()[:400]}
    else:
        locator = None
    if locator is None:
        return None
    return {"count_role": PER_VARIANT_ROLES[field], "evidence_locator": locator}


def _binding_supports_count(
    quote: str,
    target: VariantGap,
    gene: str,
    field: str,
    value: int,
) -> bool:
    return _count_evidence(quote, target, gene, field, value) is not None


def _binding_rejection_reason(
    quote: str,
    target: VariantGap,
    gene: str,
    field: str,
    value: int,
) -> str:
    """Distinguish attribution failures from absent/invalid count mentions."""
    if _quote_has_table_row(quote):
        if not any(
            _text_mentions_target(line, target, gene)
            for line in quote.splitlines()
            if "|" in line
        ):
            return "quote does not uniquely bind the count to this variant"
        return "number not stated in the target variant's table row"
    if not quote_binds_variant(quote, target, gene, value=value):
        return "quote does not uniquely bind the count to this variant"
    if not quote_supports(quote, value, field=field):
        return "number not stated in the quote as the requested count role"
    return "quote does not provide attributable count evidence"


def validate_paper_response(
    gap: PaperGap,
    payload: Sequence[Mapping[str, Any]],
    source_text: str,
) -> PaperResult:
    """Keep only quote-grounded counts for variants we actually asked about."""
    result = PaperResult(gap.gene, gap.pmid)
    by_notation = {v.notation: v for v in gap.variants}
    norm_source = _normalize(source_text)
    seen: set[tuple[int, str]] = set()

    for entry in payload:
        notation = str(entry.get("variant") or "").strip()
        target = by_notation.get(notation)
        if target is None:  # tolerate whitespace/case drift before giving up
            target = next(
                (
                    v
                    for k, v in by_notation.items()
                    if _normalize(k) == _normalize(notation)
                ),
                None,
            )
        if target is None:
            result.rejected.append(
                {"variant": notation, "reason": "not a variant we asked about"}
            )
            continue

        quote = str(entry.get("quote") or "").strip()
        quote_ok = bool(quote) and _normalize(quote) in norm_source
        declared_role = str(entry.get("count_role") or "").strip().lower() or None
        model_locator = entry.get("evidence_locator")

        for fname in target.missing:
            raw = entry.get(fname)
            if raw is None:
                continue
            # The model's own admission that the number is not a per-variant
            # count is decisive, even when the deterministic gate would pass.
            if declared_role in NON_PER_VARIANT_ROLES:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "count_role": declared_role,
                        "reason": f"model declared a non-per-variant role: {declared_role}",
                    }
                )
                continue
            if (
                declared_role is not None
                and declared_role not in PER_VARIANT_ROLES.values()
            ):
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "count_role": declared_role,
                        "reason": f"unrecognized count_role: {declared_role}",
                    }
                )
                continue
            if declared_role is not None and declared_role != PER_VARIANT_ROLES[fname]:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "count_role": declared_role,
                        "reason": (
                            f"count_role {declared_role} does not match the "
                            f"requested field {fname}"
                        ),
                    }
                )
                continue
            if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": "not a nonnegative int",
                    }
                )
                continue
            if raw > PLAUSIBLE_COUNT_CEILING:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": f"exceeds plausible ceiling {PLAUSIBLE_COUNT_CEILING}",
                    }
                )
                continue
            if not quote_ok:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": "quote missing or not verbatim in source",
                    }
                )
                continue
            if len(quote) > GROUNDING_QUOTE_MAX_CHARS:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": "quote is too long to provide local attribution",
                    }
                )
                continue
            evidence = _count_evidence(quote, target, gap.gene, fname, raw)
            if evidence is None:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": _binding_rejection_reason(
                            quote, target, gap.gene, fname, raw
                        ),
                    }
                )
                continue
            key = (target.variant_id, fname)
            if key in seen:
                result.rejected.append(
                    {
                        "variant": notation,
                        "field": fname,
                        "value": raw,
                        "reason": "duplicate answer for the same variant/field",
                    }
                )
                continue
            seen.add(key)
            locator = dict(evidence["evidence_locator"])
            if isinstance(model_locator, Mapping):
                locator["model_declared"] = json_safe_locator(model_locator)
            result.accepted.append(
                RecoveredCount(
                    target.variant_id,
                    notation,
                    fname,
                    raw,
                    quote,
                    count_role=evidence["count_role"],
                    evidence_locator=locator,
                    model_declared_role=declared_role,
                )
            )
    return result


def json_safe_locator(value: Mapping[str, Any]) -> dict[str, Any]:
    """Flatten a model-declared locator to JSON scalars, bounded in size."""
    out: dict[str, Any] = {}
    for key, item in list(value.items())[:12]:
        if isinstance(item, (str, int, float, bool)) or item is None:
            out[str(key)[:60]] = item if not isinstance(item, str) else item[:300]
        else:
            out[str(key)[:60]] = json.dumps(item, default=str)[:300]
    return out


def _ensure_log_table(cur: sqlite3.Cursor) -> None:
    cur.execute(
        """CREATE TABLE IF NOT EXISTS count_recovery_log(
            variant_id INTEGER, pmid TEXT, field TEXT, column_name TEXT,
            value INTEGER, quote TEXT, model TEXT, reasoning_effort TEXT,
            version TEXT, recovered_at TEXT, count_role TEXT,
            evidence_locator TEXT, model_declared_role TEXT
        )"""
    )
    # Additive migration for logs written by v1/v2 before the role columns.
    existing = {row[1] for row in cur.execute("PRAGMA table_info(count_recovery_log)")}
    for column, decl in (
        ("count_role", "TEXT"),
        ("evidence_locator", "TEXT"),
        ("model_declared_role", "TEXT"),
    ):
        if column not in existing:
            cur.execute(f"ALTER TABLE count_recovery_log ADD COLUMN {column} {decl}")
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_count_recovery_log_pmid "
        "ON count_recovery_log(pmid, variant_id)"
    )
    # penetrance_data has no index at all on (variant_id, pmid), so the write
    # path's per-count SELECT is a full scan. Non-unique on purpose: whether one
    # variant-paper may legitimately carry several cohort rows is an open
    # modelling question (see docs/ARCHITECTURE.md), and UNIQUE would prejudge it.
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_penetrance_variant_pmid "
        "ON penetrance_data(variant_id, pmid)"
    )


def _table_columns(cur: sqlite3.Cursor, table: str) -> set[str]:
    return {row[1] for row in cur.execute(f"PRAGMA table_info({table})")}


def backup_database(
    db: Path | str, *, suffix: str = ".before_count_recovery.db"
) -> Optional[Path]:
    """Copy the DB beside itself before mutation, as the recovery layers do.

    Returns the backup path, or None when the source is missing. A pre-existing
    backup is kept: the first copy is the one that predates any recovery write.
    """
    import shutil

    source = Path(db)
    if not source.is_file():
        return None
    destination = source.with_name(source.name + suffix)
    if destination.exists():
        logger.info(
            "count recovery: reusing existing pre-mutation backup %s", destination
        )
        return destination
    shutil.copy2(source, destination)
    logger.info("count recovery: pre-mutation backup -> %s", destination)
    return destination


def _merge_count_provenance(
    cur: sqlite3.Cursor,
    variant_id: int,
    pmid: str,
    recovered: RecoveredCount,
    *,
    model: str,
) -> bool:
    """Write role + evidence into the provenance the trust gate actually reads.

    ``pipeline.trust_gate`` gates on ``<field>_count_type`` and
    ``<field>_column_label`` read out of ``variant_papers.count_provenance``.
    A recovered value that leaves those NULL fires no rule and stays ``trusted``;
    a value that inherits a stale ``cohort_total`` is spuriously stamped
    ``count_is_total``. Merging (rather than replacing) fixes both directions
    while leaving other fields' provenance intact.
    """
    if "count_provenance" not in _table_columns(cur, "variant_papers"):
        return False
    row = cur.execute(
        "SELECT count_provenance FROM variant_papers WHERE variant_id=? AND pmid=?",
        (variant_id, pmid),
    ).fetchone()
    if row is None:
        return False
    provenance: dict[str, Any] = {}
    if row[0]:
        try:
            parsed = json.loads(row[0])
            if isinstance(parsed, dict):
                provenance = parsed
        except (TypeError, ValueError, json.JSONDecodeError):
            provenance = {}
    locator = recovered.evidence_locator or {}
    label = str(locator.get("header") or locator.get("label") or "")
    provenance[f"{recovered.field}_count_type"] = recovered.count_role
    provenance[f"{recovered.field}_column_label"] = label
    provenance[f"{recovered.field}_source"] = "count_recovery"
    provenance[f"{recovered.field}_recovery"] = {
        "version": COUNT_RECOVERY_VERSION,
        "model": model,
        "evidence_kind": locator.get("kind"),
        "model_declared_role": recovered.model_declared_role,
        "quote_sha256": hashlib.sha256(
            (recovered.quote or "").encode("utf-8")
        ).hexdigest(),
        "evidence_locator": locator,
    }
    cur.execute(
        "UPDATE variant_papers SET count_provenance=? WHERE variant_id=? AND pmid=?",
        (json.dumps(provenance, sort_keys=True), variant_id, pmid),
    )
    return True


def _quarantine_recovered_row(
    cur: sqlite3.Cursor,
    penetrance_id: int,
    trust_columns: set[str],
) -> None:
    """Land a recovered count as ``quarantine`` and let Step 3.7 promote it.

    ``penetrance_data.trust_tier`` has DDL default ``'trusted'``, so an inserted
    or filled row is trusted the moment it exists — before any rule has looked
    at it. Landing in quarantine makes the trust gate the thing that promotes,
    not the thing that would have had to demote.
    """
    if "trust_tier" not in trust_columns:
        return
    # Column names are validated against PRAGMA table_info; values are bound.
    assignments = ["trust_tier=?"]
    parameters: list[Any] = ["quarantine"]
    for column, value in (
        ("trust_reasons", json.dumps(["count_recovery_pending_review"])),
        ("trust_sources", json.dumps(["count_recovery"])),
        ("trust_rule_version", COUNT_RECOVERY_VERSION),
    ):
        if column in trust_columns:
            assignments.append(f"{column}=?")
            parameters.append(value)
    parameters.append(penetrance_id)
    cur.execute(
        f"UPDATE penetrance_data SET {', '.join(assignments)} WHERE penetrance_id=?",
        parameters,
    )


def write_recovered_counts(
    db: Path | str,
    results: Sequence[PaperResult],
    *,
    model: str,
    reasoning_effort: Optional[str],
    dry_run: bool = False,
) -> dict:
    """Fill NULL count columns only. Existing values are left untouched.

    The UPDATE carries its own ``IS NULL`` predicate, so a concurrent or
    earlier write always wins and re-running the pass is idempotent.

    ``dry_run`` reports exactly what a real run would do: it opens the database
    read-only and runs the same already-populated check, so ``counts_written``
    is not inflated by counts that would land on an occupied slot. Each paper
    commits under its own ``SAVEPOINT``, so an interrupt after N papers of paid
    model calls keeps those N papers' writes.
    """
    written = skipped = inserted_rows = provenance_written = quarantined = 0
    now = datetime.now(timezone.utc).isoformat()
    path = Path(db)
    con = (
        sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        if dry_run
        else sqlite3.connect(str(path))
    )
    try:
        con.execute("PRAGMA busy_timeout = 30000")
        cur = con.cursor()
        trust_columns: set[str] = set()
        if not dry_run:
            _ensure_log_table(cur)
            con.commit()
        trust_columns = _table_columns(cur, "penetrance_data")
        for res in results:
            if not dry_run and res.accepted:
                cur.execute("SAVEPOINT count_recovery_paper")
            try:
                for rc in res.accepted:
                    column = RECOVERABLE_FIELDS[rc.field]
                    existing_rows = cur.execute(
                        f"SELECT penetrance_id, {column} FROM penetrance_data "
                        "WHERE variant_id=? AND pmid=? ORDER BY penetrance_id",
                        (rc.variant_id, res.pmid),
                    ).fetchall()
                    if existing_rows and any(
                        row[1] is not None for row in existing_rows
                    ):
                        skipped += 1
                        continue
                    if dry_run:
                        written += 1
                        if not existing_rows:
                            inserted_rows += 1
                        continue
                    if not existing_rows:
                        cur.execute(
                            f"INSERT INTO penetrance_data (variant_id, pmid, {column}) "
                            "VALUES (?,?,?)",
                            (rc.variant_id, res.pmid, rc.value),
                        )
                        inserted_rows += 1
                        penetrance_id = cur.lastrowid
                        changed = 1
                    else:
                        penetrance_id = existing_rows[0][0]
                        cur.execute(
                            f"UPDATE penetrance_data SET {column}=? "
                            f"WHERE penetrance_id=? AND {column} IS NULL",
                            (rc.value, penetrance_id),
                        )
                        changed = cur.rowcount
                    if not changed:
                        skipped += 1
                        continue
                    written += changed
                    if _merge_count_provenance(
                        cur, rc.variant_id, res.pmid, rc, model=model
                    ):
                        provenance_written += 1
                    if "trust_tier" in trust_columns:
                        _quarantine_recovered_row(cur, penetrance_id, trust_columns)
                        quarantined += 1
                    cur.execute(
                        "INSERT INTO count_recovery_log (variant_id, pmid, field, "
                        "column_name, value, quote, model, reasoning_effort, version, "
                        "recovered_at, count_role, evidence_locator, "
                        "model_declared_role) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (
                            rc.variant_id,
                            res.pmid,
                            rc.field,
                            column,
                            rc.value,
                            rc.quote,
                            model,
                            reasoning_effort,
                            COUNT_RECOVERY_VERSION,
                            now,
                            rc.count_role,
                            json.dumps(
                                rc.evidence_locator, sort_keys=True, default=str
                            ),
                            rc.model_declared_role,
                        ),
                    )
            except Exception:
                if not dry_run and res.accepted:
                    cur.execute("ROLLBACK TO count_recovery_paper")
                    cur.execute("RELEASE count_recovery_paper")
                raise
            if not dry_run and res.accepted:
                cur.execute("RELEASE count_recovery_paper")
                con.commit()
    finally:
        con.close()
    return {
        "counts_written": written,
        "already_populated_skipped": skipped,
        "penetrance_rows_inserted": inserted_rows,
        "count_provenance_written": provenance_written,
        "landed_in_quarantine": quarantined,
        "dry_run": dry_run,
    }


#: Source renderings recovery will read, richest-first by size.
SOURCE_BASENAMES = ("{pmid}_FULL_CONTEXT.md", "{pmid}_CLEANED.md")


def corpus_root() -> Path:
    """The shared corpus directory, honouring ``GVF_CORPUS_DIR``."""
    override = os.environ.get("GVF_CORPUS_DIR")
    if override:
        return Path(override)
    # utils/ and pipeline/ are siblings of the repo root in a checkout; in an
    # installed layout there is no repo root, so fall back to the CWD's corpus.
    return Path(__file__).resolve().parents[1] / "corpus"


def default_source_roots(db: Path | str) -> list[Path]:
    """Run-dir ``pmc_fulltext`` first, then the run dir, then the shared corpus."""
    roots: list[Path] = []
    run_dir = Path(db).parent
    for candidate in (run_dir / "pmc_fulltext", run_dir):
        if candidate.is_dir():
            roots.append(candidate)
    roots.append(corpus_root())
    return roots


def make_source_resolver(
    gene: str, source_roots: Sequence[Path]
) -> Callable[[str], Optional[str]]:
    """Return ``pmid -> text``, preferring the richest available rendering.

    Extraction takes the first candidate that clears a size floor, which is why
    a paper can be read from a 25 KB ``_FULL_CONTEXT.md`` while a 37 KB
    ``_CLEANED.md`` sits unread beside it (KCNQ1 17470695 is the known case).
    Recovery deliberately picks the largest candidate instead: a count that only
    appears in the richer rendering is exactly what this pass is for.

    Lives here rather than in ``scripts/`` because ``cli/gvf_run.py`` imports it
    on the Step 3.55 hot path and ``scripts`` is not in the installed wheel.
    """

    def resolve(pmid: str) -> Optional[str]:
        candidates: list[Path] = []
        for root in source_roots:
            for pattern in SOURCE_BASENAMES:
                name = pattern.format(pmid=pmid)
                for path in (
                    root / name,
                    root / gene / pmid / name,
                    root / pmid / name,
                ):
                    if path.is_file():
                        candidates.append(path)
        if not candidates:
            return None
        best = max(candidates, key=lambda p: p.stat().st_size)
        logger.debug("PMID %s source: %s (%d B)", pmid, best, best.stat().st_size)
        return best.read_text(errors="replace")

    return resolve


def result_to_dict(result: PaperResult) -> dict[str, Any]:
    """JSON-serializable view of one paper's recovery outcome.

    ``stats["results"]`` used to hand back live dataclasses, so every consumer
    had to hand-serialise them and any ``json.dumps(stats)`` raised.
    """
    return {
        "gene": result.gene,
        "pmid": result.pmid,
        "error": result.error,
        "accepted": [
            {
                "variant_id": c.variant_id,
                "variant": c.notation,
                "field": c.field,
                "value": c.value,
                "quote": c.quote,
                "count_role": c.count_role,
                "evidence_locator": c.evidence_locator,
                "model_declared_role": c.model_declared_role,
            }
            for c in result.accepted
        ],
        "rejected": list(result.rejected),
    }


def recover_counts(
    db: Path | str,
    gene: str,
    *,
    source_for_pmid: Callable[[str], Optional[str]],
    llm_caller: Callable[..., Any],
    model: str,
    reasoning_effort: Optional[str] = "high",
    pmids: Optional[Iterable[str]] = None,
    fields: Sequence[str] = ("carriers",),
    max_variants_per_call: int = 40,
    max_source_chars: int = 120_000,
    dry_run: bool = False,
    paper_derived_only: bool = True,
    backup: bool = True,
) -> dict:
    """Find count gaps, ask the model to fill them, write grounded answers.

    ``llm_caller`` matches the ``litellm.completion`` shape used elsewhere in
    the pipeline; it is a parameter so the unit suite can drive this offline.
    ``source_for_pmid`` returns the paper text, or None when it is unavailable
    (that paper is then reported as skipped rather than silently dropped).
    ``paper_derived_only=False`` inspects the full gap set including
    ClinVar/PubTator linkage rows — a diagnostic mode, not a recovery mode.

    Every model call runs inside a ``count_recovery`` trace scope carrying gene,
    PMID and attempt, and each paper emits a ``count_recovery_decision`` event
    linking the accepted response. This is the judgment-heaviest write stage in
    the pipeline; it must be the most auditable, not the least.
    """
    from utils.llm_trace import (
        OUTCOME_PARSE_FAILED,
        OUTCOME_PARSED,
        attempt_link_summary,
        llm_attempt_ledger,
        llm_trace_scope,
        note_llm_outcome,
        record_trace_event,
    )

    gaps = find_count_gaps(
        db, gene, pmids=pmids, fields=fields, paper_derived_only=paper_derived_only
    )
    results: list[PaperResult] = []
    stats: dict[str, Any] = {
        "papers_with_gaps": len(gaps),
        "gaps_found": sum(g.missing_total for g in gaps),
        "papers_attempted": 0,
        "papers_without_source": 0,
        "papers_failed": 0,
        "batch_failures": 0,
        "llm_calls": 0,
        "backup_path": None,
    }
    failed_pmids: set[str] = set()

    if backup and not dry_run and gaps:
        backup_path = backup_database(db)
        stats["backup_path"] = str(backup_path) if backup_path else None

    for gap in gaps:
        source = source_for_pmid(gap.pmid)
        if not source:
            stats["papers_without_source"] += 1
            logger.info(
                "count recovery: PMID %s - no source available; %d gap(s) left as-is",
                gap.pmid,
                gap.missing_total,
            )
            continue
        stats["papers_attempted"] += 1
        for index, start in enumerate(
            range(0, len(gap.variants), max_variants_per_call), start=1
        ):
            batch = gap.variants[start : start + max_variants_per_call]
            chunk = PaperGap(gap.gene, gap.pmid, batch)
            prompt = build_prompt(
                gene,
                gap.pmid,
                source,
                batch,
                fields=fields,
                max_source_chars=max_source_chars,
            )
            with (
                llm_trace_scope(
                    gene=gene,
                    pmid=gap.pmid,
                    stage="count_recovery",
                    component="count_recovery",
                    operation="fill_missing_counts",
                    batch=index,
                    attempt=index,
                ),
                llm_attempt_ledger(),
            ):
                try:
                    stats["llm_calls"] += 1
                    raw, call_trace_id = _call_text(
                        llm_caller, model, prompt, reasoning_effort
                    )
                    try:
                        payload = parse_response(raw)
                    except Exception:
                        # A provider success whose content will not parse is NOT
                        # an accepted response. Marking it accepted inside
                        # _call_text let the failure decision below still claim
                        # an accepted_response_trace_id.
                        note_llm_outcome(call_trace_id, OUTCOME_PARSE_FAILED)
                        raise
                    note_llm_outcome(call_trace_id, OUTCOME_PARSED)
                except Exception as exc:  # noqa: BLE001 - one bad paper must not stop the pass
                    stats["batch_failures"] += 1
                    failed_pmids.add(gap.pmid)
                    results.append(
                        PaperResult(
                            gene, gap.pmid, error=f"{type(exc).__name__}: {exc}"
                        )
                    )
                    logger.warning(
                        "count recovery: PMID %s - %s: %s",
                        gap.pmid,
                        type(exc).__name__,
                        exc,
                    )
                    record_trace_event(
                        "count_recovery_decision",
                        {
                            "batch": index,
                            "gaps_requested": chunk.missing_total,
                            "counts_accepted": 0,
                            "counts_rejected": 0,
                            "failure": f"{type(exc).__name__}: {exc}",
                            "decision_source": "batch_failed",
                            **attempt_link_summary(),
                        },
                    )
                    continue
                # Only a quote from source text actually shown to the model can
                # ground its answer. Validating against the unseen truncated tail
                # would let a lucky hallucination pass as verbatim evidence.
                res = validate_paper_response(chunk, payload, source[:max_source_chars])
                results.append(res)
                logger.info(
                    "count recovery: PMID %s - %d/%d gap(s) grounded (%d rejected)",
                    gap.pmid,
                    len(res.accepted),
                    chunk.missing_total,
                    len(res.rejected),
                )
                record_trace_event(
                    "count_recovery_decision",
                    {
                        "batch": index,
                        "gaps_requested": chunk.missing_total,
                        "counts_accepted": len(res.accepted),
                        "counts_rejected": len(res.rejected),
                        "accepted": [
                            {
                                "variant": c.notation,
                                "field": c.field,
                                "value": c.value,
                                "count_role": c.count_role,
                                "evidence_kind": (c.evidence_locator or {}).get("kind"),
                                "count_rationale": c.quote,
                            }
                            for c in res.accepted
                        ],
                        "rejection_reasons": sorted(
                            {str(r.get("reason")) for r in res.rejected}
                        ),
                        "decision_source": "quote_grounded_validation",
                        "selection_policy": (
                            "Only counts whose verbatim quote locally binds the "
                            "value to the requested variant AND proves a "
                            "per-variant role are written; everything else is "
                            "rejected, never stored with a caveat."
                        ),
                        **attempt_link_summary(),
                    },
                )

    stats["papers_failed"] = len(failed_pmids)
    write_stats = write_recovered_counts(
        db, results, model=model, reasoning_effort=reasoning_effort, dry_run=dry_run
    )
    stats.update(write_stats)
    stats["counts_accepted"] = sum(len(r.accepted) for r in results)
    stats["counts_rejected"] = sum(len(r.rejected) for r in results)
    stats["model"] = model
    stats["reasoning_effort"] = reasoning_effort
    stats["version"] = COUNT_RECOVERY_VERSION
    stats["paper_derived_only"] = paper_derived_only
    # `results` is JSON-serializable; `result_objects` keeps the dataclasses for
    # in-process callers that want them.
    stats["results"] = [result_to_dict(r) for r in results]
    stats["result_objects"] = results
    stats["all_batches_failed"] = bool(
        stats["papers_attempted"] and not stats["counts_accepted"] and failed_pmids
    )
    return stats


def _call_text(
    llm_caller: Callable[..., Any],
    model: str,
    prompt: str,
    reasoning_effort: Optional[str],
) -> tuple[str, Optional[str]]:
    """Run one recovery call and return ``(raw_text, trace_id)``.

    Deliberately records only the PROVIDER outcome (``error``, or nothing yet on
    success). Whether the response is an accepted decision input depends on
    parsing and on quote/role validation, which happen in the caller — so the
    caller sets ``parsed`` / ``parse_failed``.
    """
    from utils.llm_trace import (
        OUTCOME_ERROR,
        last_llm_trace,
        note_llm_attempt,
        note_llm_outcome,
    )
    from utils.llm_utils import build_reasoning_effort_kwargs

    kwargs: dict[str, Any] = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.0,
        "max_tokens": 8192,
    }
    kwargs.update(build_reasoning_effort_kwargs(model, reasoning_effort))
    try:
        response = llm_caller(**kwargs)
    except BaseException:
        trace_id = note_llm_attempt(last_llm_trace(), role="count_recovery")
        note_llm_outcome(trace_id, OUTCOME_ERROR)
        raise
    trace_id = note_llm_attempt(last_llm_trace(), role="count_recovery")
    return response.choices[0].message.content, trace_id
