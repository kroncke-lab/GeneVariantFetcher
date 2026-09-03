"""Cross-route variant identity: when two rows are one allele, and when not.

The problem
-----------
Several routes observe the same paper, and ``variants`` rows were keyed on the
**raw** notation string. ``normalize_variant`` already proves ``p.Q403*`` and
``p.Gln403Ter`` are the same call, yet both were inserted as separate
identities. In one 50-paper database that produced 134 paper-groups where a
single cDNA mapped to more than one ``variant_id``, and 9 where a single
protein did::

    c.1766A>G  ->  {p.Y589C, p.Tyr589Cys, <no protein>}   3 identities, 1 allele
    c.1207C>T  ->  {p.Q403*, p.Gln403Ter}
    p.A154G    ->  {c.461C>G, <no cDNA>}

A curator meeting the same allele three times splits their attention and
learns to merge by hand, which is exactly the judgement the pipeline must not
delegate.

The rule
--------
Folding is allowed **only for a spelling difference in a fully specified
identity**. It is never allowed as a *completion* of a partial one. The same
database contains::

    c.1156G>A  ->  {p.Glu386Gln, p.Glu386Lys}

which is a genuine source contradiction (``Glu386Gln`` would be ``c.1156G>C``).
A relation keyed on cDNA alone would merge those two, and a third row with no
protein at all is compatible with *both*, so a relation that is transitive
through a missing field would silently pull them together through it. That is
the "NULL bridge", and it is why :func:`fold_decision` is pairwise, requires
both sides to carry the field being compared, and refuses sparse/rich pairs.

Counts are never a reason to fold and never travel across a fold. Two rows that
disagree on any non-null count for the same paper are a contradiction to
report, not an average to take, and a null count is never filled from a
neighbouring row.

Everything that is not a provably safe fold is reported by
:func:`find_identity_conflicts` as a worklist, not resolved.
"""

from __future__ import annotations

import logging
import re
import sqlite3
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

from utils.variant_normalizer import (
    normalize_variant,
    protein_substitution_frameshift_alias,
)

logger = logging.getLogger(__name__)

__all__ = [
    "FoldVerdict",
    "FoldDecision",
    "VariantIdentity",
    "fold_decision",
    "canonical_identity_key",
    "candidate_rows_by_position",
    "find_identity_conflicts",
]


class FoldVerdict(str, Enum):
    """Outcome of comparing two candidate identities."""

    #: Same allele, different spelling of a fully specified identity.
    FOLD = "fold"
    #: Same gene, but a field both sides specify disagrees. Keep both, flag.
    REFUSE_CONFLICT = "refuse_conflict"
    #: One side omits a field the other supplies. Completing it would invent
    #: a transcript event, so keep both.
    REFUSE_SPARSE = "refuse_sparse"
    #: Different genes never fold, whatever the notation says.
    REFUSE_GENE = "refuse_gene"
    #: Nothing comparable (no shared fully specified field).
    REFUSE_UNCOMPARABLE = "refuse_uncomparable"


@dataclass(frozen=True)
class FoldDecision:
    verdict: FoldVerdict
    reason: str

    @property
    def folds(self) -> bool:
        return self.verdict is FoldVerdict.FOLD


@dataclass(frozen=True)
class VariantIdentity:
    """The identity-bearing fields of one ``variants`` row."""

    gene_symbol: Optional[str] = None
    cdna_notation: Optional[str] = None
    protein_notation: Optional[str] = None
    genomic_position: Optional[str] = None
    legacy_notation: Optional[str] = None
    variant_class: Optional[str] = None
    structural_description: Optional[str] = None

    @classmethod
    def from_row(cls, row: Any) -> "VariantIdentity":
        get = row.get if hasattr(row, "get") else (lambda key, _r=row: _r[key])
        return cls(
            gene_symbol=_clean(get("gene_symbol")),
            cdna_notation=_clean(get("cdna_notation")),
            protein_notation=_clean(get("protein_notation")),
            genomic_position=_clean(get("genomic_position")),
            legacy_notation=_clean(get("legacy_notation")),
            variant_class=_clean(get("variant_class")),
            structural_description=_clean(get("structural_description")),
        )


def _clean(value: Any) -> Optional[str]:
    text = str(value).strip() if value is not None else ""
    return text or None


def _canon(value: Optional[str], gene: Optional[str]) -> Optional[str]:
    """Canonical comparison form for one notation field.

    Returns ``None`` only for an absent value. A value the normalizer cannot
    parse falls back to its own stripped text rather than ``None``, because two
    unparseable values must not compare equal to each other: ``None == None``
    would have folded ``c.foo`` onto ``c.bar`` as "matching cDNA".
    """
    if not value:
        return None
    try:
        canonical = normalize_variant(value, gene or "")
    except Exception:  # noqa: BLE001 - a normalizer slip must not fold rows
        canonical = None
    return canonical or str(value).strip()


_POSITION_RE = re.compile(r"\d{1,6}")


def notation_positions(*values: Optional[str]) -> list[str]:
    """Coordinate digits appearing in these notation strings.

    Every spelling of one allele shares its coordinates: ``p.Q403*``,
    ``p.Gln403Ter`` and ``p.Gln403*`` all contain ``403``. That makes the
    position set a cheap, gene-agnostic, spelling-independent way to bound a
    candidate lookup to a handful of rows, which is what lets identity
    resolution compare *canonical* forms instead of raw strings without
    scanning the whole gene for every observation.
    """
    positions: list[str] = []
    for value in values:
        if not value:
            continue
        for match in _POSITION_RE.finditer(str(value)):
            token = match.group(0).lstrip("0") or "0"
            if token not in positions:
                positions.append(token)
    return positions


def canonical_identity_key(
    identity: VariantIdentity,
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Gene plus canonical cDNA/protein/genomic. Used for grouping only.

    Grouping is not folding: rows sharing a key still go through
    :func:`fold_decision` pairwise before anything is treated as one allele.
    """
    gene = (identity.gene_symbol or "").strip().upper() or None
    return (
        gene,
        _canon(identity.cdna_notation, gene),
        _canon(identity.protein_notation, gene),
        (identity.genomic_position or "").strip() or None,
    )


def _protein_equivalent(
    left: Optional[str], right: Optional[str], gene: Optional[str]
) -> bool:
    if not left or not right:
        return False
    if left == right:
        return True
    if _canon(left, gene) is not None and _canon(left, gene) == _canon(right, gene):
        return True
    # A truncated substitution-looking label folded onto the explicit
    # frameshift it abbreviates, proven by the shared ref+position+alt rule.
    try:
        return bool(protein_substitution_frameshift_alias(left, right))
    except Exception:  # noqa: BLE001
        return False


def fold_decision(left: VariantIdentity, right: VariantIdentity) -> FoldDecision:
    """Decide whether two identities are one allele spelled two ways.

    Implements the full conjunction; any single failure refuses. A fold
    requires that **every** field either both sides specify and agree on, or
    neither side specifies. A field present on one side only is
    ``REFUSE_SPARSE``: completing it is inventing evidence, not deduplicating.
    """
    left_gene = (left.gene_symbol or "").strip().upper()
    right_gene = (right.gene_symbol or "").strip().upper()
    if not left_gene or not right_gene or left_gene != right_gene:
        return FoldDecision(
            FoldVerdict.REFUSE_GENE,
            f"different or missing gene ({left_gene or '?'} vs {right_gene or '?'})",
        )
    gene = left_gene

    comparable = False

    # cDNA
    if left.cdna_notation and right.cdna_notation:
        if _canon(left.cdna_notation, gene) != _canon(right.cdna_notation, gene):
            return FoldDecision(
                FoldVerdict.REFUSE_CONFLICT,
                f"cDNA conflict ({left.cdna_notation} vs {right.cdna_notation})",
            )
        comparable = True
    elif bool(left.cdna_notation) != bool(right.cdna_notation):
        return FoldDecision(
            FoldVerdict.REFUSE_SPARSE,
            "one side has no cDNA; completing it would assert a transcript "
            "event the source did not supply",
        )

    # Protein
    if left.protein_notation and right.protein_notation:
        if not _protein_equivalent(left.protein_notation, right.protein_notation, gene):
            return FoldDecision(
                FoldVerdict.REFUSE_CONFLICT,
                "protein conflict "
                f"({left.protein_notation} vs {right.protein_notation})",
            )
        comparable = True
    elif bool(left.protein_notation) != bool(right.protein_notation):
        return FoldDecision(
            FoldVerdict.REFUSE_SPARSE,
            "one side has no protein consequence; a missing field is not a "
            "proven match",
        )

    # Genomic position: exact when both sides carry one, never normalized
    # across assemblies. A one-sided coordinate is NOT a sparse identity: the
    # allele is already pinned by gene plus cDNA plus consequence, so the
    # coordinate is a derived annotation and its absence is an annotation gap
    # rather than an unknown identity. Treating it as sparse kept rows with
    # byte-identical cDNA *and* protein apart, which is the opposite of the
    # duplicate this relation exists to prevent.
    if left.genomic_position and right.genomic_position:
        if left.genomic_position != right.genomic_position:
            return FoldDecision(
                FoldVerdict.REFUSE_CONFLICT,
                "genomic position conflict "
                f"({left.genomic_position} vs {right.genomic_position})",
            )
        comparable = True

    # Legacy (BIC-style) labels are a distinct identity class and never alias
    # onto a coordinate-bearing identity.
    if left.legacy_notation and right.legacy_notation:
        if left.legacy_notation != right.legacy_notation:
            return FoldDecision(
                FoldVerdict.REFUSE_CONFLICT,
                "legacy label conflict "
                f"({left.legacy_notation} vs {right.legacy_notation})",
            )
        comparable = True
    elif bool(left.legacy_notation) != bool(right.legacy_notation):
        return FoldDecision(FoldVerdict.REFUSE_SPARSE, "one side has no legacy label")

    # A declared class disagreement means the two rows describe different
    # events even when a notation string coincides.
    if (
        left.variant_class
        and right.variant_class
        and left.variant_class.lower() != right.variant_class.lower()
    ):
        return FoldDecision(
            FoldVerdict.REFUSE_CONFLICT,
            f"variant class conflict ({left.variant_class} vs {right.variant_class})",
        )

    if not comparable:
        return FoldDecision(
            FoldVerdict.REFUSE_UNCOMPARABLE,
            "no fully specified field is shared by both rows",
        )

    return FoldDecision(
        FoldVerdict.FOLD, "identical fully specified identity, different spelling"
    )


# --- Shared cross-route resolver --------------------------------------------

_CANDIDATE_COLUMNS = (
    "variant_id",
    "gene_symbol",
    "cdna_notation",
    "protein_notation",
    "genomic_position",
    "legacy_notation",
    "variant_class",
    "structural_description",
)


def candidate_rows_by_position(
    con: Any, gene_symbol: str, incoming: VariantIdentity
) -> list[dict[str, Any]]:
    """Stored rows at the same coordinates as ``incoming``, as dicts.

    Shared by the migration writer and the linkage resolver so both offer the
    predicate the same candidate set.
    """
    positions = notation_positions(incoming.cdna_notation, incoming.protein_notation)
    if not positions:
        return []
    clauses = []
    params: list[Any] = [gene_symbol]
    for position in positions[:8]:
        clauses.append("cdna_notation GLOB ? OR protein_notation GLOB ?")
        params.extend([f"*{position}*", f"*{position}*"])
    rows = con.execute(
        f"""
        SELECT {", ".join(_CANDIDATE_COLUMNS)}
        FROM variants
        WHERE gene_symbol = ? AND ({" OR ".join(clauses)})
        ORDER BY variant_id
        """,
        params,
    ).fetchall()
    return [dict(zip(_CANDIDATE_COLUMNS, row)) for row in rows]


def _fold_candidates(con: Any, gene_symbol: str, incoming: VariantIdentity) -> set[int]:
    """Stored ``variant_id``s that are this identity spelled differently.

    Candidates are selected by shared coordinate digits rather than by raw
    string equality. A raw-equality pre-filter can only ever find a duplicate
    whose *other* field happens to be byte-identical, so the commonest real
    duplicate -- the same call written ``p.Q403*`` on one route and
    ``p.Gln403Ter`` on another, with no cDNA on either -- was never offered to
    the predicate and a second identity was minted anyway.

    Positions bound the scan to a handful of rows while being independent of
    spelling, so the comparison itself can use canonical forms.
    """
    return {
        int(record["variant_id"])
        for record in candidate_rows_by_position(con, gene_symbol, incoming)
        if fold_decision(VariantIdentity.from_row(record), incoming).folds
    }


def resolve_variant_identity(
    con: Any,
    gene: str,
    cdna: Optional[str],
    protein: Optional[str],
    *,
    genomic: Optional[str] = None,
) -> int:
    """Return the ``variant_id`` for one observation, creating it if new.

    This is the single identity gate for every route that writes ``variants``.
    Each database-linkage ingest previously carried its own copy of an
    ``ensure_variant`` helper that matched on the **raw** notation string::

        SELECT variant_id FROM variants
        WHERE gene_symbol=? AND cdna_notation IS ? AND protein_notation IS ?

    So a linkage row spelling a call ``p.Gln403Ter`` minted a second identity
    beside the paper row's ``p.Q403*``, and a curator met one allele twice.
    Linkage was the largest contributor to the duplicate groups measured in a
    shipped database.

    Folding obeys the same invariants as the extraction path: a spelling
    difference in a fully specified identity folds, and nothing else does. A
    linkage row carrying only a protein consequence keeps its own identity
    rather than being completed with a neighbouring row's cDNA, and an
    ambiguous or conflicting match is refused and inserted separately.
    """
    gene_symbol = (gene or "").strip().upper()
    # Exact match first: cheap, and the overwhelmingly common case.
    row = con.execute(
        """
        SELECT variant_id FROM variants
        WHERE gene_symbol = ?
          AND cdna_notation IS ?
          AND protein_notation IS ?
          AND genomic_position IS ?
        """,
        (gene_symbol, cdna, protein, genomic),
    ).fetchone()
    if row:
        return int(row[0])

    incoming = VariantIdentity(
        gene_symbol=gene_symbol,
        cdna_notation=_clean(cdna),
        protein_notation=_clean(protein),
        genomic_position=_clean(genomic),
    )
    if incoming.cdna_notation or incoming.protein_notation:
        matches = _fold_candidates(con, gene_symbol, incoming)
        # Ambiguity is a refusal: two stored rows both claiming this spelling
        # means picking a winner arbitrarily.
        if len(matches) == 1:
            return next(iter(matches))

    cursor = con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation, "
        "genomic_position) VALUES (?, ?, ?, ?)",
        (gene_symbol, cdna, protein, genomic),
    )
    return int(cursor.lastrowid)


# --- Detector ---------------------------------------------------------------


def _counts_conflict(rows: list[sqlite3.Row]) -> Optional[str]:
    """Report a per-paper count disagreement across candidate identities."""
    by_pmid: dict[str, list[tuple]] = {}
    for row in rows:
        by_pmid.setdefault(str(row["pmid"]), []).append(
            (
                row["total_carriers_observed"],
                row["affected_count"],
                row["unaffected_count"],
            )
        )
    for pmid, entries in sorted(by_pmid.items()):
        for index, field in enumerate(("carriers", "affected", "unaffected")):
            values = {entry[index] for entry in entries if entry[index] is not None}
            if len(values) > 1:
                return (
                    f"PMID {pmid} reports conflicting {field} "
                    f"({sorted(values)}) across these identities"
                )
    return None


def find_identity_conflicts(db_path: str) -> dict[str, Any]:
    """Report duplicate-identity and contradiction groups. Never mutates.

    Returns a worklist: ``foldable`` groups are pure spelling duplicates that a
    current-code run would no longer create; ``conflicts`` are genuine source
    disagreements a human must adjudicate; ``sparse`` groups are partial
    observations the pipeline deliberately refuses to complete.
    """
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            """
            SELECT v.variant_id, v.gene_symbol, v.cdna_notation, v.protein_notation,
                   v.genomic_position, v.legacy_notation, v.variant_class,
                   v.structural_description
            FROM variants v
            """
        ).fetchall()
        penetrance = con.execute(
            """
            SELECT variant_id, pmid, total_carriers_observed, affected_count,
                   unaffected_count
            FROM penetrance_data
            """
        ).fetchall()
    finally:
        con.close()

    pen_by_variant: dict[int, list[sqlite3.Row]] = {}
    for row in penetrance:
        pen_by_variant.setdefault(int(row["variant_id"]), []).append(row)

    # Group by any shared fully specified field, then decide pairwise.
    buckets: dict[tuple, list[sqlite3.Row]] = {}
    for row in rows:
        gene = (row["gene_symbol"] or "").strip().upper()
        for field in ("cdna_notation", "protein_notation"):
            value = _clean(row[field])
            if not value:
                continue
            buckets.setdefault((gene, field, _canon(value, gene) or value), []).append(
                row
            )

    foldable: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    sparse: list[dict[str, Any]] = []
    seen_pairs: set[tuple[int, int]] = set()

    for (gene, field, key), group in sorted(
        buckets.items(),
        key=lambda item: (str(item[0][0]), str(item[0][1]), str(item[0][2])),
    ):
        if len(group) < 2:
            continue
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                left_row, right_row = group[i], group[j]
                pair = (
                    min(int(left_row["variant_id"]), int(right_row["variant_id"])),
                    max(int(left_row["variant_id"]), int(right_row["variant_id"])),
                )
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                decision = fold_decision(
                    VariantIdentity.from_row(left_row),
                    VariantIdentity.from_row(right_row),
                )
                entry = {
                    "gene_symbol": gene,
                    "grouped_by": field,
                    "key": key,
                    "variant_ids": list(pair),
                    "left": dict(left_row),
                    "right": dict(right_row),
                    "verdict": decision.verdict.value,
                    "reason": decision.reason,
                }
                if decision.verdict is FoldVerdict.FOLD:
                    count_note = _counts_conflict(
                        pen_by_variant.get(pair[0], [])
                        + pen_by_variant.get(pair[1], [])
                    )
                    if count_note:
                        entry["verdict"] = "refuse_count_conflict"
                        entry["reason"] = count_note
                        conflicts.append(entry)
                    else:
                        foldable.append(entry)
                elif decision.verdict is FoldVerdict.REFUSE_SPARSE:
                    sparse.append(entry)
                elif decision.verdict is FoldVerdict.REFUSE_CONFLICT:
                    conflicts.append(entry)

    return {
        "database": str(db_path),
        "variants_examined": len(rows),
        "foldable_spelling_duplicates": len(foldable),
        "identity_conflicts": len(conflicts),
        "sparse_partial_observations": len(sparse),
        "foldable": foldable,
        "conflicts": conflicts,
        "sparse": sparse,
    }
