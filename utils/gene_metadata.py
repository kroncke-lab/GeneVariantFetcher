"""Shared gene metadata for discovery, triage, and specificity QC.

The pipeline should not need one-off alias/protein-length maps in every stage.
This module centralizes conservative built-ins and augments them with a local
VariantFeatures SQLite database when available.
"""

from __future__ import annotations

import logging
import os
import re
import sqlite3
import string
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Iterable, Mapping, Optional

from utils.env_utils import local_data_discovery_disabled

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class VariantFeaturesResidue:
    gene_symbol: str
    position: int
    reference_residues: tuple[str, ...] = ()
    alternate_residues: tuple[str, ...] = ()
    transcripts: tuple[str, ...] = ()
    matched_hgvs_p: bool = False
    matched_hgvs_c: bool = False


@dataclass(frozen=True)
class GeneMetadata:
    symbol: str
    aliases: tuple[str, ...] = ()
    query_aliases: tuple[str, ...] = ()
    protein_length: Optional[int] = None
    canonical_transcript: Optional[str] = None
    refseq_transcripts: tuple[str, ...] = ()
    ensembl_id: Optional[str] = None
    ncbi_id: Optional[str] = None
    protein_ids: tuple[str, ...] = ()
    sources: tuple[str, ...] = field(default_factory=tuple)
    variantfeatures_db: Optional[str] = None


BUILTIN_GENE_METADATA: dict[str, GeneMetadata] = {
    "APOE": GeneMetadata(
        symbol="APOE",
        aliases=("APOE", "ApoE", "APO-E", "apolipoprotein E", "apolipoprotein-E"),
        query_aliases=("APOE", "ApoE", "APO-E", "apolipoprotein E"),
        protein_length=317,
        sources=("builtin",),
    ),
    "BRCA1": GeneMetadata(
        symbol="BRCA1",
        aliases=("BRCA1", "BRCC1", "FANCS", "breast cancer 1"),
        query_aliases=("BRCA1", "BRCAI", "BRCC1", "FANCS", "IRIS", "breast cancer 1"),
        protein_length=1863,
        sources=("builtin",),
    ),
    "MYBPC3": GeneMetadata(
        symbol="MYBPC3",
        aliases=(
            "MYBPC3",
            "MYPBC3",
            "MYBP-C",
            "cMyBP-C",
            "cardiac myosin-binding protein C",
            "cardiac myosin binding protein C",
            "myosin-binding protein C cardiac",
            "myosin binding protein C cardiac",
        ),
        query_aliases=(
            "MYBPC3",
            "MYPBC3",
            "MYBP-C",
            "cMyBP-C",
            "cardiac myosin-binding protein C",
            "cardiac myosin binding protein C",
        ),
        protein_length=1274,
        sources=("builtin",),
    ),
    "KCNH2": GeneMetadata(
        symbol="KCNH2",
        aliases=("KCNH2", "hERG", "HERG", "HERG1", "Kv11.1", "KVLQT2"),
        query_aliases=("KCNH2", "hERG", "HERG", "HERG1", "Kv11.1", "KVLQT2", "LQT2"),
        protein_length=1159,
        sources=("builtin",),
    ),
    "KCNQ1": GeneMetadata(
        symbol="KCNQ1",
        aliases=("KCNQ1", "Kv7.1", "KVLQT1"),
        query_aliases=("KCNQ1", "Kv7.1", "KVLQT1", "LQT1"),
        protein_length=676,
        sources=("builtin",),
    ),
    "SCN5A": GeneMetadata(
        symbol="SCN5A",
        aliases=("SCN5A", "Nav1.5", "NaV1.5"),
        query_aliases=("SCN5A", "Nav1.5", "NaV1.5", "LQT3"),
        protein_length=2016,
        sources=("builtin",),
    ),
    "RYR2": GeneMetadata(
        symbol="RYR2",
        aliases=("RYR2", "RyR2", "ryanodine receptor 2"),
        query_aliases=("RYR2", "RyR2", "ryanodine receptor 2", "CPVT1"),
        protein_length=4967,
        sources=("builtin",),
    ),
    "KCNE1": GeneMetadata(
        symbol="KCNE1",
        aliases=("KCNE1", "minK", "ISK"),
        query_aliases=("KCNE1", "minK", "ISK", "LQT5"),
        protein_length=129,
        sources=("builtin",),
    ),
    "KCNE2": GeneMetadata(
        symbol="KCNE2",
        aliases=("KCNE2", "MiRP1", "MIRP1"),
        query_aliases=("KCNE2", "MiRP1", "MIRP1", "LQT6"),
        protein_length=123,
        sources=("builtin",),
    ),
    "KCNJ2": GeneMetadata(
        symbol="KCNJ2",
        aliases=("KCNJ2", "Kir2.1", "IRK1"),
        query_aliases=("KCNJ2", "Kir2.1", "IRK1", "Andersen-Tawil"),
        protein_length=427,
        sources=("builtin",),
    ),
    "CACNA1C": GeneMetadata(
        symbol="CACNA1C",
        aliases=("CACNA1C", "CaV1.2", "Cav1.2"),
        query_aliases=("CACNA1C", "CaV1.2", "Cav1.2", "LQT8", "Timothy syndrome"),
        protein_length=2221,
        sources=("builtin",),
    ),
    "LDLR": GeneMetadata(
        symbol="LDLR",
        aliases=("LDLR", "LDL receptor", "low density lipoprotein receptor"),
        query_aliases=(
            "LDLR",
            "LDL receptor",
            "low density lipoprotein receptor",
            "FH",
        ),
        protein_length=860,
        sources=("builtin",),
    ),
    "BMPR2": GeneMetadata(
        symbol="BMPR2",
        aliases=(
            "BMPR2",
            "BMPR-2",
            "BMPR-II",
            "BMPRII",
            "BMPR3",
            "BMR2",
            "BRK-3",
            "bone morphogenetic protein receptor type 2",
            "bone morphogenetic protein receptor type II",
        ),
        # PPH1/POVD1 are the historical disease-locus names and T-ALK an older
        # receptor alias; they retrieve real BMPR2 papers but are too
        # disease-shaped to treat as gene mentions during specificity checks.
        query_aliases=(
            "BMPR2",
            "BMPR-II",
            "BMPR3",
            "BMR2",
            "BRK-3",
            "POVD1",
            "PPH1",
            "T-ALK",
            "bone morphogenetic protein receptor type 2",
        ),
        protein_length=1038,
        canonical_transcript="NM_001204.7",
        refseq_transcripts=("NM_001204.7",),
        ncbi_id="659",
        protein_ids=("NP_001195.2",),
        sources=("builtin",),
    ),
}


def normalize_gene_symbol(gene_symbol: str) -> str:
    return (gene_symbol or "").strip().upper()


def default_variantfeatures_db_path() -> Optional[Path]:
    """Return a likely VariantFeatures SQLite path, if configured or present.

    An explicit ``VARIANTFEATURES_DB`` / ``VARIANT_FEATURES_DB`` always wins.
    Absent that, a sibling checkout is only *guessed* at when local-data
    discovery is enabled — see ``local_data_discovery_disabled``.
    """

    for key in ("VARIANTFEATURES_DB", "VARIANT_FEATURES_DB"):
        value = os.environ.get(key)
        if not value:
            continue
        configured = Path(value).expanduser()
        if configured.is_file():
            return configured
        # Explicitly configured but unreadable — typically the external volume
        # holding it is not mounted. Falling through to built-in metadata is the
        # documented behaviour, but doing it silently hides a real misconfiguration.
        parts = configured.parts
        volume = parts[parts.index("Volumes") + 1] if "Volumes" in parts[:-1] else None
        hint = (
            f" That path is on the volume '{volume}' — mount it at "
            f"/Volumes/{volume} and re-run."
            if volume
            else ""
        )
        logger.warning(
            "%s=%s is not a readable file — falling back to built-in gene metadata.%s",
            key,
            configured,
            hint,
        )

    if local_data_discovery_disabled():
        return None

    repo_root = Path(__file__).resolve().parents[1]
    candidates = (
        repo_root / "data" / "variants.db",
        repo_root.parent / "variantFeatures" / "data" / "variants.db",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


@lru_cache(maxsize=128)
def get_gene_metadata(
    gene_symbol: str,
    variantfeatures_db: str | None = None,
) -> GeneMetadata:
    """Return merged metadata for one gene.

    Built-ins are used for safe aliases and known canonical protein lengths.
    VariantFeatures can add transcript IDs, Ensembl/NCBI IDs, and protein length
    for genes not yet in the built-ins. For built-in genes, a suspiciously short
    VariantFeatures slice will not overwrite the curated length.
    """

    gene = normalize_gene_symbol(gene_symbol)
    base = BUILTIN_GENE_METADATA.get(
        gene, GeneMetadata(symbol=gene, aliases=(gene,), query_aliases=(gene,))
    )
    vf = _load_variantfeatures_metadata(gene, variantfeatures_db)
    if not vf:
        return base

    aliases = _dedupe([*base.aliases, *vf.aliases])
    query_aliases = _dedupe([*base.query_aliases, *vf.query_aliases, *aliases])
    protein_length = _choose_protein_length(base.protein_length, vf.protein_length)
    sources = _dedupe([*base.sources, *vf.sources])
    return GeneMetadata(
        symbol=gene,
        aliases=aliases,
        query_aliases=query_aliases,
        protein_length=protein_length,
        canonical_transcript=vf.canonical_transcript or base.canonical_transcript,
        refseq_transcripts=_dedupe([*base.refseq_transcripts, *vf.refseq_transcripts]),
        ensembl_id=vf.ensembl_id or base.ensembl_id,
        ncbi_id=vf.ncbi_id or base.ncbi_id,
        protein_ids=_dedupe([*base.protein_ids, *vf.protein_ids]),
        sources=sources,
        variantfeatures_db=vf.variantfeatures_db,
    )


def clear_gene_metadata_cache() -> None:
    """Clear cached metadata lookups.

    Tests and one-off scripts that change VARIANTFEATURES_DB at runtime need a
    clean cache so the next lookup sees the new database path. The stored-spelling
    census is keyed by file path, so rewriting a database in place -- which every
    fixture-based test does -- also requires this call.
    """

    get_gene_metadata.cache_clear()
    gene_alias_regex.cache_clear()
    lookup_variantfeatures_residue.cache_clear()
    _stored_gene_spellings.cache_clear()


def get_gene_aliases(
    gene_symbol: str,
    *,
    include_query_aliases: bool = False,
    variantfeatures_db: str | None = None,
) -> tuple[str, ...]:
    metadata = get_gene_metadata(gene_symbol, variantfeatures_db)
    if include_query_aliases:
        return _dedupe([*metadata.aliases, *metadata.query_aliases])
    return metadata.aliases or (metadata.symbol,)


def known_gene_aliases(
    *, include_query_aliases: bool = False
) -> dict[str, tuple[str, ...]]:
    """Return built-in aliases for caption/scope detection."""

    return {
        gene: get_gene_aliases(gene, include_query_aliases=include_query_aliases)
        for gene in BUILTIN_GENE_METADATA
    }


@lru_cache(maxsize=256)
def gene_alias_regex(
    gene_symbol: str,
    *,
    include_query_aliases: bool = False,
    variantfeatures_db: str | None = None,
) -> re.Pattern[str]:
    aliases = get_gene_aliases(
        gene_symbol,
        include_query_aliases=include_query_aliases,
        variantfeatures_db=variantfeatures_db,
    )
    parts = [re.escape(alias).replace(r"\ ", r"\s+") for alias in aliases if alias]
    if not parts:
        parts = [re.escape(normalize_gene_symbol(gene_symbol))]
    return re.compile(
        r"(?<![A-Za-z0-9])(?:" + "|".join(parts) + r")(?![A-Za-z0-9])",
        re.IGNORECASE,
    )


@lru_cache(maxsize=8192)
def lookup_variantfeatures_residue(
    gene_symbol: str,
    *,
    position: int,
    protein_notation: str | None = None,
    cdna_notation: str | None = None,
    variantfeatures_db: str | None = None,
) -> Optional[VariantFeaturesResidue]:
    """Look up residue/transcript support for one protein position.

    Returns None when VariantFeatures is unavailable or the gene/position is not
    represented. The lookup is intentionally read-only and best-effort.

    Memoized because ``pipeline.target_gene_specificity`` calls this once per
    parsed protein position, and papers repeat positions across variants and
    tables. The return value is a frozen dataclass, so callers share it safely.
    Cleared by :func:`clear_gene_metadata_cache`.
    """

    gene = normalize_gene_symbol(gene_symbol)
    db_path = _resolve_variantfeatures_db(variantfeatures_db)
    if not db_path:
        return None
    try:
        with _connect_readonly(db_path) as conn:
            if not _has_table(conn, "variant_consequences"):
                return _lookup_legacy_variantfeatures_residue(
                    conn,
                    gene,
                    position=position,
                    protein_notation=protein_notation,
                    cdna_notation=cdna_notation,
                )
            rows = _gene_rows(
                conn,
                """
                SELECT transcript_id, hgvs_p, hgvs_c, aa_ref, aa_alt
                FROM variant_consequences
                WHERE {gene_predicate} AND aa_pos = ?
                LIMIT 500
                """,
                "variant_consequences",
                "gene_symbol",
                gene,
                (position,),
            )
    except sqlite3.Error:
        return None
    if not rows:
        return None

    protein_norm = _normalize_notation(protein_notation)
    cdna_norm = _normalize_notation(cdna_notation)
    reference_residues = _dedupe(str(row[3]).strip().upper() for row in rows if row[3])
    alternate_residues = _dedupe(str(row[4]).strip().upper() for row in rows if row[4])
    transcripts = _dedupe(str(row[0]).strip() for row in rows if row[0])
    matched_hgvs_p = any(
        protein_norm and _normalize_notation(row[1]).endswith(protein_norm)
        for row in rows
    )
    matched_hgvs_c = any(
        cdna_norm and _normalize_notation(row[2]).endswith(cdna_norm) for row in rows
    )
    return VariantFeaturesResidue(
        gene_symbol=gene,
        position=position,
        reference_residues=reference_residues,
        alternate_residues=alternate_residues,
        transcripts=transcripts,
        matched_hgvs_p=matched_hgvs_p,
        matched_hgvs_c=matched_hgvs_c,
    )


_IDENTIFIER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def resolve_variantfeatures_gene_symbols(
    conn: sqlite3.Connection,
    gene_symbol: str,
    *,
    table: str = "variant_consequences",
    column: str = "gene_symbol",
) -> tuple[str, ...]:
    """Return the stored symbol casing(s) to match with an indexed equality test.

    :func:`_gene_rows` fixes the index-defeat per query, which is right when a
    caller runs one lookup. A caller that runs *several* gene-scoped queries
    against the same connection should resolve the casing once instead: the
    ``UPPER()`` fallback is what costs the full scan, so paying it once on the
    small table beats letting every later query fall back on its own. Callers
    then bind the result with ``column IN (?, ...)``, which SQLite still answers
    from an index such as ``idx_consequences_gene``.

    An empty result means the gene is absent from this database, which lets the
    caller skip its remaining queries outright -- the common case for a gene
    that has no VariantFeatures slice yet.

    ``table``/``column`` are SQL identifiers and cannot be bound as parameters,
    so they are checked against a conservative identifier pattern.
    """

    for identifier in (table, column):
        if not _IDENTIFIER_RE.fullmatch(identifier):
            raise ValueError(f"unsafe SQL identifier: {identifier!r}")

    gene = normalize_gene_symbol(gene_symbol)
    if not gene:
        return ()
    row = conn.execute(
        f"SELECT 1 FROM {table} WHERE {column} = ? LIMIT 1", (gene,)
    ).fetchone()
    if row:
        return (gene,)
    # Mixed-case-only symbols stay reachable: HGNC keeps a lowercase ``orf`` in
    # names such as ``C19orf25``, and real slices carry those rows verbatim.
    return _dedupe_exact(
        str(value).strip()
        for (value,) in conn.execute(
            f"SELECT DISTINCT {column} FROM {table} WHERE UPPER({column}) = ?", (gene,)
        )
        if value
    )


def _dedupe_exact(values: Iterable[str]) -> tuple[str, ...]:
    """Dedupe preserving case, unlike :func:`_dedupe` which folds it.

    Two casings of one symbol are distinct values to match here, so folding them
    would drop rows.
    """

    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return tuple(out)


def _load_variantfeatures_metadata(
    gene: str, variantfeatures_db: str | None = None
) -> Optional[GeneMetadata]:
    db_path = _resolve_variantfeatures_db(variantfeatures_db)
    if not db_path:
        return None
    try:
        with _connect_readonly(db_path) as conn:
            return _read_variantfeatures_metadata(conn, gene, db_path)
    except sqlite3.Error:
        return None


def _resolve_variantfeatures_db(variantfeatures_db: str | None) -> Optional[Path]:
    if variantfeatures_db:
        path = Path(variantfeatures_db).expanduser()
        return path if path.is_file() else None
    return default_variantfeatures_db_path()


def _connect_readonly(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def _read_variantfeatures_metadata(
    conn: sqlite3.Connection, gene: str, db_path: Path
) -> Optional[GeneMetadata]:
    aliases: list[str] = [gene]
    refseq_transcripts: list[str] = []
    protein_ids: list[str] = []
    canonical_transcript: str | None = None
    ensembl_id: str | None = None
    ncbi_id: str | None = None
    protein_length: int | None = None

    if _has_table(conn, "genes"):
        row = _first_gene_row(
            conn,
            "SELECT symbol, canonical_transcript, ncbi_id, ensembl_id FROM genes WHERE {gene_predicate}",
            "genes",
            "symbol",
            gene,
        )
        if row:
            aliases.append(str(row[0]))
            canonical_transcript = row[1]
            ncbi_id = row[2]
            ensembl_id = row[3]

    if _has_table(conn, "transcripts"):
        rows = _gene_rows(
            conn,
            """
            SELECT transcript_id, refseq_match, protein_id, cds_length
            FROM transcripts
            WHERE {gene_predicate}
            ORDER BY is_mane_select DESC, is_canonical DESC, transcript_id
            LIMIT 20
            """,
            "transcripts",
            "gene_symbol",
            gene,
        )
        for transcript_id, refseq_match, protein_id, cds_length in rows:
            if transcript_id and not canonical_transcript:
                canonical_transcript = str(transcript_id)
            if refseq_match:
                refseq_transcripts.append(str(refseq_match))
            if protein_id:
                protein_ids.append(str(protein_id))
            derived = _protein_length_from_cds(cds_length)
            if derived and not protein_length:
                protein_length = derived

    if _has_table(conn, "variant_consequences"):
        row = _first_gene_row(
            conn,
            """
            SELECT MAX(aa_pos)
            FROM variant_consequences
            WHERE {gene_predicate}
            """,
            "variant_consequences",
            "gene_symbol",
            gene,
        )
        if row and row[0]:
            # Saturation tables include the stop codon as the final aa_pos.
            protein_length = max(int(row[0]) - 1, 1)
        rows = _gene_rows(
            conn,
            """
            SELECT DISTINCT transcript_id
            FROM variant_consequences
            WHERE {gene_predicate} AND transcript_id IS NOT NULL
            LIMIT 20
            """,
            "variant_consequences",
            "gene_symbol",
            gene,
        )
        for (transcript_id,) in rows:
            if transcript_id and not canonical_transcript:
                canonical_transcript = str(transcript_id)

    elif _has_table(conn, "variants"):
        row = _first_gene_row(
            conn,
            "SELECT MAX(resnum) FROM variants WHERE {gene_predicate}",
            "variants",
            "gene",
            gene,
        )
        if row and row[0]:
            protein_length = int(row[0])
        rows = _gene_rows(
            conn,
            "SELECT DISTINCT uniprot_id FROM variants WHERE {gene_predicate} AND uniprot_id IS NOT NULL LIMIT 20",
            "variants",
            "gene",
            gene,
        )
        protein_ids.extend(str(row[0]) for row in rows if row[0])

    if not any(
        [canonical_transcript, ensembl_id, ncbi_id, protein_length, protein_ids]
    ):
        return None

    return GeneMetadata(
        symbol=gene,
        aliases=_dedupe(aliases),
        query_aliases=_dedupe(aliases),
        protein_length=protein_length,
        canonical_transcript=canonical_transcript,
        refseq_transcripts=_dedupe(refseq_transcripts),
        ensembl_id=ensembl_id,
        ncbi_id=ncbi_id,
        protein_ids=_dedupe(protein_ids),
        sources=("variantfeatures",),
        variantfeatures_db=str(db_path),
    )


def _lookup_legacy_variantfeatures_residue(
    conn: sqlite3.Connection,
    gene: str,
    *,
    position: int,
    protein_notation: str | None,
    cdna_notation: str | None,
) -> Optional[VariantFeaturesResidue]:
    if not _has_table(conn, "variants"):
        return None
    rows = _gene_rows(
        conn,
        """
        SELECT var, var_hgvs_p, var_hgvs_c, wt_aa, mut_aa, uniprot_id
        FROM variants
        WHERE {gene_predicate} AND resnum = ?
        LIMIT 500
        """,
        "variants",
        "gene",
        gene,
        (position,),
    )
    if not rows:
        return None
    protein_norm = _normalize_notation(protein_notation)
    cdna_norm = _normalize_notation(cdna_notation)
    return VariantFeaturesResidue(
        gene_symbol=gene,
        position=position,
        reference_residues=_dedupe(
            str(row[3]).strip().upper() for row in rows if row[3]
        ),
        alternate_residues=_dedupe(
            str(row[4]).strip().upper() for row in rows if row[4]
        ),
        transcripts=_dedupe(str(row[5]).strip() for row in rows if row[5]),
        matched_hgvs_p=any(
            protein_norm
            and (
                _normalize_notation(row[0]).endswith(protein_norm)
                or _normalize_notation(row[1]).endswith(protein_norm)
            )
            for row in rows
        ),
        matched_hgvs_c=any(
            cdna_norm and _normalize_notation(row[2]).endswith(cdna_norm)
            for row in rows
        ),
    )


def _choose_protein_length(
    builtin_length: Optional[int], variantfeatures_length: Optional[int]
) -> Optional[int]:
    if not builtin_length:
        return variantfeatures_length
    if not variantfeatures_length:
        return builtin_length
    if variantfeatures_length >= int(builtin_length * 0.8):
        return variantfeatures_length
    return builtin_length


def _protein_length_from_cds(cds_length: object) -> Optional[int]:
    try:
        length = int(cds_length)
    except (TypeError, ValueError):
        return None
    if length <= 0:
        return None
    return max(length // 3 - 1, 1)


def _has_table(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE name = ? AND type IN ('table', 'view')",
        (name,),
    ).fetchone()
    return bool(row)


def _gene_rows(
    conn: sqlite3.Connection,
    sql: str,
    table: str,
    column: str,
    gene: str,
    extra_params: tuple[object, ...] = (),
) -> list[tuple]:
    """Run a gene-scoped query using only index-friendly equality predicates.

    ``sql`` must contain a single ``{gene_predicate}`` placeholder. ``gene`` is
    already uppercased by :func:`normalize_gene_symbol`, so the exact form
    matches whenever the database stores uppercase symbols and can then use an
    index such as ``idx_consequences_gene``. Wrapping the column in ``UPPER()``
    makes the column non-indexable and degrades every lookup into a full index
    scan, which costs seconds per query on a multi-gigabyte VariantFeatures
    slice.

    Genes stored only in mixed case must still resolve -- HGNC keeps a lowercase
    ``orf`` in symbols such as ``C19orf25``, and real databases carry those rows
    verbatim. Rather than reach for ``UPPER()``, the exact-match miss consults a
    cached census of the spellings the column actually stores
    (:func:`_stored_gene_spellings`) and re-queries those spellings by equality.
    Because the census enumerates every stored value, ``column IN (<spellings>)``
    selects exactly the rows ``UPPER(column) = ?`` would have, and a gene with no
    stored spelling at all resolves to an in-memory miss with no second query --
    which is the whole point, since a novel gene is by definition absent and the
    scan used to be paid once per protein position.

    An aggregate such as ``SELECT MAX(aa_pos)`` always returns one row, so an
    empty result cannot drive the fallback on its own; a row of nothing but
    ``NULL`` is treated as "no match" for that reason.
    """

    rows = conn.execute(
        sql.format(gene_predicate=f"{column} = ?"), (gene, *extra_params)
    ).fetchall()
    if _has_payload(rows):
        return rows

    spellings = _upper_matched_spellings(conn, table, column, gene)
    if spellings is None:
        # No filename to key a census on (an in-memory or temporary database).
        # Keep the original scanning fallback rather than assume the gene is
        # absent: guessing wrong here silently zeroes a gene's metadata.
        return conn.execute(
            sql.format(gene_predicate=f"UPPER({column}) = ?"), (gene, *extra_params)
        ).fetchall()
    if not spellings or spellings == (gene,):
        # Either nothing is stored under this symbol, or the only stored spelling
        # is the one the exact probe already tried. An aggregate over an empty
        # set yields the same all-NULL row for any predicate, so the rows from
        # that probe are what the scanning fallback would have returned.
        return rows

    placeholders = ", ".join("?" * len(spellings))
    return conn.execute(
        sql.format(gene_predicate=f"{column} IN ({placeholders})"),
        (*spellings, *extra_params),
    ).fetchall()


def _first_gene_row(
    conn: sqlite3.Connection,
    sql: str,
    table: str,
    column: str,
    gene: str,
    extra_params: tuple[object, ...] = (),
) -> Optional[tuple]:
    rows = _gene_rows(conn, sql, table, column, gene, extra_params)
    return rows[0] if rows else None


def _has_payload(rows: list[tuple]) -> bool:
    """True when ``rows`` carry at least one non-NULL value.

    ``SELECT MAX(aa_pos)`` always returns one row, so ``if not rows`` cannot
    distinguish "matched nothing" from "matched rows whose columns are NULL".
    """

    return any(value is not None for row in rows for value in row)


def _upper_matched_spellings(
    conn: sqlite3.Connection, table: str, column: str, gene: str
) -> Optional[tuple[str, ...]]:
    """Stored spellings of ``gene``, or ``None`` when no census is possible.

    An empty tuple means the census ran and the gene is genuinely absent.
    """

    db_path = _connection_filename(conn)
    if not db_path:
        return None
    return _stored_gene_spellings(db_path, table, column).get(gene, ())


def _connection_filename(conn: sqlite3.Connection) -> str:
    """Path backing the ``main`` schema, or ``""`` for in-memory databases."""

    try:
        for _, name, filename in conn.execute("PRAGMA database_list"):
            if name == "main":
                return filename or ""
    except sqlite3.Error:
        return ""
    return ""


@lru_cache(maxsize=32)
def _stored_gene_spellings(
    db_path: str, table: str, column: str
) -> Mapping[str, tuple[str, ...]]:
    """Map SQLite-``UPPER()`` of each stored gene symbol to its spellings.

    Built lazily -- only an exact-match miss asks for it -- so a database that
    stores uppercase symbols never pays for the census. On the real 38 GB slice
    ``variant_consequences`` holds 3.47M rows but only 200 distinct symbols, and
    the walk in :func:`_distinct_column_values` reaches them in 0.034s cold /
    0.001s warm.

    Cached for the life of the process and dropped by
    :func:`clear_gene_metadata_cache`; a database rewritten in place without that
    call will keep answering from the old census.
    """

    groups: dict[str, list[str]] = {}
    try:
        with _connect_readonly(Path(db_path)) as conn:
            for value in _distinct_column_values(conn, table, column):
                groups.setdefault(_sqlite_upper(value), []).append(value)
    except sqlite3.Error:
        return MappingProxyType({})
    return MappingProxyType({key: tuple(values) for key, values in groups.items()})


def _distinct_column_values(
    conn: sqlite3.Connection, table: str, column: str
) -> list[str]:
    """Enumerate the distinct non-NULL values of one column."""

    if _has_leading_index(conn, table, column):
        # Loose index scan: one index seek per distinct value instead of walking
        # every row. On variant_consequences that is ~200 seeks rather than a
        # 3.47M-entry covering-index scan (0.034s versus 2.4-3.4s).
        sql = f"""
            WITH RECURSIVE walk(value) AS (
                SELECT MIN({column}) FROM {table}
                UNION ALL
                SELECT (
                    SELECT MIN({column}) FROM {table} WHERE {column} > walk.value
                )
                FROM walk
                WHERE walk.value IS NOT NULL
            )
            SELECT value FROM walk WHERE value IS NOT NULL
        """
    else:
        # Without an index the walk would degrade to one full scan per distinct
        # value, so take a single pass instead. Still amortized: this replaces a
        # scan paid on every lookup with one paid once.
        sql = f"SELECT DISTINCT {column} FROM {table} WHERE {column} IS NOT NULL"
    return [str(row[0]) for row in conn.execute(sql) if row[0] is not None]


def _has_leading_index(conn: sqlite3.Connection, table: str, column: str) -> bool:
    """True when some complete index has ``column`` as its first column.

    Partial indexes do not count: SQLite only uses one when the query implies its
    WHERE clause, which the unconstrained ``MIN()`` walk never does, so the walk
    would silently fall back to a full scan per step.
    """

    try:
        indexes = conn.execute(f"PRAGMA index_list({table})").fetchall()
    except sqlite3.Error:
        return False
    for row in indexes:
        name = row[1]
        partial = row[4] if len(row) > 4 else 0
        if partial:
            continue
        info = conn.execute(f"PRAGMA index_info({name})").fetchall()
        if info and info[0][2] == column:
            return True
    return False


_ASCII_UPPER = str.maketrans(string.ascii_lowercase, string.ascii_uppercase)


def _sqlite_upper(value: str) -> str:
    """Mirror SQLite's ``upper()``, which only folds ASCII.

    ``normalize_gene_symbol`` uses Python's Unicode-aware ``str.upper()``, so
    matching the census with that instead would diverge from the predicate this
    replaces on any non-ASCII symbol.
    """

    return value.translate(_ASCII_UPPER)


def _dedupe(values: Iterable[object]) -> tuple[str, ...]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = str(value).strip()
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return tuple(out)


def _normalize_notation(value: object) -> str:
    text = str(value or "").strip()
    if ":" in text:
        text = text.rsplit(":", 1)[-1]
    return re.sub(r"\s+", "", text).lower()
