"""Shared gene metadata for discovery, triage, and specificity QC.

The pipeline should not need one-off alias/protein-length maps in every stage.
This module centralizes conservative built-ins and augments them with a local
VariantFeatures SQLite database when available.
"""

from __future__ import annotations

import os
import re
import sqlite3
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Optional


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
}


def normalize_gene_symbol(gene_symbol: str) -> str:
    return (gene_symbol or "").strip().upper()


def default_variantfeatures_db_path() -> Optional[Path]:
    """Return a likely VariantFeatures SQLite path, if configured or present."""

    for key in ("VARIANTFEATURES_DB", "VARIANT_FEATURES_DB"):
        value = os.environ.get(key)
        if value and Path(value).expanduser().is_file():
            return Path(value).expanduser()

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
    clean cache so the next lookup sees the new database path.
    """

    get_gene_metadata.cache_clear()
    gene_alias_regex.cache_clear()


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
            rows = conn.execute(
                """
                SELECT transcript_id, hgvs_p, hgvs_c, aa_ref, aa_alt
                FROM variant_consequences
                WHERE UPPER(gene_symbol) = ? AND aa_pos = ?
                LIMIT 500
                """,
                (gene, position),
            ).fetchall()
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
        row = _first_row(
            conn,
            "SELECT symbol, canonical_transcript, ncbi_id, ensembl_id FROM genes WHERE UPPER(symbol) = ?",
            (gene,),
        )
        if row:
            aliases.append(str(row[0]))
            canonical_transcript = row[1]
            ncbi_id = row[2]
            ensembl_id = row[3]

    if _has_table(conn, "transcripts"):
        rows = conn.execute(
            """
            SELECT transcript_id, refseq_match, protein_id, cds_length
            FROM transcripts
            WHERE UPPER(gene_symbol) = ?
            ORDER BY is_mane_select DESC, is_canonical DESC, transcript_id
            LIMIT 20
            """,
            (gene,),
        ).fetchall()
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
        row = _first_row(
            conn,
            """
            SELECT MAX(aa_pos)
            FROM variant_consequences
            WHERE UPPER(gene_symbol) = ?
            """,
            (gene,),
        )
        if row and row[0]:
            # Saturation tables include the stop codon as the final aa_pos.
            protein_length = max(int(row[0]) - 1, 1)
        rows = conn.execute(
            """
            SELECT DISTINCT transcript_id
            FROM variant_consequences
            WHERE UPPER(gene_symbol) = ? AND transcript_id IS NOT NULL
            LIMIT 20
            """,
            (gene,),
        ).fetchall()
        for (transcript_id,) in rows:
            if transcript_id and not canonical_transcript:
                canonical_transcript = str(transcript_id)

    elif _has_table(conn, "variants"):
        row = _first_row(
            conn,
            "SELECT MAX(resnum) FROM variants WHERE UPPER(gene) = ?",
            (gene,),
        )
        if row and row[0]:
            protein_length = int(row[0])
        rows = conn.execute(
            "SELECT DISTINCT uniprot_id FROM variants WHERE UPPER(gene) = ? AND uniprot_id IS NOT NULL LIMIT 20",
            (gene,),
        ).fetchall()
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
    rows = conn.execute(
        """
        SELECT var, var_hgvs_p, var_hgvs_c, wt_aa, mut_aa, uniprot_id
        FROM variants
        WHERE UPPER(gene) = ? AND resnum = ?
        LIMIT 500
        """,
        (gene, position),
    ).fetchall()
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


def _first_row(
    conn: sqlite3.Connection, sql: str, params: tuple[object, ...]
) -> Optional[tuple]:
    return conn.execute(sql, params).fetchone()


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
