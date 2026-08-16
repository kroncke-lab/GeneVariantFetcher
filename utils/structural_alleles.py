"""Generic structural-allele equivalence.

Maps literature spellings (``EXON 3 DELETION``, ``ΔKPQ``, ``INS G 1682``) onto
protein-range events using the gene's cached RefSeq protein and the NCBI exon
map for its canonical transcript. No per-variant or per-gene special cases.

Fail closed when the gene has no map/sequence, an exon span is not an in-frame
coding block, or a delta peptide is missing or occurs more than once.
"""

from __future__ import annotations

import json
import logging
import re
from functools import lru_cache
from importlib.resources import files
from typing import Any, Iterable, Optional

logger = logging.getLogger(__name__)

_AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

_EXON_EVENT_RE = re.compile(
    r"(?:(?P<op>del(?:etion)?|dup(?:lication)?)\s+(?:of\s+)?exons?\s*"
    r"(?P<e1>\d+)(?:\s*(?:[-–—to]+)\s*(?P<e2>\d+))?"
    r"|exons?\s*(?P<e1b>\d+)(?:\s*(?:[-–—to]+)\s*(?P<e2b>\d+))?\s*"
    r"(?P<opb>del(?:etion)?|dup(?:lication)?)"
    r"|(?P<opc>del|dup):exon(?P<e1c>\d+)(?:-(?P<e2c>\d+))?)",
    re.IGNORECASE,
)
_DELTA_RE = re.compile(
    r"(?:Δ|delta)\s*-?\s*([A-Z]{2,}|[A-Z][a-z]{2}(?:[A-Z][a-z]{2})+)",
    re.IGNORECASE,
)
_RANGE_DEL_RE = re.compile(
    r"(?:p\.)?(?:(?P<r1>[A-Z][a-z]{2})|(?P<s1>[A-Z]))(?P<p1>\d+)"
    r"[_-](?:(?P<r2>[A-Z][a-z]{2})|(?P<s2>[A-Z]))(?P<p2>\d+)del",
    re.IGNORECASE,
)
_INS_POS_RE = re.compile(
    r"\bins\s+[ACGT]+\s+(\d+)\b",
    re.IGNORECASE,
)


def parse_ncbi_feature_table(
    text: str,
) -> tuple[list[tuple[int, int]], Optional[tuple[int, int]]]:
    """Return ``(exon_spans, cds_span)`` from an NCBI ``rettype=ft`` table."""
    exons: list[tuple[int, int]] = []
    cds: Optional[tuple[int, int]] = None
    for raw in text.splitlines():
        if raw.startswith(">") or raw.startswith("\t"):
            continue
        parts = raw.split("\t")
        if len(parts) < 3 or not parts[0] or not parts[1]:
            continue
        try:
            start, end = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        kind = parts[2]
        if kind == "exon":
            exons.append((start, end))
        elif kind == "CDS" and cds is None:
            cds = (start, end)
    return exons, cds


def build_exon_records(
    exons: list[tuple[int, int]],
    cds: tuple[int, int],
    sequence: str,
) -> list[dict[str, Any]]:
    """Intersect transcript exons with CDS and attach in-frame protein spans."""
    cds_start, cds_end = cds
    protein_nt = 3 * len(sequence)
    coding_offset = 0
    records: list[dict[str, Any]] = []
    for idx, (estart, eend) in enumerate(exons, start=1):
        ov_start = max(estart, cds_start)
        ov_end = min(eend, cds_end)
        if ov_start > ov_end:
            records.append(
                {
                    "exon": idx,
                    "transcript_start": estart,
                    "transcript_end": eend,
                    "coding_nt": 0,
                    "protein_nt": 0,
                    "codon_aligned": False,
                    "in_frame": False,
                }
            )
            continue
        n_nt = ov_end - ov_start + 1
        protein_here = max(0, min(n_nt, protein_nt - coding_offset))
        aligned = coding_offset % 3 == 0
        in_frame = aligned and protein_here > 0 and protein_here % 3 == 0
        rec: dict[str, Any] = {
            "exon": idx,
            "transcript_start": estart,
            "transcript_end": eend,
            "coding_nt": n_nt,
            "protein_nt": protein_here,
            "cds_nt_start": coding_offset + 1,
            "cds_nt_end": coding_offset + n_nt,
            "codon_aligned": aligned,
            "in_frame": in_frame,
        }
        if in_frame:
            aa_start = coding_offset // 3 + 1
            aa_end = (coding_offset + protein_here) // 3
            rec["aa_start"] = aa_start
            rec["aa_end"] = aa_end
            rec["aa_start_res"] = sequence[aa_start - 1]
            rec["aa_end_res"] = sequence[aa_end - 1]
            rec["protein_del"] = protein_range_key(sequence, aa_start, aa_end)
        records.append(rec)
        coding_offset += n_nt
    return records


def protein_range_key(sequence: str, aa_start: int, aa_end: int) -> str:
    return f"{sequence[aa_start - 1]}{aa_start}_{sequence[aa_end - 1]}{aa_end}del"


@lru_cache(maxsize=32)
def load_exon_maps() -> dict[str, Any]:
    path = files("gvf_data").joinpath("reference_sequences/exon_maps.json")
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("exon map load failed: %s", exc)
        return {}
    genes = payload.get("genes") if isinstance(payload, dict) else None
    return genes if isinstance(genes, dict) else {}


def exon_map_for(gene: str) -> Optional[dict[str, Any]]:
    if not gene:
        return None
    return load_exon_maps().get(gene.upper())


def _aa1(token: Optional[str]) -> Optional[str]:
    if not token:
        return None
    if len(token) == 1:
        return token.upper()
    return _AA3_TO_1.get(token.upper())


def parse_exon_event(text: str) -> Optional[tuple[str, int, int]]:
    """Return ``(del|dup, first_exon, last_exon)`` or None."""
    if not text:
        return None
    match = _EXON_EVENT_RE.search(str(text))
    if not match:
        return None
    op = (match.group("op") or match.group("opb") or match.group("opc") or "").lower()
    op = "dup" if op.startswith("dup") else "del"
    start = match.group("e1") or match.group("e1b") or match.group("e1c")
    end = match.group("e2") or match.group("e2b") or match.group("e2c") or start
    if not start:
        return None
    first, last = int(start), int(end)
    if first > last:
        first, last = last, first
    return op, first, last


def parse_delta_peptide(text: str) -> Optional[str]:
    if not text:
        return None
    match = _DELTA_RE.search(str(text))
    if not match:
        return None
    raw = match.group(1)
    if raw.isupper() and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", raw):
        return raw
    residues = re.findall(r"[A-Z][a-z]{2}", raw)
    letters = [_AA3_TO_1.get(res.upper()) for res in residues]
    if letters and all(letters):
        return "".join(letters)  # type: ignore[arg-type]
    return None


def parse_protein_range_del(text: str) -> Optional[tuple[int, int]]:
    if not text:
        return None
    match = _RANGE_DEL_RE.search(str(text).replace(" ", ""))
    if not match:
        return None
    start = int(match.group("p1"))
    end = int(match.group("p2"))
    if start > end:
        start, end = end, start
    return start, end


def parse_ins_position(text: str) -> Optional[int]:
    if not text:
        return None
    match = _INS_POS_RE.search(str(text))
    return int(match.group(1)) if match else None


def _exon_protein_span(
    records: Iterable[dict[str, Any]], first: int, last: int
) -> Optional[tuple[int, int]]:
    chosen = [rec for rec in records if first <= int(rec.get("exon") or 0) <= last]
    chosen.sort(key=lambda rec: int(rec["exon"]))
    if not chosen or int(chosen[0]["exon"]) != first or int(chosen[-1]["exon"]) != last:
        return None
    if any(not rec.get("in_frame") for rec in chosen):
        return None
    expected = int(chosen[0]["aa_start"])
    for rec in chosen:
        if int(rec["aa_start"]) != expected:
            return None
        expected = int(rec["aa_end"]) + 1
    return int(chosen[0]["aa_start"]), int(chosen[-1]["aa_end"])


def _unique_peptide_span(sequence: str, peptide: str) -> Optional[tuple[int, int]]:
    if not sequence or not peptide:
        return None
    hits = [
        idx
        for idx in range(len(sequence) - len(peptide) + 1)
        if sequence[idx : idx + len(peptide)] == peptide
    ]
    if len(hits) != 1:
        return None
    start = hits[0] + 1
    return start, start + len(peptide) - 1


def _iter_in_frame_spans(
    records: Iterable[dict[str, Any]],
) -> Iterable[tuple[int, int, int, int]]:
    """Yield ``(first_exon, last_exon, aa_start, aa_end)`` for contiguous in-frame runs."""
    run: list[dict[str, Any]] = []
    for rec in sorted(records, key=lambda item: int(item.get("exon") or 0)):
        if rec.get("in_frame"):
            if run and int(rec["aa_start"]) != int(run[-1]["aa_end"]) + 1:
                run = []
            run.append(rec)
            for i in range(len(run)):
                yield (
                    int(run[i]["exon"]),
                    int(rec["exon"]),
                    int(run[i]["aa_start"]),
                    int(rec["aa_end"]),
                )
        else:
            run = []


def expand_structural_keys(variant: str, gene: str) -> set[str]:
    """Return comparable keys for a structural / named-indel spelling."""
    if not variant or not gene:
        return set()
    keys: set[str] = set()
    gene_map = exon_map_for(gene)
    records = list((gene_map or {}).get("exons") or [])

    from pipeline.reference_validation import load_reference_protein

    sequence = load_reference_protein(gene)

    event = parse_exon_event(variant)
    if event:
        op, first, last = event
        suffix = f"{first}-{last}" if first != last else str(first)
        keys.add(f"{op}:exon{suffix}")
        keys.add(f"EXON{first}{f'_{last}' if first != last else ''}{op.upper()}")
        span = _exon_protein_span(records, first, last) if records else None
        if span and op == "del" and sequence:
            keys.add(protein_range_key(sequence, *span))

    peptide = parse_delta_peptide(variant)
    if peptide:
        keys.add(f"DELTA{peptide}")
        span = _unique_peptide_span(sequence or "", peptide) if sequence else None
        if span and sequence:
            keys.add(protein_range_key(sequence, *span))

    range_del = parse_protein_range_del(variant)
    if range_del:
        start, end = range_del
        if sequence and 1 <= start <= end <= len(sequence):
            keys.add(protein_range_key(sequence, start, end))
        else:
            keys.add(f"{start}_{end}del")
        for first, last, aa_start, aa_end in _iter_in_frame_spans(records):
            if (aa_start, aa_end) == (start, end):
                suffix = f"{first}_{last}" if first != last else str(first)
                keys.add(f"EXON{suffix}DEL")
                keys.add(f"del:exon{first}{f'-{last}' if first != last else ''}")

    ins_pos = parse_ins_position(variant)
    if ins_pos and sequence and 1 <= ins_pos <= len(sequence):
        keys.add(f"{sequence[ins_pos - 1]}{ins_pos}ins")

    return {key for key in keys if key}


def structural_keys_match(left: str, right: str, gene: str) -> bool:
    if not left or not right or not gene:
        return False
    return bool(
        expand_structural_keys(left, gene) & expand_structural_keys(right, gene)
    )
