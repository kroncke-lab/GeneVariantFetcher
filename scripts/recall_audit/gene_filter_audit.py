#!/usr/bin/env python3
"""Flag PMIDs whose extracted rows look off-target for a single-gene run."""

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

try:
    from scripts.recall_audit.common import (
        query_pipeline_rows,
        resolve_gene_db,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import query_pipeline_rows, resolve_gene_db, write_csv_rows

GENE_TOKEN_RE = re.compile(r"\b[A-Z][A-Z0-9-]{2,11}\b")
IGNORE_TOKENS = {
    "ALL",
    "DNA",
    "RNA",
    "ECG",
    "QT",
    "QTC",
    "LQT",
    "LQTS",
    "CPVT",
    "BRGDA",
    "SCD",
    "WT",
    "NA",
    "PMID",
    "GRCH37",
    "GRCH38",
    "HG19",
    "HG38",
    "ACMG",
    "ARVC",
    "BRUGADA",
    "CICR",
    "CRDS",
    "CTD",
    "EXPERIMENTAL",
    "GOF",
    "HD1",
    "HUMAN",
    "ISO",
    "LOF",
    "MAF",
    "MAPP",
    "MODELS",
    "MVP",
    "NTD",
    "SIFT",
    "SNV",
    "SNP",
    "SOICR",
    "SUD",
    "SUDY",
    "TCH",
    "VSC",
    "VUS",
    "WES",
}
PROTEIN_LIKE_RE = re.compile(
    r"^(?:p\.)?(?:(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
    r"Phe|Pro|Ser|Thr|Trp|Tyr|Val)|[ACDEFGHIKLMNPQRSTVWY])\d{1,5}"
    r"(?:[A-Z*?]|[a-z]{2}|fs|del|dup|ins)",
    re.IGNORECASE,
)


def looks_like_gene_symbol(token: str, target_tokens: set[str]) -> bool:
    token = token.strip().upper()
    if token in target_tokens or token in IGNORE_TOKENS:
        return False
    if (
        token.startswith("IVS")
        or re.match(r"^[A-Z]\d", token)
        or PROTEIN_LIKE_RE.match(token)
    ):
        return False
    letters = sum(ch.isalpha() for ch in token)
    return letters >= 2 and bool(GENE_TOKEN_RE.fullmatch(token))


def target_tokens(gene: str) -> set[str]:
    tokens = {gene.upper()}
    # Common aliases that otherwise look like off-target symbols.
    if gene.upper() == "KCNH2":
        tokens.update({"HERG", "HERG1", "ERG", "ERG1"})
    return tokens


def extract_gene_tokens(text: str, target: set[str]) -> set[str]:
    return {
        token.upper()
        for token in GENE_TOKEN_RE.findall(text or "")
        if looks_like_gene_symbol(token, target)
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--db", help="Explicit SQLite DB path")
    parser.add_argument("--run-dir", help="Validation run directory containing dbs/")
    parser.add_argument("--min-off-target-rows", type=int, default=5)
    parser.add_argument("--min-distinct-symbols", type=int, default=3)
    parser.add_argument("--out", help="Write CSV here instead of stdout")
    args = parser.parse_args()

    gene = args.gene.upper()
    rows = query_pipeline_rows(resolve_gene_db(gene, args.db, args.run_dir))
    target = target_tokens(gene)

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get("pmid") or "")].append(row)

    flagged: list[dict[str, object]] = []
    for pmid, pmid_rows in grouped.items():
        token_counts: Counter[str] = Counter()
        examples: list[str] = []
        off_target_rows = 0
        for row in pmid_rows:
            row_text = " ".join(
                str(row.get(key) or "")
                for key in (
                    "gene_symbol",
                    "cdna_notation",
                    "protein_notation",
                    "genomic_position",
                    "source_location",
                    "additional_notes",
                )
            )
            tokens = extract_gene_tokens(row_text, target)
            if tokens:
                off_target_rows += 1
                token_counts.update(tokens)
                if len(examples) < 3:
                    examples.append(
                        f"{row.get('protein_notation') or row.get('cdna_notation')}: {','.join(sorted(tokens))}"
                    )

        if (
            off_target_rows >= args.min_off_target_rows
            or len(token_counts) >= args.min_distinct_symbols
        ):
            flagged.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "title": pmid_rows[0].get("title") or "",
                    "total_rows": len(pmid_rows),
                    "off_target_rows": off_target_rows,
                    "distinct_off_target_symbols": len(token_counts),
                    "top_symbols": ";".join(
                        f"{token}:{count}"
                        for token, count in token_counts.most_common(10)
                    ),
                    "examples": " | ".join(examples),
                }
            )

    flagged.sort(
        key=lambda row: (
            -int(row["off_target_rows"]),
            -int(row["total_rows"]),
            row["pmid"],
        )
    )
    write_csv_rows(
        flagged,
        [
            "gene",
            "pmid",
            "title",
            "total_rows",
            "off_target_rows",
            "distinct_off_target_symbols",
            "top_symbols",
            "examples",
        ],
        Path(args.out) if args.out else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
