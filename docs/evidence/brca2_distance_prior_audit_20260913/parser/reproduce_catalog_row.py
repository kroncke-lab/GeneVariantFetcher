"""Read-only reproduction of a catalog row becoming an affected carrier.

No model calls, source writes, or production extraction run are performed.
Controls alter only the in-memory snippet to isolate header/placeholder effects.
"""

import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
os.environ["LITELLM_LOCAL_MODEL_COST_MAP"] = "True"
sys.path.insert(0, str(REPO))
from pipeline.extraction import ExpertExtractor
from pipeline.table_router import (
    _infer_column_mapping_from_headers,
    _looks_like_row_level_clinical_list,
    _split_pipe_row,
    enumerate_markdown_tables,
)


def pipe_row(cells):
    return "| " + " | ".join(cells) + " |"


def exercise(extractor, header, separator, row):
    text = "\n".join([header, separator, row]) + "\n"
    cells = _split_pipe_row(header)
    records = extractor._parse_markdown_table_variants(text, "BRCA2")
    keep = [
        "protein_notation",
        "clinical_significance",
        "penetrance_data",
        "count_provenance",
    ]
    tables = enumerate_markdown_tables(text)
    assert len(tables) == 1
    mapping = _infer_column_mapping_from_headers(tables[0], target_gene="BRCA2")
    prototype = {
        "protein": cells.index("Protein change"),
        "gene": cells.index("Gene(s)"),
    }
    return {
        "header_cells": len(cells),
        "row_cells": len(_split_pipe_row(row)),
        "has_row_subject_header": any(
            any(term in cell.lower() for term in extractor.VARIANT_ROW_LEVEL_HEADERS)
            for cell in cells
        ),
        "clinical_header_triggers": [
            [cell, term]
            for cell in cells
            for term in extractor.VARIANT_CLINICAL_CONTEXT_HEADERS
            if term in cell.lower()
        ],
        "legacy_row_level_clinical_header": extractor._looks_like_row_level_clinical_header(
            cells
        ),
        "legacy_parser_output": [
            {key: value[key] for key in keep} for value in records
        ],
        "router_deterministic_mapping": mapping,
        "router_row_level_test_with_notation_mapping": _looks_like_row_level_clinical_list(
            tables[0], prototype, False
        ),
    }


def main():
    source_dir = (
        REPO / "results/grant_e2e_20260909/shards/BRCA2_0/BRCA2/20260909_115038"
    )
    source_path = source_dir / "pmc_fulltext/40664060_CLEANED.md"
    full_text = source_path.read_text()
    lines = full_text.splitlines()
    target = next(i for i, row in enumerate(lines) if "A1043V" in row)
    header = max(i for i in range(target) if lines[i].startswith("| Name | Gene(s) |"))
    original = (lines[header], lines[header + 1], lines[target])
    (HERE / "source_excerpt.md").write_text("\n".join(original) + "\n")
    extractor = ExpertExtractor.__new__(ExpertExtractor)
    results = {"exact_source": exercise(extractor, *original)}
    aligned_row = original[2].replace(
        "Hereditary breast ovarian cancer syndrome|Hereditary cancer-predisposing syndrome",
        "Hereditary breast ovarian cancer syndrome;Hereditary cancer-predisposing syndrome",
    )
    results["condition_delimiter_repaired_control"] = exercise(
        extractor, original[0], original[1], aligned_row
    )
    headers, cells = _split_pipe_row(original[0]), _split_pipe_row(aligned_row)
    assert len(headers) == len(cells)
    clinical = [i for i, name in enumerate(headers) if "clinical" in name.lower()]
    known_missing = cells.copy()
    for i in clinical:
        known_missing[i] = "unknown"
    results["clinical_nan_to_unknown_control"] = exercise(
        extractor, original[0], original[1], pipe_row(known_missing)
    )
    without_annotation = [i for i in range(len(headers)) if i not in clinical]
    results["somatic_clinical_annotation_columns_removed_control"] = exercise(
        extractor,
        pipe_row([headers[i] for i in without_annotation]),
        pipe_row(["---"] * len(without_annotation)),
        pipe_row([cells[i] for i in without_annotation]),
    )
    exact = results["exact_source"]["legacy_parser_output"]
    assert len(exact) == 1 and exact[0]["protein_notation"] == "A1043V"
    assert exact[0]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 1,
        "unaffected_count": None,
    }
    with sqlite3.connect(f"file:{source_dir / 'BRCA2.db'}?mode=ro", uri=True) as conn:
        conn.row_factory = sqlite3.Row
        archived = dict(
            conn.execute(
                "SELECT * FROM variant_papers WHERE variant_id=? AND pmid=?",
                (12496, "40664060"),
            ).fetchone()
        )
        facts = [
            dict(row)
            for row in conn.execute(
                "SELECT fact_type,fact_value,count_type,source_layer,provenance_kind FROM fact_provenance WHERE variant_id=? AND pmid=?",
                (12496, "40664060"),
            )
        ]
    receipt = {
        "pmid": "40664060",
        "archive_variant_id": 12496,
        "source_file": str(source_path.relative_to(REPO)),
        "source_sha256": hashlib.sha256(full_text.encode()).hexdigest(),
        "header_line": header + 1,
        "target_line": target + 1,
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "source_files_modified": subprocess.check_output(
            [
                "git",
                "status",
                "--porcelain",
                "--",
                "pipeline/extraction.py",
                "pipeline/table_router.py",
            ],
            cwd=REPO,
            text=True,
        ).strip(),
        "current_source_hashes": {
            name: hashlib.sha256((REPO / name).read_bytes()).hexdigest()
            for name in ["pipeline/extraction.py", "pipeline/table_router.py"]
        },
        "reproduction_scope": "Current deterministic markdown parser; no live model or complete extraction pipeline rerun",
        "current_parser_bug_reproduced": True,
        "results": results,
        "archived_count_provenance": json.loads(archived["count_provenance"]),
        "archived_facts": facts,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE / "reproduction.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in receipt.items()
                if key not in ["archived_facts"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
