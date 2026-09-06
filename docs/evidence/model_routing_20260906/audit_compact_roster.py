"""Mechanical transcription audit for this source diagnostic; no gold counts."""

import hashlib
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    result = json.loads((HERE / "compact_roster_result.json").read_text())
    source = (
        ROOT
        / "validation_runs/model_routing_20260906/frozen_sources/MYBPC3/20433692/20433692_FULL_CONTEXT.md"
    )
    assert hashlib.sha256(source.read_bytes()).hexdigest() == result["source_sha256"]
    lines = source.read_text().splitlines()
    payload = result["parsed"]
    assert payload["columns"] == [
        "line",
        "variant",
        "family",
        "person",
        "mut",
        "diagnosis_cell",
        "index_case",
        "inherited_group",
    ]
    expected = {}
    variant = family = None
    fragments = []
    for number in range(result["source_lines"][0], result["source_lines"][1] + 1):
        line = lines[number - 1]
        if (
            not line.startswith("|")
            or line.startswith("|---")
            or "Mutation | Fam." in line
        ):
            continue
        cells = [s.strip() for s in line.strip("|").split("|")]
        positions = [
            i for i, v in enumerate(cells[:3]) if re.fullmatch(r"[IVX]+[:-]\d+\*?", v)
        ]
        if not positions:
            fragments.append(number)
            continue
        assert len(positions) == 1
        p = positions[0]
        if p == 2:
            variant, family = cells[:2]
        elif p == 1:
            assert re.fullmatch(r"H\d+", cells[0])
            family = cells[0]
        assert variant and family
        expected[number] = [
            number,
            variant,
            family,
            cells[p].rstrip("*"),
            cells[p + 1],
            cells[p + 3],
            cells[p].endswith("*"),
            p < 2,
        ]
    observed = {r[0]: r for r in payload["rows"]}
    assert len(observed) == len(payload["rows"]), "duplicate output line"
    differences = [
        {
            "line": n,
            "expected_transcription": expected.get(n),
            "observed": observed.get(n),
        }
        for n in sorted(set(expected) | set(observed))
        if expected.get(n) != observed.get(n)
    ]
    identities = [(r[2], r[3]) for r in payload["rows"]]
    carriers = [r for r in payload["rows"] if r[4] == "Y"]
    output = {
        "classification": "source-text transcription consistency only; no reference score or clinical acceptance",
        "source_sha256": result["source_sha256"],
        "source_person_rows": len(expected),
        "returned_rows": len(observed),
        "differences": differences,
        "unique_family_person_ids": len(set(identities)),
        "ignored_fragments_match": fragments == payload["ignored_fragments"],
        "mutation_status": dict(Counter(r[4] for r in payload["rows"])),
        "index_people": sum(r[6] for r in payload["rows"]),
        "carrier_diagnosis_cell_groups": {
            "numeric_diagnosis_age": sum(r[5].isdigit() for r in carriers),
            "No_Dx": sum(r[5].replace(" ", "") == "NoDx" for r in carriers),
            "unknown": sum(r[5] == "¿" for r in carriers),
        },
        "per_variant_carrier_row_counts": dict(Counter(r[1] for r in carriers)),
        "qualification": "Audit uses the same text conversion and explicit continuation-row carry-forward convention. It does not validate the original DOC geometry or resolve phenotype definitions, contradictory family IDs in prose, or gold agreement. No Dx is not promoted to asymptomatic/healthy, and numeric diagnosis age is not a new accepted affected count. These totals never enter benchmark predictions.",
    }
    (HERE / "compact_roster_audit.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
