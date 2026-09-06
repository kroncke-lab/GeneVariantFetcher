import hashlib
import json
from pathlib import Path
from bs4 import BeautifulSoup

r = Path(__file__).resolve().parents[3]
e = r / "docs/evidence/model_routing_20260906"
s = r / "validation_runs/model_routing_20260906/source_layout"
doc = (
    r
    / "corpus/MYBPC3/20433692/20433692_supplements/PMC2880974_supplements/1471-2350-11-67-S1.DOC"
)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


h = s / "1471-2350-11-67-S1.html"
soup = BeautifulSoup(h.read_text(), "html.parser")
table = soup.find("table")
grid = []
spans = {}
for ri, tr in enumerate(table.find_all("tr")):
    row = {c: text for c, (until, text) in spans.items() if until >= ri}
    col = 0
    for cell in tr.find_all(["td", "th"], recursive=False):
        while col in row:
            col += 1
        assert int(cell.get("colspan", 1)) == 1
        txt = " ".join(cell.get_text(" ", strip=True).split())
        row[col] = txt
        spans[col] = (ri + int(cell.get("rowspan", 1)) - 1, txt)
        col += 1
    assert len(row) == 17, (ri, row)
    grid.append([row[c] for c in range(17)])
result = json.loads((e / "compact_roster_result.json").read_text())
records = [dict(zip(result["parsed"]["columns"], x)) for x in result["parsed"]["rows"]]
assert len(records) == len(grid) - 1 == 48
errors = []
family_blank = []
last_family = None
for i, (cells, record) in enumerate(zip(grid[1:], records), 1):
    variant, family, person, mut, sex, diagnosis = cells[:6]
    if not family:
        family_blank.append(
            {
                "original_person_row": i,
                "person": person,
                "carried_family": last_family,
                "note": "Family cell is explicitly blank, not a merged H49 cell; prose line 120 supplies H49 for these noncarriers.",
            }
        )
        family = last_family
    else:
        last_family = family
    expected = dict(
        variant=variant,
        family=family,
        person=person.replace("*", ""),
        mut=mut,
        diagnosis_cell=diagnosis,
        index_case="*" in person,
    )
    for field, value in expected.items():
        if record[field] != value:
            errors.append(
                {"row": i, "field": field, "original": value, "model": record[field]}
            )
body = " ".join(soup.get_text(" ", strip=True).split())
assert "No Dx: Unaffected or healthy" in body
out = {
    "classification": "Source-only original-DOC layout cross-check, not gold/count accuracy",
    "original_doc": str(doc.relative_to(r)),
    "original_doc_sha256": sha(doc),
    "html_sha256": sha(h),
    "rendered_pdf_sha256": sha(s / "1471-2350-11-67-S1.pdf"),
    "conversion": "LibreOffice read-only HTML and PDF exports; original corpus file unchanged",
    "visual_pages_inspected": [1, 2, 3],
    "person_rows": 48,
    "mechanical_comparison_fields": [
        "variant",
        "family",
        "person",
        "mut",
        "diagnosis_cell",
        "index_case",
    ],
    "differences": errors,
    "blank_family_cells_carried_with_prose_support": family_blank,
    "footnote": "No Dx: Unaffected or healthy",
    "endpoint_qualification": "The table includes 9 No Dx carriers; prose distinguishes 6 healthy carriers from 4 with suggestive ECG but no diagnostic HCM (3 No Dx plus 1 unknown diagnosis). Do not equate all table No Dx with ECG-normal or asymptomatic. This source audit does not choose the benchmark endpoint or modify predictions.",
}
(e / "original_layout_audit.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
