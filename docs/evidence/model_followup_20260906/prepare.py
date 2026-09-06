"""Freeze bounded, source-only packets and source-adjudicated references."""

import json
import re
import time
from pathlib import Path
from bs4 import BeautifulSoup
from client import HERE, write, sha

ROOT = HERE.parents[2]
FROZEN = ROOT / "validation_runs/model_routing_20260906/frozen_sources"
COLS = [
    "sources",
    "variants",
    "family",
    "person",
    "genotype",
    "phenotype",
    "index_case",
    "endpoint_status",
]
PROMPT = """Read only the supplied current-study source packet. Return a compact patient roster, not scalar totals. Use only individually identified people in the target table; do not manufacture people from aggregate counts, cells, healthy-control denominators, repeated tests or earlier publications. A person carrying multiple target-gene variants is one record with multiple variants. Index cases already occur in the roster. Include explicitly listed noncarriers with genotype N. Canonical genotype values Y/N/?; missing family or index status is JSON null. Preserve literal phenotype evidence (including No Dx or unknown marks). Strip index stars and footnote letters from person IDs, not digits. Preserve variant spelling except remove gene prefix and p. prefix and whitespace; use the target variant column specified below. Source references must include the actual table row and any necessary footnote/cohort/join line; never cite an unrelated row. Record identifiers only from supplied source. No external facts or medical inference. Map endpoint_status to positive/negative/unknown only according to the defined endpoint and source; asymptomatic does not generally mean disease negative, and follow-up status does not erase baseline disease. Preserve uncertainty, don't silently reconcile contradictions. If the packet has no individually identified patient roster, return rows:[] and explain why in limitations; aggregate counts are not person records. Keep explanations short. Return JSON only: {"columns":["sources","variants","family","person","genotype","phenotype","index_case","endpoint_status"],"rows":[...],"limitations":[]}. Each row is an 8-element array; sources and variants are arrays of strings. Do not repeat column names or provide scalar counts.\n"""
packets = []
refs = {}


def load(g, p):
    path = FROZEN / g / p / (p + "_FULL_CONTEXT.md")
    return path, path.read_text().splitlines()


def add(name, g, p, path, lines, selected, task, reference, cap=4096):
    evidence = "\n".join(f"{k}: {v}" for k, v in selected)
    packet = dict(
        name=name,
        gene=g,
        pmid=p,
        source_path=str(path.relative_to(ROOT)),
        source_sha256=sha(path),
        source_units=dict(selected),
        task=task,
        output_cap=cap,
        timeout=180,
        prompt=PROMPT
        + "\nTARGET: "
        + g
        + "; PMID "
        + p
        + "\nTASK: "
        + task
        + "\nSOURCE:\n"
        + evidence,
    )
    packets.append(packet)
    refs[name] = reference


def row(src, vs, fam, person, mut, ph, idx, status):
    return dict(zip(COLS, [src, vs, fam, person, mut, ph, idx, status]))


# Original DOC layout supplies independent cell reference, before model outputs.
g, p = "MYBPC3", "20433692"
path, lines = load(g, p)
html = (
    ROOT
    / "validation_runs/model_routing_20260906/source_layout/1471-2350-11-67-S1.html"
)
soup = BeautifulSoup(html.read_text(), "html.parser")
grid = []
spans = {}
origins = {}
inherited = []
for ri, tr in enumerate(soup.find("table").find_all("tr")):
    values = {c: t for c, (until, t) in spans.items() if until >= ri}
    inherited.append(sorted(values))
    c = 0
    for cell in tr.find_all(["td", "th"], recursive=False):
        while c in values:
            c += 1
        assert int(cell.get("colspan", 1)) == 1
        t = " ".join(cell.get_text(" ", strip=True).split())
        values[c] = t
        spans[c] = (ri + int(cell.get("rowspan", 1)) - 1, t)
        c += 1
    assert len(values) == 17
    grid.append([values[c] for c in range(17)])
mdrows = [
    i for i in range(287, 338) if re.search(r"\| (?:I{1,3}|IV)[:-]\d+\*? \|", lines[i])
]
assert len(mdrows) == 48 == len(grid) - 1
reference = []
for i, (cells, li) in enumerate(zip(grid[1:], mdrows), 1):
    v, f, person, m, _, ph = cells[:6]
    sources = [f"L{li + 1}"]
    f = f or "H49"
    if not cells[1]:
        sources += ["L120"]
    reference.append(
        row(
            sources,
            [v],
            f,
            person.replace("*", ""),
            "Y" if m == "Y" else "N",
            ph,
            "*" in person,
            "positive"
            if ph.isdigit()
            else "negative"
            if ph.replace(" ", "").lower() == "nodx"
            else "unknown",
        )
    )
task = "Roster is Additional file 1 only. Copy Age diagnosis as phenotype. Endpoint is table-defined diagnosis at recorded diagnosis/last observation: numeric diagnosis age is positive, No Dx is table-defined negative, unknown mark remains unknown. This is not an ECG-normal endpoint. Table footnote defines No Dx. Prose is available to support blank-family joins and expose contradictions; do not add people mentioned only in prose. Shifted Markdown continuation cells require care; variant/family can inherit only supported grouping. Two original family cells for IV-1 and IV-5 are blank rather than merged; use prose only if it supplies the family."
selected = [("L120", lines[119])] + [(f"L{i + 1}", lines[i]) for i in range(284, 339)]
add("myb204_markdown", g, p, path, lines, selected, task, reference)
structured = [("L120", lines[119]), ("H0", " | ".join(grid[0]))]
for i, cells in enumerate(grid[1:], 1):
    structured.append(
        (
            f"H{i}",
            " | ".join(cells)
            + f" [columns inherited from actual rowspan: {inherited[i]}]",
        )
    )
structured.append(("L339", lines[338]))
ref2 = json.loads(json.dumps(reference))
for i, r in enumerate(ref2, 1):
    r["sources"][0] = f"H{i}"
add(
    "myb204_structured",
    g,
    p,
    path,
    lines,
    structured,
    task
    + " Source H rows are an original-DOC HTML grid with only real rowspans expanded. Empty cells remain empty.",
    ref2,
)
packets[-1]["original_html_sha256"] = sha(html)
# Multi-variant roster: scope explicitly excludes relatives in family-history cells.
g, p = "RYR2", "30403697"
path, lines = load(g, p)
reference = []
for li in range(106, 121):
    cells = [x.strip() for x in lines[li].strip("|").split("|")]
    person = cells[0]
    vs = []
    for cell in cells[5:7]:
        if "RYR2" in cell:
            vs.append(re.sub(r"\s+", "", cell.split("p.")[1]))
    src = [f"L{li + 1}"]
    if person == "14":
        vs += ["c.3599-9delT", "c.14091-11dupT"]
        src += ["L124"]
    status = "negative" if person in ["6", "12"] else "positive"
    reference.append(
        row(src, vs, None, person, "Y", cells[10], cells[4] == "Yes", status)
    )
assert len(reference) == 15
add(
    "ryr304_multivariant",
    g,
    p,
    path,
    lines,
    [(f"L{i + 1}", lines[i]) for i in range(100, 127)],
    "Roster is numbered subjects in Table 1 only; do not add relatives from Family History. Copy Symptoms as phenotype. Include extra target-gene variants from the subject-14 footnote. Endpoint is documented VT and/or sudden cardiac arrest at presentation, NOT all CPVT diagnoses or any symptoms. The cohort prose specifies the exceptions. Family is null because no family IDs are supplied. Keep subject 5 self-sibling wording as a limitation; do not merge people.",
    reference,
)
# Flat original letter, chunked at literal person-ID boundaries.
g, p = "RYR2", "25435091"
path, lines = load(g, p)
text = path.read_text()
start = text.index("Table Clinical Characteristics")
end = text.index("Scientific letter", start)
table = text[start:end]
ids = ["II:1", "II:3", "II:6", "II:8", "II:9", "II:15", "III:4", "III:9"]
positions = [table.index(i + " ") for i in ids]
foot = table.index("Big, bigeminy;")
units = [("T0", table[: positions[0]])]
phs = [
    "VE, Big, D, NSVT",
    "VE, Big, NSVT",
    "VE, Big, NSVT",
    "VE, Big, D",
    "VE",
    "VE, Big, Trig",
    "VE, Big, D",
    "VE, Big, D",
]
reference = []
for i, person in enumerate(ids):
    t = table[positions[i] : positions[i + 1] if i + 1 < len(ids) else foot].strip()
    units.append((f"T{i + 1}", t))
    reference.append(
        row(
            [f"T{i + 1}", "T9", "N1"],
            ["C2277R"],
            None,
            person,
            "Y",
            phs[i],
            i == 0,
            "negative" if person == "II:9" else "positive",
        )
    )
units.append(("T9", table[foot:]))
ns = text.index("Our objective is to describe a kindred")
ne = text.index("FUNDING", ns)
units.append(("N1", text[ns:ne]))
add(
    "ryr254_narrative",
    g,
    p,
    path,
    lines,
    units,
    "Roster is the eight individually identified members in the clinical table. Copy Maximum arrhythmia in the initial EST as phenotype, omit superscript footnote letters. Endpoint is initial diagnostic CPVT on EST/adrenaline, not previous symptoms or later treatment response. Family null (no identifier). Proband true only when indicated; others false. Use table footnote and narrative for diagnosis threshold.",
    reference,
)
# Enumerated patient IDs; includes related people, not an unrelated-proband total.
g, p = "MYBPC3", "21302287"
path, lines = load(g, p)
reference = []
for li in range(119, 137):
    c = [s.strip() for s in lines[li].strip("|").split("|")]
    assert len(c) == 6
    for person in c[5].split(","):
        person = re.sub("[a-z]$", "", person.strip())
        reference.append(
            row(
                [f"L{li + 1}", "L13"],
                [re.sub(r"\s+", "", c[3])],
                None,
                person,
                "Y",
                "HCM",
                None,
                "positive",
            )
        )
assert len(reference) == 22
add(
    "myb213_patient_ids",
    g,
    p,
    path,
    lines,
    [(f"L{i + 1}", lines[i]) for i in [12, 84] + list(range(112, 138)) + [166, 174]],
    "Roster is Table 3 Panel (a) only, target variant is the Protein mutation column (copy reported spelling, even if notation appears inconsistent). Every listed patient ID is a person; suffix c/d is a footnote marker, not another patient. Include all listed relatives, not just unrelated probands. Copy HCM as phenotype from the enrollment sentence. Endpoint is enrolled clinical HCM at ascertainment. Do not guess family IDs or proband status. Do not include Panel (b).",
    reference,
)
# Negative roster cases: source contains aggregate evidence, but no identified persons.
for name, g, p, indices, task in [
    (
        "ryr189_missing_table",
        "RYR2",
        "18929323",
        list(range(18, 31)),
        "This packet includes study-subject prose but not Table 1 itself. Do not turn cohort totals or repeated Holters into individual people. No roster can be reconstructed without individual identifiers.",
    ),
    (
        "ryr193_abstract",
        "RYR2",
        "19398417",
        None,
        "Abstract only. Copy individually identified people only; aggregate family summaries do not identify individual patient rows.",
    ),
    (
        "ryr258_functional",
        "RYR2",
        "25814417",
        list(range(34, 42)),
        "Only a functional-assay section is supplied. Cells, alleles, channels, experimental replicates and other nonhuman units are not patients. This is a source-section negative control, not a claim that the whole paper lacks patient counts.",
    ),
]:
    path, lines = load(g, p)
    selected = (
        [("A1", json.loads(path.read_text())["abstract"])]
        if indices is None
        else [(f"L{i + 1}", lines[i]) for i in indices]
    )
    add(name, g, p, path, lines, selected, task, [], cap=1024)
assert not (HERE / "prepared.json").exists()
(HERE / "packets").mkdir(exist_ok=True)
for packet in packets:
    write(HERE / "packets" / (packet["name"] + ".json"), packet)
write(
    HERE / "source_reference.json",
    dict(
        classification="Investigator-prepared source adjudication, previously opened papers, no gold access in preparation; not independent human validation",
        columns=COLS,
        rows=refs,
    ),
)
write(
    HERE / "prepared.json",
    dict(
        prepared_unix=time.time(),
        packets={
            p["name"]: sha(HERE / "packets" / (p["name"] + ".json")) for p in packets
        },
        reference_sha256=sha(HERE / "source_reference.json"),
        prompt_sha256=__import__("hashlib").sha256(PROMPT.encode()).hexdigest(),
        source_preparation_sha256=sha(Path(__file__)),
        head=__import__("subprocess")
        .check_output(["git", "rev-parse", "HEAD"], text=True)
        .strip(),
        arms=["astra_low", "astra_medium"],
        order_seed=20260906,
        notes="Eight source packets from seven opened papers. Structured-vs-Markdown pair is separate representation ablation. One draw per cell, no population inference. All new outputs locked before grading.",
    ),
)
print(
    json.dumps(
        {
            p["name"]: {
                "people": len(refs[p["name"]]),
                "prompt_bytes": len(p["prompt"].encode()),
                "cap": p["output_cap"],
            }
            for p in packets
        },
        indent=2,
    )
)
