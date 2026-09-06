"""Freeze a paired source-reading decision panel. No API calls or gold access."""

import json
import random
import re
import time
from pathlib import Path

from client import HERE, sha, write

ROOT = HERE.parents[2]
FROZEN = ROOT / "validation_runs/model_routing_20260906/frozen_sources"
PROMPT = """Extract only the requested human per-variant counts from this source packet.
The queries define the study population and clinical endpoint, not the answers.
Return one claim per query ID. Use null when the source cannot support an exact
count for that variant, population and endpoint. Unknown is not zero. A list of
variants does not establish one carrier per row. Do not turn cohort enrollment,
population allele frequencies, prior-literature cases, relatives without stated
genotype, repeat tests or percentages rounded in prose into exact patient counts.
Include explicitly documented carrier relatives when the requested population
includes them; respect variant-specific inheritance and overlapping people.
Use the paper's actual phenotype labels and timepoints. Disease-negative,
symptom-free, no diagnostic test result, and unknown are different states.
If the query requests a specific endpoint, do not substitute a different one.
An explicit count is a literal published number bound to this variant and role.
A derived count is a transparent aggregation or join of source-supported people
or disjoint groups. For derived counts state the arithmetic or membership in
reason. Do not label arithmetic as explicit. sources must use exact provided
IDs, including all needed headers, footnotes and cohort/join evidence. quote is
a short verbatim contiguous excerpt (at most 600 characters), or null when no
count is supplied; additional support can be cited in sources. Do not invent
source IDs, clean up a quotation, or consult external information. Keep reasons
concise. Return only JSON matching the schema; report unresolved limitations.
"""


def source(gene, pmid, ranges):
    path = FROZEN / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    lines = path.read_text().splitlines()
    units = {f"L{i}": lines[i - 1] for a, b in ranges for i in range(a, b + 1)}
    return path, units


def queries(variants, fields=("carriers", "affected", "unaffected")):
    return [
        dict(id=f"q{i + 1}_{f}", variant=v, field=f)
        for i, v in enumerate(variants)
        for f in fields
    ]


def main():
    assert not (HERE / "plan.json").exists(), "Frozen plan already exists"
    packets = []
    references = {}

    def add(
        name,
        gene,
        pmid,
        ranges,
        variants,
        values,
        scope,
        endpoint,
        fields=("carriers", "affected", "unaffected"),
        custom=None,
    ):
        path, units = source(gene, pmid, ranges)
        if custom is not None:
            units = custom(path.read_text())
        qs = queries(variants, fields)
        assert len(qs) == len(values)
        packet = dict(
            name=name,
            gene=gene,
            pmid=pmid,
            source_path=str(path.relative_to(ROOT)),
            source_sha256=sha(path),
            source_units=units,
            scope=scope,
            endpoint=endpoint,
            queries=qs,
        )
        packet["prompt"] = (
            PROMPT
            + f"\nGENE: {gene}; PMID: {pmid}\n"
            + f"POPULATION: {scope}\nENDPOINT: {endpoint}\n"
            + "QUERIES: "
            + json.dumps(qs)
            + "\nSOURCE:\n"
            + "\n".join(f"{k}: {v}" for k, v in units.items())
        )
        write(HERE / "packets" / f"{name}.json", packet)
        references[name] = {q["id"]: value for q, value in zip(qs, values)}
        packets.append(packet)

    add(
        "scn200_table",
        "SCN5A",
        "20031634",
        [(109, 175), (184, 216)],
        ["p.Gly1408Arg", "c.3963+2T>C", "p.Ala665GlyfsX16"],
        [14, 4, 9, 10, 2, 8, 9, 6, 2],
        "The current study families, using the genotype-positive people in Table 1.",
        "BrS ECG positive versus explicitly Not BrS ECG positive, as defined by this paper. Preserve undetermined phenotype separately.",
    )
    add(
        "ryr304_relatives",
        "RYR2",
        "30403697",
        [(93, 124)],
        ["p.R417L", "p.R2028H", "p.Y4721C", "p.G4772S"],
        [3, 2, 2, 2],
        "All distinct people explicitly documented as carrying each variant in this current-study packet, including stated carrier relatives outside the numbered patient rows. Count documented people, not an estimate of unreported family size.",
        "Carrier counts only; no affected/unaffected query because family-history endpoints are not uniform.",
        fields=("carriers",),
    )

    def letter(text):
        text = text[text.index("Figure A: Family tree") :]
        text = text[: text.index("FUNDING")]
        # Preserve exact text and all scientific table/prose content. Split on
        # sentence boundaries; never manufacture a new table or count label.
        spans = re.split(r"(?<=[.!?]) (?=[A-Z0-9])", text)
        return {f"S{i + 1}": span for i, span in enumerate(spans)}

    add(
        "ryr254_timepoint",
        "RYR2",
        "25435091",
        [],
        ["p.C2277R"],
        [8, 7, 1],
        "Current genotyped family members, at initial clinical assessment.",
        "CPVT diagnostic phenotype on the initial exercise/adrenaline testing; do not substitute previous symptoms or response at follow-up.",
        custom=letter,
    )
    add(
        "ryr189_aggregate",
        "RYR2",
        "18929323",
        [(19, 29)],
        ["p.P2328S", "p.V4653F"],
        [13, None, None, 6, None, None],
        "The patients enrolled in the present Holter study, counted once per person.",
        "Clinical manifestation status at enrollment (symptomatic/phenotype positive versus explicitly asymptomatic/phenotype negative). Require variant-specific attribution of the phenotype split.",
    )
    add(
        "ryr193_abstract",
        "RYR2",
        "19398417",
        [],
        ["p.W4645R"],
        [4, 2, 2],
        "The family clinically and genetically evaluated in this report.",
        "Reported clinical symptom status among carriers; this query is symptoms, not an assertion of diagnostic test negativity.",
        custom=lambda text: {
            f"A{i + 1}": v
            for i, v in enumerate(json.loads(text)["abstract"].splitlines())
            if v
        },
    )
    add(
        "scn325_prior_counts",
        "SCN5A",
        "32533946",
        [(222, 255), (345, 350), (1099, 1117)],
        ["p.Thr220Ile", "p.Gly752Arg", "p.Glu1784Lys"],
        [None] * 9,
        "Human clinical participants newly studied in this publication; count no prior-literature or population-database individuals.",
        "Current-study Brugada phenotype-positive/negative carriers; in-vitro measurements are not patient phenotypes.",
    )
    add(
        "ryr258_roster",
        "RYR2",
        "25814417",
        [(27, 33), (49, 61), (225, 408)],
        ["p.G357S"],
        [179, 45, 133],
        "Living mutation-positive subjects at baseline in Supplementary Table 3. Keep other deceased and later-follow-up cohorts separate.",
        "Previous symptoms at baseline, as recorded in the Previous symptoms column: reported symptoms versus explicit No. Preserve missing values, and do not substitute VA or CVA results.",
    )
    add(
        "scn251_variant_list",
        "SCN5A",
        "25163546",
        [(1180, 1184), (1374, 1397)],
        ["c.704-2A>G", "p.T220N"],
        [None] * 6,
        "Current-study human participants carrying the requested variant.",
        "Per-variant DCM-affected versus explicitly unaffected people; the source must establish a patient count, not just an entry in a variant list.",
    )
    write(HERE / "reference_values.json", references)
    plans = []
    rng = random.Random(2026090607)
    for p in packets:
        props = dict(
            id={"type": "string", "enum": [q["id"] for q in p["queries"]]},
            value={"type": ["integer", "null"]},
            basis={"type": "string", "enum": ["explicit", "derived", "unknown"]},
            sources={
                "type": "array",
                "items": {"type": "string", "enum": list(p["source_units"])},
            },
            quote={"type": ["string", "null"]},
            reason={"type": "string"},
        )
        schema = dict(
            type="object",
            additionalProperties=False,
            required=["claims", "limitations"],
            properties=dict(
                claims=dict(
                    type="array",
                    items=dict(
                        type="object",
                        additionalProperties=False,
                        required=list(props),
                        properties=props,
                    ),
                ),
                limitations=dict(type="array", items={"type": "string"}),
            ),
        )
        arms = ["grok_first", "grok_repeat", "astra_low"]
        rng.shuffle(arms)
        for arm in arms:
            body = dict(
                model="gpt-6-astra" if arm == "astra_low" else "grok-4.6",
                messages=[dict(role="user", content=p["prompt"])],
                reasoning_effort="low",
                max_completion_tokens=4096,
                response_format=dict(
                    type="json_schema",
                    json_schema=dict(name="source_counts", strict=True, schema=schema),
                ),
            )
            plans.append(
                dict(name=p["name"] + "__" + arm, packet=p["name"], arm=arm, body=body)
            )
    write(
        HERE / "plan.json",
        dict(
            prepared_unix=time.time(),
            seed=2026090607,
            classification="Selected opened-source component diagnostic; not a held-out, full-paper or production benchmark.",
            budget_limit_usd=2.17,
            timeout_seconds=180,
            retries=0,
            decision="Complete the fixed panel before inspecting outcomes. Routine Astra inclusion requires incremental correct current-validator accepted fields, no new accepted wrong fields and benefit beyond a same-cost-cheaper Grok repeat. Zero unique benefit supports no routine inclusion now. A selective rescue pilot requires source-valid unique wins on at least two papers, replicated on the same generic contract; it still cannot bypass the derived-count acceptance gate. This is a pragmatic decision rule, not a population significance test.",
            stopping="Run every initial cell if reserves allow; no gold-driven early stopping. Report all failures and undispatched cells. After lock and score, use remaining allowance only for same-prompt replication of Astra-only wins. Never reveal reference answers to an API reader.",
            packet_sha256={
                p["name"]: sha(HERE / "packets" / (p["name"] + ".json"))
                for p in packets
            },
            reference_sha256=sha(HERE / "reference_values.json"),
            plans=plans,
        ),
    )
    print(
        json.dumps(
            dict(
                papers=len(packets),
                queries=sum(len(p["queries"]) for p in packets),
                calls=len(plans),
            )
        )
    )


if __name__ == "__main__":
    main()
