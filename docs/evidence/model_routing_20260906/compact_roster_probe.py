"""Single source-only transcription diagnostic; never read gold or merge counts."""

import hashlib
import json
import os
import time
from pathlib import Path
from campaign import environment, HERE, FROZEN
from budget_guard import install


def main():
    target = HERE / "compact_roster_result.json"
    assert not target.exists(), "Preserve previous diagnostic"
    source = FROZEN / "MYBPC3/20433692/20433692_FULL_CONTEXT.md"
    lines = source.read_text().splitlines()
    start = next(
        i
        for i, s in enumerate(lines)
        if s.startswith("Additional file 1- Clinical characteristics")
    )
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("ABPR:"))
    evidence = "\n".join(f"L{i + 1}: {lines[i]}" for i in range(start, end + 1))
    prompt = (
        """Transcribe this current-study MYBPC3 clinical supplement into compact person records. This is a transcription task, not scalar count extraction or clinical diagnosis. Only this article source is evidence. No external facts.
The DOC-to-Markdown conversion omitted merged variant/family cells on continuation rows and shifted columns left. Preserve evidence line IDs. Inherit variant/family only when table grouping is unambiguous; mark uncertain rows rather than guess. Exclude wrapped clinical-cell continuation fragments that are not people, but list those line IDs in ignored_fragments. Keep carriers (Y) and non-carriers (N/no) separate. Do not convert No Dx into healthy: copy the diagnosis cell literally. Preserve unknown/question-mark cells. Index cases have a star; they are already people in the table, not additional people to add later. Do not invent people or split a homozygote into two people. Return one record per table person, deduplicated by family plus person ID. Do not emit scalar totals, long explanations, or repeated table quotes; source line references bind the evidence.
Return JSON: {"columns":["line","variant","family","person","mut","diagnosis_cell","index_case","inherited_group"],"rows":[[290,"D75N","H73","II:3","N","NoDx",false,false]],"ambiguous_rows":[],"ignored_fragments":[]}.
The example is a format illustration; transcribe actual source line IDs. Output only JSON.

SOURCE LINES:
"""
        + evidence
    )
    (HERE / "compact_roster_prompt.txt").write_text(prompt + "\n")
    os.environ.update(environment("astra_medium_verified"))
    os.environ["GVF_EXPERIMENT_ARM"] = "compact_roster_low_diagnostic"
    install()
    from utils.llm_utils import litellm_completion
    from utils.llm_trace import configure_llm_tracing, llm_trace_scope

    configure_llm_tracing(
        HERE / "compact_roster_traces", run_id="compact_roster_low_diagnostic"
    )
    started = time.monotonic()
    result = {
        "classification": "unvalidated source-roster transcription diagnostic, no scalar score",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_lines": [start + 1, end + 1],
        "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
    }
    try:
        with llm_trace_scope(
            gene="MYBPC3", pmid="20433692", stage="compact_person_roster"
        ):
            r = litellm_completion(
                model="azure_ai/gpt-6-astra",
                messages=[{"role": "user", "content": prompt}],
                reasoning_effort="low",
                max_tokens=8192,
                response_format={"type": "json_object"},
                timeout=240,
                max_retries=0,
                num_retries=0,
            )
        result.update(
            status="returned",
            response_model=r.model,
            usage=r.usage.model_dump(),
            finish_reason=r.choices[0].finish_reason,
            raw_text=r.choices[0].message.content,
        )
        result["parsed"] = json.loads(result["raw_text"])
    except Exception as exc:
        result.update(status="failed", error=str(exc)[:1000])
    result["seconds"] = time.monotonic() - started
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in {"raw_text", "parsed"}}
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
