"""Final bounded availability diagnostic, with source-free inputs and no retries."""

import concurrent.futures
import json
import os
import time
import httpx
from campaign import environment, HERE
from budget_guard import install

os.environ.update(environment("grok46_verified"))
os.environ["GVF_EXPERIMENT_ARM"] = "grok_final_health"
install()
from utils.llm_utils import azure_connection_for_model
from utils.llm_trace import configure_llm_tracing, capture_llm_call

configure_llm_tracing(HERE / "grok_final_health_traces", run_id="grok_final_health")


def probe(kind):
    base, key = azure_connection_for_model("azure_ai/grok-4.6")
    body = {"model": "grok-4.6"}
    if kind == "responses_low":
        body.update(
            input="Return exactly OK.",
            max_output_tokens=32000,
            reasoning={"effort": "low"},
        )
        path = "/responses"
    else:
        body.update(
            messages=[{"role": "user", "content": "Return exactly OK."}],
            max_tokens=32000,
        )
        path = "/chat/completions"
    started = time.monotonic()
    try:

        def call():
            with httpx.Client(timeout=55) as client:
                r = client.post(
                    base + path, headers={"Authorization": "Bearer " + key}, json=body
                )
            r.raise_for_status()
            return r.json()

        response, ref = capture_llm_call(
            provider="azure_http_probe",
            requested_model="azure_ai/grok-4.6",
            resolved_model="grok-4.6",
            request=body,
            call=call,
        )
        result = {"status": "returned", "response": response, "trace_id": ref}
    except Exception as exc:
        result = {"status": "failed", "error": str(exc)[:500]}
    result.update(kind=kind, seconds=time.monotonic() - started)
    (HERE / f"grok_final_health_{kind}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        list(pool.map(probe, ["responses_low", "chat_default"]))
