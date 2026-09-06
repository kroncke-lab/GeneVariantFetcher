"""Bounded 2x2 diagnostic: reasoning effort and JSON mode on the new deployment."""

import concurrent.futures
import json
import os
import time
from campaign import environment, HERE
from budget_guard import install

os.environ.update(environment("grok46_verified"))
os.environ["GVF_EXPERIMENT_ARM"] = "grok_mode_probe"
install()
from utils.llm_utils import litellm_completion
from utils.llm_trace import configure_llm_tracing, llm_trace_scope

configure_llm_tracing(
    HERE / "grok_mode_probe_traces", run_id="20260906_grok_mode_probe"
)


def probe(pair):
    effort, json_mode = pair
    started = time.monotonic()
    kwargs = {"response_format": {"type": "json_object"}} if json_mode else {}
    try:
        with llm_trace_scope(
            stage="deployment_mode_probe", reasoning_effort=effort, json_mode=json_mode
        ):
            response = litellm_completion(
                model="azure_ai/grok-4.6",
                messages=[
                    {"role": "user", "content": 'Return exactly JSON {"ok":true}.'}
                ],
                max_tokens=1024,
                reasoning_effort=effort,
                temperature=0,
                timeout=65,
                num_retries=0,
                max_retries=0,
                **kwargs,
            )
        result = {
            "status": "returned",
            "usage": response.usage.model_dump(),
            "text": response.choices[0].message.content,
        }
    except Exception as exc:
        result = {
            "status": "failed",
            "error": type(exc).__name__ + ": " + str(exc)[:600],
        }
    result.update(
        effort=effort, json_mode=json_mode, seconds=time.monotonic() - started
    )
    (HERE / f"grok_mode_{effort}_{json_mode}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)


with concurrent.futures.ThreadPoolExecutor(4) as pool:
    list(pool.map(probe, [(e, j) for e in ("medium", "high") for j in (False, True)]))
