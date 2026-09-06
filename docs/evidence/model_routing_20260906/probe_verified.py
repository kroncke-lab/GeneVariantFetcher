"""Bounded live probe after verifying actual SDK transmission offline."""

import os
import json
from pathlib import Path
from campaign import environment, HERE
from budget_guard import install

os.environ.update(environment("astra_medium_verified"))
os.environ["GVF_EXPERIMENT_ARM"] = "verified_smoke"
install()
from utils.llm_utils import litellm_completion
from utils.llm_trace import configure_llm_tracing, llm_trace_scope

configure_llm_tracing(HERE / "verified_smoke_traces", run_id="20260906_verified_smoke")
results = []
for model, effort in [
    ("azure_ai/gpt-6-astra", "medium"),
    ("azure_ai/grok-4.6", "high"),
]:
    with llm_trace_scope(stage="deployment_compatibility_smoke"):
        response = litellm_completion(
            model=model,
            messages=[{"role": "user", "content": 'Return exactly JSON {"ok":true}.'}],
            reasoning_effort=effort,
            max_tokens=32000,
            temperature=0,
            response_format={"type": "json_object"},
            timeout=120,
            num_retries=0,
        )
    assert json.loads(response.choices[0].message.content) == {"ok": True}
    result = {
        "requested_model": model,
        "transmitted_effort_verified_by_offline_sdk_transport_test": effort,
        "response_model": response.model,
        "usage": response.usage.model_dump(),
        "status": "passed",
    }
    results.append(result)
    (HERE / "verified_smoke.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(result), flush=True)
