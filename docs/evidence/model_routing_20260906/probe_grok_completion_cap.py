"""Check the reasoning-aware output-cap spelling on Grok's Azure endpoint."""

import concurrent.futures
import json
import os
import time
from campaign import environment, HERE
from budget_guard import install

os.environ.update(environment("grok46_verified"))
os.environ["GVF_EXPERIMENT_ARM"] = "grok_completion_cap_probe"
install()
from utils.llm_utils import litellm_completion
from utils.llm_trace import configure_llm_tracing

configure_llm_tracing(
    HERE / "grok_completion_cap_traces", run_id="grok_completion_cap_probe"
)


def probe(mode):
    start = time.monotonic()
    kw = {"response_format": {"type": "json_object"}} if mode else {}
    try:
        r = litellm_completion(
            model="azure_ai/grok-4.6",
            messages=[{"role": "user", "content": 'Return exactly JSON {"ok":true}.'}],
            max_completion_tokens=1024,
            reasoning_effort="high",
            timeout=65,
            num_retries=0,
            max_retries=0,
            **kw,
        )
        data = {
            "status": "returned",
            "usage": r.usage.model_dump(),
            "text": r.choices[0].message.content,
        }
    except Exception as e:
        data = {"status": "failed", "error": str(e)[:500]}
    data.update(
        json_mode=mode,
        seconds=time.monotonic() - start,
        cap_parameter="max_completion_tokens",
    )
    (HERE / f"grok_completion_cap_{mode}.json").write_text(
        json.dumps(data, indent=2) + "\n"
    )
    print(json.dumps(data), flush=True)


with concurrent.futures.ThreadPoolExecutor(2) as pool:
    list(pool.map(probe, [False, True]))
