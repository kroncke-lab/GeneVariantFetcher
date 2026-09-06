"""Compare minimal direct HTTP with SDK requests that omit the stream field."""

import concurrent.futures
import json
import os
import time
import httpx
from campaign import environment, HERE
from budget_guard import install

os.environ.update(environment("grok46_verified"))
os.environ["GVF_EXPERIMENT_ARM"] = "grok_stream_field_probe"
install()
from utils.llm_utils import litellm_completion, azure_connection_for_model
from utils.llm_trace import configure_llm_tracing, capture_llm_call

configure_llm_tracing(
    HERE / "grok_stream_field_traces", run_id="grok_stream_field_probe"
)


def probe(kind):
    started = time.monotonic()
    try:
        if kind == "direct":
            base, key = azure_connection_for_model("azure_ai/grok-4.6")
            body = {
                "model": "grok-4.6",
                "messages": [
                    {"role": "user", "content": 'Return exactly JSON {"ok":true}.'}
                ],
                "max_completion_tokens": 1024,
                "reasoning_effort": "high",
            }

            def call():
                with httpx.Client(timeout=65) as client:
                    r = client.post(
                        base + "/chat/completions",
                        headers={"Authorization": "Bearer " + key},
                        json=body,
                    )
                r.raise_for_status()
                return r.json()

            response, _ = capture_llm_call(
                provider="azure_http_probe",
                requested_model="azure_ai/grok-4.6",
                resolved_model="grok-4.6",
                request=body,
                call=call,
            )
            result = {
                "status": "returned",
                "usage": response.get("usage"),
                "text": response["choices"][0]["message"]["content"],
            }
        else:
            kw = (
                {"response_format": {"type": "json_object"}, "max_tokens": 1024}
                if kind == "sdk_json"
                else {"max_completion_tokens": 1024}
            )
            r = litellm_completion(
                model="azure_ai/grok-4.6",
                messages=[
                    {"role": "user", "content": 'Return exactly JSON {"ok":true}.'}
                ],
                reasoning_effort="high",
                stream=None,
                timeout=65,
                num_retries=0,
                max_retries=0,
                **kw,
            )
            result = {
                "status": "returned",
                "usage": r.usage.model_dump(),
                "text": r.choices[0].message.content,
            }
    except Exception as exc:
        result = {"status": "failed", "error": str(exc)[:500]}
    result.update(kind=kind, seconds=time.monotonic() - started)
    (HERE / f"grok_stream_field_{kind}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)


with concurrent.futures.ThreadPoolExecutor(3) as pool:
    list(pool.map(probe, ["direct", "sdk_plain", "sdk_json"]))
