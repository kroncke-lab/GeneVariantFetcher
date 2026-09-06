"""Live production routing/serialization check with a retry-free SDK client."""

import json
import os
import sys
import time
from client import HERE, BASE, key, sha, write, update

sys.path.insert(0, str(HERE.parents[2]))
secret = key()
os.environ["AZURE_AI_API_BASE"] = BASE
os.environ["AZURE_AI_API_KEY"] = secret
os.environ["AZURE_AI_MODEL_ROUTES"] = ""
import httpx
from openai import OpenAI
from utils.llm_utils import litellm_completion

name = "grok_production_sdk_low"
reserve = 0.01
assert not (HERE / "responses" / (name + ".json")).exists()


def begin(d):
    assert sum(c["accounted_usd"] for c in d["calls"]) + reserve <= d["limit_usd"]
    assert not any(c["name"] == name for c in d["calls"])
    d["calls"].append(
        dict(
            name=name,
            model="grok-4.6",
            status="reserved",
            accounted_usd=reserve,
            reserved_usd=reserve,
            started_unix=time.time(),
        )
    )


update(begin)
sent = []


def capture(r):
    sent.append(json.loads(r.content))


start = time.monotonic()
result = dict(name=name, retries=0, reserved_usd=reserve, accounted_usd=reserve)
try:
    with httpx.Client(
        transport=httpx.HTTPTransport(retries=0),
        timeout=55,
        event_hooks={"request": [capture]},
    ) as http:
        with OpenAI(
            api_key=secret, base_url=BASE, max_retries=0, http_client=http
        ) as sdk:
            r = litellm_completion(
                model="azure_ai/grok-4.6",
                client=sdk,
                messages=[dict(role="user", content="Return exactly OK.")],
                max_tokens=256,
                reasoning_effort="low",
                temperature=0,
                response_format=None,
                num_retries=0,
                max_retries=0,
                timeout=55,
            )
    result.update(
        status="returned", response=r.model_dump(), usage=r.usage.model_dump()
    )
    u = result["usage"]
    ni = u["prompt_tokens"]
    no = max(u["completion_tokens"], u["total_tokens"] - ni)
    result.update(
        api_proxy_usd=(ni * 2 + no * 6) / 1e6, accounted_usd=(ni * 2.5 + no * 6) / 1e6
    )
except Exception as e:
    result.update(status="failed", error=str(e)[:300])
result.update(seconds=time.monotonic() - start, serialized_requests=sent)
write(HERE / "responses" / (name + ".json"), result)


def finish(d):
    c = next(c for c in d["calls"] if c["name"] == name)
    c.update({k: result[k] for k in ["status", "accounted_usd", "seconds"]})
    c["usage_known"] = "usage" in result
    c["response_sha256"] = sha(HERE / "responses" / (name + ".json"))
    if "api_proxy_usd" in result:
        c["api_proxy_usd"] = result["api_proxy_usd"]


update(finish)
print(json.dumps({k: v for k, v in result.items() if k != "response"}, indent=2))
