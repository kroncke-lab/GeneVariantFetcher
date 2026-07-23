import json

import pytest
import requests

from utils.responses_api import (
    REASONING_EFFORTS,
    ResponsesAPIClient,
    ResponsesLoopLimits,
)


class FakeResponse:
    def __init__(self, status_code=200, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload
        self.text = text

    def json(self):
        if isinstance(self._payload, Exception):
            raise self._payload
        return self._payload


class FakeTransport:
    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.calls = []

    def __call__(self, url, **kwargs):
        self.calls.append((url, kwargs))
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


class StepClock:
    def __init__(self, step=0.25):
        self.value = 0.0
        self.step = step

    def __call__(self):
        value = self.value
        self.value += self.step
        return value


def _message_response(response_id, value, *, usage=None):
    return FakeResponse(
        payload={
            "id": response_id,
            "status": "completed",
            "output": [
                {
                    "type": "message",
                    "content": [{"type": "output_text", "text": json.dumps(value)}],
                }
            ],
            "usage": usage or {},
        }
    )


def _tool_response(response_id, calls, *, usage=None):
    return FakeResponse(
        payload={
            "id": response_id,
            "status": "completed",
            "output": [
                {
                    "type": "function_call",
                    "call_id": call_id,
                    "name": name,
                    "arguments": arguments,
                }
                for call_id, name, arguments in calls
            ],
            "usage": usage or {},
        }
    )


def _client(transport, **kwargs):
    kwargs.setdefault("backoff_s", 0)
    return ResponsesAPIClient.azure(
        "https://example.services.ai.azure.com/openai/v1",
        "secret",
        transport=transport,
        **kwargs,
    )


def _run(client, **kwargs):
    arguments = {
        "model": "gpt-5.6-sol",
        "reasoning_effort": "max",
        "instructions": "Return grounded JSON.",
        "initial_input": "paper overview",
        "tools": [
            {
                "type": "function",
                "name": "read_table",
                "description": "Read one table",
                "parameters": {"type": "object", "properties": {}},
            }
        ],
        "tool_handlers": {},
        "text": {
            "format": {
                "type": "json_schema",
                "name": "result",
                "strict": True,
                "schema": {"type": "object"},
            }
        },
        "max_output_tokens": 4096,
    }
    arguments.update(kwargs)
    return client.run_tool_loop(**arguments)


def test_executes_all_tool_calls_and_repeats_continuation_contract():
    transport = FakeTransport(
        [
            _tool_response(
                "resp_1",
                [
                    ("call_a", "read_table", '{"table_id":"T1"}'),
                    ("call_b", "search", '{"query":"carrier"}'),
                ],
                usage={
                    "input_tokens": 100,
                    "input_tokens_details": {
                        "cached_tokens": 20,
                        "cache_write_tokens": 7,
                    },
                    "output_tokens": 40,
                    "output_tokens_details": {"reasoning_tokens": 30},
                    "total_tokens": 140,
                },
            ),
            _message_response(
                "resp_2",
                {"variants": []},
                usage={
                    "input_tokens": 50,
                    "cached_input_tokens": 10,
                    "cache_creation_input_tokens": 3,
                    "output_tokens": 12,
                    "reasoning_tokens": 8,
                    "total_tokens": 62,
                },
            ),
        ]
    )
    calls = []
    client = _client(transport, clock=StepClock())

    result = _run(
        client,
        tool_handlers={
            "read_table": lambda args: (
                calls.append(("read_table", args)) or {"rows": [1]}
            ),
            "search": lambda args: calls.append(("search", args)) or "one match",
        },
    )

    assert result.ok is True
    assert result.status == "completed"
    assert result.value == {"variants": []}
    assert result.response_id == "resp_2"
    assert calls == [
        ("read_table", {"table_id": "T1"}),
        ("search", {"query": "carrier"}),
    ]

    assert len(transport.calls) == 2
    first_url, first_request = transport.calls[0]
    second_url, second_request = transport.calls[1]
    assert first_url == (
        "https://example.services.ai.azure.com/openai/v1/responses?api-version=v1"
    )
    assert second_url == first_url
    assert first_request["headers"]["api-key"] == "secret"
    assert "Authorization" not in first_request["headers"]
    assert first_request["json"]["reasoning"] == {"effort": "max"}
    assert first_request["json"]["input"] == "paper overview"
    assert "previous_response_id" not in first_request["json"]

    continuation = second_request["json"]
    assert continuation["previous_response_id"] == "resp_1"
    assert continuation["instructions"] == first_request["json"]["instructions"]
    assert continuation["tools"] == first_request["json"]["tools"]
    assert continuation["text"] == first_request["json"]["text"]
    assert continuation["max_output_tokens"] == 4096
    assert continuation["store"] is True
    assert continuation["input"] == [
        {
            "type": "function_call_output",
            "call_id": "call_a",
            "output": '{"rows":[1]}',
        },
        {
            "type": "function_call_output",
            "call_id": "call_b",
            "output": "one match",
        },
    ]

    assert result.usage.input_tokens == 150
    assert result.usage.cached_input_tokens == 30
    assert result.usage.cache_write_input_tokens == 10
    assert result.usage.output_tokens == 52
    assert result.usage.reasoning_tokens == 38
    assert result.usage.total_tokens == 202
    assert result.telemetry.api_calls == 2
    assert result.telemetry.rounds == 2
    assert result.telemetry.retries == 0
    assert result.telemetry.latency_seconds == pytest.approx(0.5)
    assert result.telemetry.tool_latency_seconds == pytest.approx(0.5)
    assert result.telemetry.tool_calls_requested == 2
    assert result.telemetry.tool_calls_executed == 2
    assert result.telemetry.tool_errors == 0


@pytest.mark.parametrize("effort", REASONING_EFFORTS)
def test_every_reasoning_effort_is_sent_exactly_without_alias(effort):
    transport = FakeTransport([_message_response("resp", {"ok": True})])

    result = _run(_client(transport), reasoning_effort=effort)

    assert result.ok is True
    assert transport.calls[0][1]["json"]["reasoning"] == {"effort": effort}


@pytest.mark.parametrize(
    "model",
    ["azure_ai/gpt-5.6-sol", "openai/gpt-5.6-sol", "gpt-5.6-sol"],
)
def test_provider_prefix_is_removed_from_wire_model(model):
    transport = FakeTransport([_message_response("resp", {"ok": True})])

    result = _run(_client(transport), model=model)

    assert result.ok is True
    assert transport.calls[0][1]["json"]["model"] == "gpt-5.6-sol"


def test_rejects_reasoning_effort_outside_exact_allowlist():
    transport = FakeTransport([])

    with pytest.raises(ValueError, match="reasoning_effort"):
        _run(_client(transport), reasoning_effort="minimal")

    assert transport.calls == []


def test_openai_endpoint_uses_bearer_auth():
    transport = FakeTransport([_message_response("resp", {"ok": True})])
    client = ResponsesAPIClient.openai(
        "openai-secret",
        transport=transport,
        backoff_s=0,
    )

    assert _run(client).ok is True

    url, request = transport.calls[0]
    assert url == "https://api.openai.com/v1/responses"
    assert request["headers"]["Authorization"] == "Bearer openai-secret"
    assert "api-key" not in request["headers"]


def test_unknown_malformed_and_failing_tools_return_safe_matching_outputs():
    transport = FakeTransport(
        [
            _tool_response(
                "resp_1",
                [
                    ("unknown", "delete_everything", "{}"),
                    ("malformed", "read_table", "not json"),
                    ("failure", "read_table", "{}"),
                ],
            ),
            _message_response("resp_2", {"variants": []}),
        ]
    )

    def fail(_args):
        raise RuntimeError("table unavailable")

    result = _run(_client(transport), tool_handlers={"read_table": fail})

    assert result.ok is True
    assert result.telemetry.tool_calls_requested == 3
    assert result.telemetry.tool_calls_executed == 1
    assert result.telemetry.tool_errors == 3
    outputs = transport.calls[1][1]["json"]["input"]
    assert [item["call_id"] for item in outputs] == [
        "unknown",
        "malformed",
        "failure",
    ]
    decoded = [json.loads(item["output"]) for item in outputs]
    assert [item["error"]["type"] for item in decoded] == [
        "unknown_tool",
        "malformed_arguments",
        "tool_execution_error",
    ]
    assert "delete_everything" in decoded[0]["error"]["message"]
    assert "RuntimeError" in decoded[2]["error"]["message"]


def test_tool_outputs_are_bounded_per_call_and_in_total():
    transport = FakeTransport(
        [
            _tool_response(
                "resp_1",
                [
                    ("a", "large", "{}"),
                    ("b", "large", "{}"),
                ],
            ),
            _message_response("resp_2", {"ok": True}),
        ]
    )
    result = _run(
        _client(transport),
        tool_handlers={"large": lambda _args: "x" * 100},
        limits=ResponsesLoopLimits(
            max_rounds=3,
            max_tool_calls=2,
            max_tool_output_chars=20,
            max_total_tool_output_chars=25,
        ),
    )

    assert result.ok is True
    outputs = transport.calls[1][1]["json"]["input"]
    assert [len(item["output"]) for item in outputs] == [20, 5]
    assert outputs[0]["output"].endswith("...[truncated]")
    assert result.telemetry.tool_output_chars == 25
    assert result.telemetry.truncated_tool_outputs == 2


def test_tool_call_limit_stops_before_executing_batch():
    transport = FakeTransport(
        [
            _tool_response(
                "resp_1",
                [("a", "tool", "{}"), ("b", "tool", "{}")],
            )
        ]
    )
    executed = []

    result = _run(
        _client(transport),
        tool_handlers={"tool": lambda args: executed.append(args)},
        limits=ResponsesLoopLimits(max_tool_calls=1),
    )

    assert result.ok is False
    assert result.status == "tool_call_limit_exceeded"
    assert executed == []
    assert len(transport.calls) == 1


def test_round_limit_stops_without_running_unreturnable_tools():
    transport = FakeTransport([_tool_response("resp_1", [("a", "tool", "{}")])])
    executed = []

    result = _run(
        _client(transport),
        tool_handlers={"tool": lambda args: executed.append(args)},
        limits=ResponsesLoopLimits(max_rounds=1),
    )

    assert result.ok is False
    assert result.status == "round_limit_exceeded"
    assert executed == []


@pytest.mark.parametrize("status_code", [408, 409, 429, 500, 503, 599])
def test_retries_only_retryable_http_statuses(status_code):
    transport = FakeTransport(
        [
            FakeResponse(status_code=status_code, payload={}, text="try later"),
            _message_response("resp", {"ok": True}),
        ]
    )
    sleeps = []
    client = _client(
        transport,
        max_retries=1,
        backoff_s=0.5,
        sleep_fn=sleeps.append,
    )

    result = _run(client)

    assert result.ok is True
    assert result.telemetry.api_calls == 2
    assert result.telemetry.retries == 1
    assert sleeps == [0.5]


@pytest.mark.parametrize("status_code", [300, 400, 401, 403, 404, 422, 600])
def test_does_not_retry_nonretryable_http_statuses(status_code):
    transport = FakeTransport(
        [FakeResponse(status_code=status_code, payload={}, text="bad request")]
    )

    result = _run(_client(transport, max_retries=3))

    assert result.ok is False
    assert result.status == "http_error"
    assert result.telemetry.api_calls == 1
    assert result.telemetry.retries == 0


def test_retries_transport_errors_but_not_bad_success_json():
    transport = FakeTransport(
        [requests.ConnectionError("offline"), _message_response("resp", {"ok": 1})]
    )
    result = _run(_client(transport, max_retries=1))

    assert result.ok is True
    assert result.telemetry.retries == 1

    invalid_json_transport = FakeTransport(
        [FakeResponse(payload=ValueError("not json"), text="not-json")]
    )
    invalid = _run(_client(invalid_json_transport, max_retries=3))
    assert invalid.ok is False
    assert invalid.status == "protocol_error"
    assert invalid.telemetry.api_calls == 1
    assert invalid.telemetry.retries == 0


def test_retry_exhaustion_returns_transport_or_http_error():
    transport = FakeTransport(
        [requests.Timeout("slow"), requests.Timeout("still slow")]
    )
    result = _run(_client(transport, max_retries=1))

    assert result.ok is False
    assert result.status == "transport_error"
    assert result.telemetry.api_calls == 2
    assert result.telemetry.retries == 1


def test_detects_incomplete_refusal_no_final_and_invalid_json():
    incomplete = FakeTransport(
        [
            FakeResponse(
                payload={
                    "id": "one",
                    "status": "incomplete",
                    "incomplete_details": {"reason": "max_output_tokens"},
                    "output": [],
                }
            )
        ]
    )
    refusal = FakeTransport(
        [
            FakeResponse(
                payload={
                    "id": "two",
                    "status": "completed",
                    "output": [
                        {
                            "type": "message",
                            "content": [
                                {"type": "refusal", "refusal": "Cannot comply"}
                            ],
                        }
                    ],
                }
            )
        ]
    )
    no_final = FakeTransport(
        [FakeResponse(payload={"id": "three", "status": "completed", "output": []})]
    )
    invalid_json = FakeTransport(
        [
            FakeResponse(
                payload={
                    "id": "four",
                    "status": "completed",
                    "output_text": "not JSON",
                }
            )
        ]
    )

    assert _run(_client(incomplete)).status == "incomplete"
    assert _run(_client(refusal)).status == "refusal"
    assert _run(_client(no_final)).status == "no_final_output"
    invalid = _run(_client(invalid_json))
    assert invalid.status == "invalid_json"
    assert invalid.output_text == "not JSON"


def test_missing_response_or_call_id_fails_protocol_safely():
    no_response_id = FakeTransport(
        [
            FakeResponse(
                payload={
                    "status": "completed",
                    "output": [
                        {
                            "type": "function_call",
                            "call_id": "call",
                            "name": "tool",
                            "arguments": "{}",
                        }
                    ],
                }
            )
        ]
    )
    no_call_id = FakeTransport(
        [
            FakeResponse(
                payload={
                    "id": "resp",
                    "status": "completed",
                    "output": [
                        {
                            "type": "function_call",
                            "name": "tool",
                            "arguments": "{}",
                        }
                    ],
                }
            )
        ]
    )

    assert _run(_client(no_response_id)).status == "protocol_error"
    assert _run(_client(no_call_id)).status == "protocol_error"


def test_repeated_call_id_is_rejected_before_duplicate_side_effect():
    transport = FakeTransport(
        [
            _tool_response(
                "resp_1",
                [("duplicate", "tool", "{}"), ("duplicate", "tool", "{}")],
            )
        ]
    )
    executed = []

    result = _run(
        _client(transport),
        tool_handlers={"tool": lambda args: executed.append(args)},
    )

    assert result.status == "protocol_error"
    assert "repeated" in result.error
    assert executed == []


def test_result_and_telemetry_have_json_ready_dicts_without_raw_by_default():
    transport = FakeTransport([_message_response("resp", {"ok": True})])

    result = _run(_client(transport))
    serialized = result.to_dict()

    assert serialized["usage"] == result.usage.to_dict()
    assert serialized["telemetry"] == result.telemetry.to_dict()
    assert "raw_responses" not in serialized
    assert (
        result.to_dict(include_raw_responses=True)["raw_responses"]
        == result.raw_responses
    )


def test_total_tokens_fall_back_to_input_plus_output():
    transport = FakeTransport(
        [
            _message_response(
                "resp",
                {"ok": True},
                usage={"prompt_tokens": 11, "completion_tokens": 4},
            )
        ]
    )

    result = _run(_client(transport))

    assert result.usage.input_tokens == 11
    assert result.usage.output_tokens == 4
    assert result.usage.total_tokens == 15
