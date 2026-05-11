"""Regression tests for tolerant LLM JSON parsing."""

from utils.llm_utils import parse_llm_json_response


def test_parse_fenced_json_with_trailing_prose():
    raw = """```json
{
  "decision": "FAIL",
  "reason": "No patient-derived variants.",
  "confidence": 0.95
}
```

The paper is a pure mechanistic pharmacology study.
"""

    parsed = parse_llm_json_response(raw)

    assert parsed["decision"] == "FAIL"
    assert parsed["confidence"] == 0.95


def test_parse_json_after_non_json_braced_preface():
    raw = """I considered the rule {PASS only for patient data}.

{
  "decision": "PASS",
  "reason": "The abstract reports genotyped patients.",
  "confidence": 0.9
}
"""

    parsed = parse_llm_json_response(raw)

    assert parsed["decision"] == "PASS"
    assert parsed["confidence"] == 0.9
