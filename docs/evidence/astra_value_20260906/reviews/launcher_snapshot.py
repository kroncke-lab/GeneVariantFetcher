"""Run source-contained CLI consultations; no extraction API calls."""

import concurrent.futures
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "docs/evidence/astra_value_20260906"
SCRATCH = ROOT / "validation_runs/astra_value_20260906"
PHASE = sys.argv[1]
PROMPT = EVIDENCE / "reviews" / f"{PHASE}_prompt.txt"
BRIEF = PROMPT.read_text()
COMMANDS = {
    "claude": [
        "/Users/kronckbm/.local/bin/claude",
        "--print",
        "--model",
        "sonnet",
        "--effort",
        "high",
        "--max-budget-usd",
        "3",
        "--tools",
        "",
        "--strict-mcp-config",
        "--safe-mode",
        "--no-session-persistence",
        "--output-format",
        "json",
    ],
    "grok": [
        "/Users/kronckbm/.local/bin/grok",
        "--prompt-file",
        str(PROMPT),
        "--verbatim",
        "--output-format",
        "json",
        "--max-turns",
        "1",
        "--reasoning-effort",
        "high",
        "--no-subagents",
        "--disable-web-search",
        "--tools",
        "",
        "--permission-mode",
        "plan",
        "--cwd",
        str(SCRATCH),
    ],
    "agy": [
        "/Users/kronckbm/.local/bin/agy",
        "--print",
        BRIEF,
        "--model",
        "gemini-3.1-pro-high",
        "--effort",
        "high",
        "--mode",
        "plan",
        "--sandbox",
        "--disable-slash-commands",
        "--output-format",
        "json",
        "--print-timeout",
        "4m",
    ],
}


def run(name):
    raw_path = EVIDENCE / "reviews" / f"{name}_{PHASE}_raw.json"
    if raw_path.exists():
        raise RuntimeError(f"Refusing to overwrite consultation: {raw_path}")
    try:
        completed = subprocess.run(
            COMMANDS[name],
            input=BRIEF if name == "claude" else None,
            cwd=SCRATCH,
            capture_output=True,
            text=True,
            timeout=270,
        )
        raw_path.write_text(completed.stdout)
        (EVIDENCE / "reviews" / f"{name}_{PHASE}_stderr.txt").write_text(
            completed.stderr
        )
        result = {"reviewer": name, "exit_code": completed.returncode}
        try:
            payload = json.loads(completed.stdout)
            value = payload.get(
                {"claude": "result", "grok": "text", "agy": "response"}[name]
            )
            if isinstance(value, dict):
                value = value.get("text")
            if isinstance(value, str):
                (EVIDENCE / "reviews" / f"{name}_{PHASE}.md").write_text(
                    value.rstrip() + "\n"
                )
                result["readable_result"] = True
            else:
                result["readable_result"] = False
        except (ValueError, TypeError):
            result["readable_result"] = False
        return result
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or b""
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        raw_path.write_text(output)
        return {"reviewer": name, "timeout_seconds": 270}


with concurrent.futures.ThreadPoolExecutor(3) as pool:
    futures = [pool.submit(run, name) for name in COMMANDS]
    results = []
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        results.append(result)
        print(json.dumps(result), flush=True)
    (EVIDENCE / "reviews" / f"{PHASE}_status.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
