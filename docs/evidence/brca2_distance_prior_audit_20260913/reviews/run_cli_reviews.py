"""Bounded independent no-source CLI consultations from the same saved brief."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import time


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCRATCH = REPO / "tmp/brca2_distance_prior_cli"


def run(name):
    prompt = (HERE / "shared_prompt.txt").read_text()
    (HERE / f"{name}_prompt.txt").write_text(prompt)
    if name == "grok":
        args = [
            "/Users/kronckbm/.local/bin/grok",
            "--model",
            "grok-4.6",
            "--reasoning-effort",
            "high",
            "--verbatim",
            "--permission-mode",
            "plan",
            "--tools",
            "",
            "--no-subagents",
            "--disable-web-search",
            "--max-turns",
            "2",
            "--system-prompt-override",
            "You are Grok, a scientific-method reviewer. Review the provided aggregate context directly. No tools, files, web, delegation, plans or progress messages. Accept fixed user choices and check the arithmetic carefully. Return a complete concise final critique.",
            "--output-format",
            "json",
            "--prompt-file",
            str(HERE / "grok_prompt.txt"),
        ]
    else:
        args = [
            "/Users/kronckbm/.local/bin/agy",
            "--print",
            prompt,
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
        ]
    start = time.monotonic()
    meta = {
        "cli": args[0],
        "requested_model": args[2] if name == "grok" else "gemini-3.1-pro-high",
        "requested_effort": "high",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "prompt": f"{name}_prompt.txt",
        "scope": "Generic mathematical consultation with synthetic examples only; no internal empirical data or source-provenance findings transmitted.",
        "timeout_seconds": 270,
    }
    with (
        (HERE / f"{name}_raw.json").open("w") as output,
        (HERE / f"{name}_stderr.log").open("w") as error,
    ):
        process = subprocess.Popen(
            args, cwd=SCRATCH, stdout=output, stderr=error, start_new_session=True
        )
        try:
            meta["exit_code"] = process.wait(timeout=270)
            meta["status"] = "returned" if meta["exit_code"] == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            meta.update(status="timeout", exit_code=process.returncode)
    meta["elapsed_seconds"] = round(time.monotonic() - start, 3)
    meta["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (HERE / f"{name}_status.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta), flush=True)
    return meta


if __name__ == "__main__":
    SCRATCH.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run, ["agy", "grok"]))
