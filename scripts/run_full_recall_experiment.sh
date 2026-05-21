#!/usr/bin/env bash
#
# Post-insttoken GVF recall re-extraction driver.
#
# What it launches:
#   1. Pre-flight checks for env vars, gold inputs, and local Python setup.
#   2. Optional cleanup of already-running gvf-run / re-extract processes.
#   3. Four-gene post-insttoken re-extraction through copied validation_runs/
#      artifacts and
#      scripts/run_insttoken_reextract_experiment.py at --max-pmids=10000.
#   4. Recall scoring plus rows/runs MAE reports through the Python driver.
#
# Data policy:
#   This script writes only ignored local artifacts under recall_metrics/ and
#   validation_runs/. Do not commit downloaded papers, SQLite DBs, logs, or .env.
#
# Usage:
#   bash scripts/run_full_recall_experiment.sh --dry-run
#   bash scripts/run_full_recall_experiment.sh --foreground
#   bash scripts/run_full_recall_experiment.sh --genes KCNH2,KCNQ1 --max-pmids 10000
#   bash scripts/run_full_recall_experiment.sh --kill-existing
#   bash scripts/run_full_recall_experiment.sh --status
#
# Detached mode uses screen when available. On machines without screen, the
# script falls back to foreground mode. On macOS, detached and foreground runs
# are wrapped in caffeinate -i when available.
#
set -euo pipefail

cd "$(dirname "$0")/.."
REPO_ROOT="$(pwd -P)"

EMAIL_CLI=""
MAX_PMIDS="${MAX_PMIDS:-10000}"
GENES="${GENES:-KCNH2,KCNQ1,RYR2,SCN5A}"
SCREEN_NAME="${SCREEN_NAME:-gvf-experiment}"
PYTHON_BIN="${PYTHON:-.venv/bin/python}"
DRY_RUN=0
STATUS_ONLY=0
FOREGROUND=0
KILL_EXISTING=0

usage() {
    sed -n '2,38p' "$0"
}

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

file_mtime() {
    local path="$1"
    if stat -f "%Sm" "$path" >/dev/null 2>&1; then
        stat -f "%Sm" "$path"
    elif stat -c "%y" "$path" >/dev/null 2>&1; then
        stat -c "%y" "$path"
    else
        date -r "$path" "+%Y-%m-%d %H:%M:%S" 2>/dev/null || echo "unknown"
    fi
}

source_env_if_present() {
    if [[ -f .env ]]; then
        set -a
        set +u
        # shellcheck disable=SC1091
        source .env
        set -u
        set +a
    else
        echo "  .env not found; relying on exported environment variables"
    fi
}

print_screen_status() {
    if ! have_cmd screen; then
        echo "screen not installed"
        return
    fi

    local output
    output="$(screen -ls 2>&1 || true)"
    if grep -q "$SCREEN_NAME" <<<"$output"; then
        grep "$SCREEN_NAME" <<<"$output"
    else
        echo "no screen named $SCREEN_NAME"
    fi
}

latest_non_dry_launch_log() {
    local candidate
    while IFS= read -r candidate; do
        [[ -z "$candidate" ]] && continue
        if ! grep -q "Dry run: True" "$candidate" 2>/dev/null; then
            echo "$candidate"
            return
        fi
    done < <(ls -t recall_metrics/_detached_logs/parallel_*full_reextract_*.log 2>/dev/null || true)
}

latest_run_dir_with_gene_logs() {
    local candidate
    while IFS= read -r candidate; do
        [[ -z "$candidate" ]] && continue
        if compgen -G "$candidate/logs/gvf_run_*.log" >/dev/null; then
            echo "$candidate"
            return
        fi
    done < <(ls -td recall_metrics/post_reextract_* 2>/dev/null || true)
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --email)         EMAIL_CLI="$2"; shift 2 ;;
        --max-pmids)     MAX_PMIDS="$2"; shift 2 ;;
        --genes)         GENES="$2"; shift 2 ;;
        --screen-name)   SCREEN_NAME="$2"; shift 2 ;;
        --python)        PYTHON_BIN="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=1; shift ;;
        --status)        STATUS_ONLY=1; shift ;;
        --foreground)    FOREGROUND=1; shift ;;
        --kill-existing) KILL_EXISTING=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)               echo "unknown arg: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# -----------------------------------------------------------------------------
# Status check (no kill, no launch)
# -----------------------------------------------------------------------------
if [[ "$STATUS_ONLY" -eq 1 ]]; then
    echo "=== screen ==="
    print_screen_status
    echo ""
    echo "=== running gvf-run processes ==="
    pgrep -af "cli gvf-run" 2>/dev/null || echo "none"
    echo ""
    echo "=== running insttoken orchestrators ==="
    pgrep -af "run_insttoken_reextract_experiment.py" 2>/dev/null || echo "none"
    echo ""
    echo "=== latest experiment log ==="
    LATEST="$(latest_non_dry_launch_log)"
    if [[ -n "$LATEST" ]]; then
        echo "$LATEST"
        tail -20 "$LATEST"
    else
        echo "none"
    fi
    echo ""
    echo "=== latest per-gene gvf_run logs (mtime, last [N/total]) ==="
    LATEST_DIR="$(latest_run_dir_with_gene_logs)"
    if [[ -n "$LATEST_DIR" ]]; then
        for g in KCNH2 KCNQ1 RYR2 SCN5A; do
            L="$LATEST_DIR/logs/gvf_run_${g}.log"
            if [[ -f "$L" ]]; then
                MTIME=$(file_mtime "$L")
                PROG=$(grep -oE "\[[0-9]+/[0-9]+\]" "$L" | tail -1 || true)
                echo "$g: ${PROG:-stage-start} mtime=$MTIME"
            fi
        done
    else
        echo "none"
    fi
    exit 0
fi

# -----------------------------------------------------------------------------
# Stage 1: pre-flight
# -----------------------------------------------------------------------------
echo "=== Stage 1: pre-flight ==="

source_env_if_present

EMAIL="${EMAIL_CLI:-${EMAIL:-${NCBI_EMAIL:-}}}"
if [[ -n "$EMAIL" ]]; then
    export NCBI_EMAIL="${NCBI_EMAIL:-$EMAIL}"
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
    if have_cmd python3.11; then
        PYTHON_BIN="$(command -v python3.11)"
    elif have_cmd python3; then
        PYTHON_BIN="$(command -v python3)"
    else
        echo "  MISSING: no Python interpreter found. Create .venv or set PYTHON=/path/to/python." >&2
        exit 6
    fi
    echo "  WARNING: .venv/bin/python not found; using $PYTHON_BIN"
fi
echo "  python: $PYTHON_BIN"

missing=0
for var in ANTHROPIC_API_KEY NCBI_EMAIL ELSEVIER_INSTTOKEN; do
    value="${!var:-}"
    if [[ -z "$value" ]]; then
        echo "  MISSING: $var" >&2
        missing=1
    else
        echo "  OK: $var (len=${#value})"
    fi
done
if [[ "$missing" -eq 1 ]]; then
    echo "Refusing to launch with missing required env vars." >&2
    echo "Create .env from .env.example or export the missing values." >&2
    exit 3
fi

for g in ${GENES//,/ }; do
    f="gene_variant_fetcher_gold_standard/normalized/${g}_recall_input.csv"
    if [[ ! -f "$f" ]]; then
        echo "  MISSING gold standard: $f" >&2
        exit 4
    fi
done
echo "  gold standards present for: $GENES"

existing_gvf="$(pgrep -af "cli gvf-run" 2>/dev/null || true)"
existing_orch="$(pgrep -af "run_insttoken_reextract_experiment.py" 2>/dev/null || true)"
if [[ -n "$existing_gvf$existing_orch" ]]; then
    echo ""
    echo "  existing GVF processes detected:"
    [[ -n "$existing_gvf" ]] && echo "$existing_gvf"
    [[ -n "$existing_orch" ]] && echo "$existing_orch"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "  dry-run: leaving existing processes untouched"
    elif [[ "$KILL_EXISTING" -eq 1 ]]; then
        echo "  --kill-existing set; stopping prior GVF processes"
        pkill -9 -f "cli gvf-run" 2>/dev/null || true
        pkill -9 -f "run_insttoken_reextract_experiment.py" 2>/dev/null || true
        pkill -9 -f "playwright/driver" 2>/dev/null || true
        pkill -9 -f "caffeinate -i" 2>/dev/null || true
        if have_cmd screen; then
            screen -S "$SCREEN_NAME" -X quit 2>/dev/null || true
        fi
        sleep 2
    else
        echo "Refusing to launch while GVF processes are already running." >&2
        echo "Use --status to inspect, --kill-existing to replace them, or --foreground after they finish." >&2
        exit 5
    fi
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
    if pgrep -f "cli gvf-run" >/dev/null 2>&1; then
        echo "  WARNING: gvf-run processes still alive after cleanup" >&2
        pgrep -af "cli gvf-run" >&2 || true
        exit 5
    fi
    echo "  clean: no conflicting gvf-run processes"
else
    echo "  dry-run: no cleanup or launch will run"
fi

# -----------------------------------------------------------------------------
# Stage 2: launch
# -----------------------------------------------------------------------------
LABEL="post_reextract_$(date -u +%Y%m%dT%H%M%SZ)"
LAUNCH_LOG="recall_metrics/_detached_logs/parallel_${MAX_PMIDS}_full_reextract_$(date -u +%Y%m%dT%H%M%SZ).log"
mkdir -p "$(dirname "$LAUNCH_LOG")"

CMD=(
    "$PYTHON_BIN" scripts/run_insttoken_reextract_experiment.py
    --email "$EMAIL"
    --full-reextract
    --parallel
    --max-pmids "$MAX_PMIDS"
    --genes "$GENES"
    --label "$LABEL"
)
if [[ "$DRY_RUN" -eq 1 ]]; then
    CMD+=(--dry-run)
fi

if [[ "$DRY_RUN" -eq 0 && "$FOREGROUND" -eq 0 ]] && ! have_cmd screen; then
    echo "  screen not installed; falling back to foreground mode"
    FOREGROUND=1
fi

echo ""
echo "=== Stage 2: launch ==="
echo "  label:      $LABEL"
echo "  genes:      $GENES"
echo "  max-pmids:  $MAX_PMIDS"
echo "  mode:       $([[ "$FOREGROUND" -eq 1 || "$DRY_RUN" -eq 1 ]] && echo foreground || echo detached-screen)"
echo "  screen:     $SCREEN_NAME"
echo "  log:        $LAUNCH_LOG"
echo ""
echo "  command:"
printf "    %q" "${CMD[@]}"
echo ""
echo ""

if [[ "$DRY_RUN" -eq 1 || "$FOREGROUND" -eq 1 ]]; then
    if have_cmd caffeinate; then
        caffeinate -i "${CMD[@]}" 2>&1 | tee "$LAUNCH_LOG"
    else
        "${CMD[@]}" 2>&1 | tee "$LAUNCH_LOG"
    fi
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo ""
        echo "  --dry-run mode: no live processes started."
    fi
    exit 0
fi

cmd_string="$(printf '%q ' "${CMD[@]}")"
repo_q="$(printf '%q' "$REPO_ROOT")"
log_q="$(printf '%q' "$LAUNCH_LOG")"
keep_awake_prefix=""
if have_cmd caffeinate; then
    keep_awake_prefix="caffeinate -i "
fi

screen -dmS "$SCREEN_NAME" bash -lc "
    set -euo pipefail
    cd $repo_q
    if [[ -f .env ]]; then
        set -a
        set +u
        source .env
        set -u
        set +a
    fi
    ${keep_awake_prefix}${cmd_string}2>&1 | tee $log_q
"

sleep 5

echo "=== launch verification ==="
echo "  screen sessions:"
print_screen_status
echo ""
echo "  child gvf-run processes:"
pgrep -af "cli gvf-run" 2>/dev/null || echo "  none yet (orchestrator may still be starting)"
echo ""
echo "  first 25 lines of launch log:"
head -25 "$LAUNCH_LOG" 2>/dev/null || true

cat <<EOF

=== Stage 3: monitor ===

Live tail (interactive):
    screen -r $SCREEN_NAME

Status check:
    bash scripts/run_full_recall_experiment.sh --status

Per-gene logs:
    tail -F recall_metrics/${LABEL}/logs/gvf_run_*.log

Final artifacts:
    recall_metrics/${LABEL}/summary.json
    recall_metrics/${LABEL}/report.md
    recall_metrics/mae/${LABEL}_rows/report.md
    recall_metrics/mae/${LABEL}_runs/report.md
    recall_metrics/${LABEL}/experiment_manifest.json

EOF
