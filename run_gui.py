#!/usr/bin/env python3
"""
Launch script for GeneVariantFetcher GUI.

Usage:
    python run_gui.py                    # Start server at localhost:8000
    python run_gui.py --port 8080        # Custom port
    python run_gui.py --host 0.0.0.0     # Allow external connections

The GUI provides:
- Web interface for running variant extraction pipelines
- Background job execution (survives browser close)
- Checkpoint/resume for interrupted jobs
- Real-time progress updates via WebSocket
"""

import argparse
import sys
import webbrowser
from pathlib import Path

# Ensure project root is in path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def check_dependencies():
    """Check that GUI dependencies are installed."""
    missing = []

    try:
        import fastapi
    except ImportError:
        missing.append("fastapi")

    try:
        import uvicorn
    except ImportError:
        missing.append("uvicorn")

    try:
        import websockets
    except ImportError:
        missing.append("websockets")

    if missing:
        print("Missing required packages for GUI:")
        print(f"  {', '.join(missing)}")
        print("\nInstall them with:")
        print(f"  pip install {' '.join(missing)}")
        print("\nOr install all GUI requirements:")
        print("  pip install -r gui/requirements.txt")
        sys.exit(1)


def check_env():
    """Check that required environment variables are set."""
    import os
    from dotenv import load_dotenv

    load_dotenv()

    warnings = []

    if not os.getenv("OPENAI_API_KEY"):
        warnings.append("OPENAI_API_KEY not set - extraction will fail")

    if not os.getenv("NCBI_EMAIL"):
        warnings.append("NCBI_EMAIL not set - literature fetching may be rate-limited")

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  - {w}")
        print()


def show_incomplete_jobs():
    """Check for and display incomplete jobs."""
    from gui.checkpoint import CheckpointManager

    manager = CheckpointManager()
    incomplete = manager.get_incomplete_jobs()

    if incomplete:
        print(f"\n{'='*60}")
        print(f"Found {len(incomplete)} incomplete job(s) that can be resumed:")
        print(f"{'='*60}")
        for job in incomplete:
            print(f"  - {job.job_id}")
            print(f"    Gene: {job.gene_symbol}")
            print(f"    Step: {job.current_step.value}")
            print(f"    Last updated: {job.updated_at}")
        print(f"{'='*60}\n")
        print("These jobs can be resumed from the web interface.\n")


def main():
    parser = argparse.ArgumentParser(
        description="Launch GeneVariantFetcher GUI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_gui.py                    # Start at localhost:8000
  python run_gui.py --port 8080        # Custom port
  python run_gui.py --no-browser       # Don't open browser automatically
  python run_gui.py --host 0.0.0.0     # Allow external connections
        """,
    )

    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind to (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to bind to (default: 8000)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development",
    )

    args = parser.parse_args()

    # Check dependencies
    check_dependencies()

    # Check environment
    check_env()

    # Show incomplete jobs
    show_incomplete_jobs()

    # Import uvicorn
    import uvicorn

    url = f"http://{args.host}:{args.port}"
    if args.host == "0.0.0.0":
        url = f"http://localhost:{args.port}"

    print(f"{'='*60}")
    print("GeneVariantFetcher GUI")
    print(f"{'='*60}")
    print(f"Starting server at {url}")
    print(f"Press Ctrl+C to stop the server")
    print(f"{'='*60}\n")

    # Open browser
    if not args.no_browser:
        webbrowser.open(url)

    # Run server
    uvicorn.run(
        "gui.server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
