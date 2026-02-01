#!/usr/bin/env python3
"""
GeneVariantFetcher - Main Entry Point

This is the primary way to use GeneVariantFetcher. It launches a web-based GUI
that provides:
- Visual pipeline configuration and execution
- Real-time progress monitoring
- Job management (pause, resume, stop)
- Environment/settings management
- Directory browser for output paths

Usage:
    python main.py              # Launch GUI (default)
    python main.py --cli GENE   # CLI mode for scripting/automation

For CLI-only usage (scripting, automation):
    python -m cli.automated_workflow GENE --email EMAIL --output OUTPUT_DIR
"""

import argparse
import sys
from pathlib import Path

# Ensure project root is in path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def main():
    parser = argparse.ArgumentParser(
        description="GeneVariantFetcher - Extract genetic variants from literature",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py                    # Launch GUI (recommended)
  python main.py --port 8080        # GUI on custom port
  python main.py --cli SCN5A        # CLI mode for automation

The GUI is the recommended way to use GeneVariantFetcher. It provides:
- Easy configuration of API keys and settings
- Visual pipeline progress monitoring
- Job management with pause/resume support
- Directory browser for selecting output paths
        """,
    )

    # Mode selection
    parser.add_argument(
        "--cli",
        metavar="GENE",
        help="Run in CLI mode with specified gene (for scripting/automation)",
    )

    # GUI options
    parser.add_argument(
        "--host",
        default="localhost",
        help="Host to bind to (default: localhost)",
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

    # CLI mode options
    parser.add_argument(
        "--email",
        help="Email for NCBI (CLI mode only, or set NCBI_EMAIL env var)",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output directory (CLI mode only)",
    )

    args = parser.parse_args()

    if args.cli:
        # CLI mode - run automated workflow
        run_cli_mode(args)
    else:
        # GUI mode (default)
        run_gui_mode(args)


def run_gui_mode(args):
    """Launch the web-based GUI."""
    # Check dependencies
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

    # Load environment (don't validate yet - GUI will help configure)
    from dotenv import load_dotenv

    load_dotenv()

    # Check for incomplete jobs
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
        print(f"{'='*60}\n")

    # Print startup info
    import uvicorn
    import webbrowser

    url = f"http://{args.host}:{args.port}"
    if args.host == "0.0.0.0":
        url = f"http://localhost:{args.port}"

    print(f"\n{'='*60}")
    print("GeneVariantFetcher")
    print(f"{'='*60}")
    print(f"Starting GUI at {url}")
    print()
    print("Features:")
    print("  - Configure API keys and settings in the Settings tab")
    print("  - Browse and select output directories")
    print("  - Monitor pipeline progress in real-time")
    print("  - Pause, resume, and manage jobs")
    print()
    print("Press Ctrl+C to stop the server")
    print(f"{'='*60}\n")

    # Open browser after a short delay to allow server to start
    if not args.no_browser:
        import threading

        def open_browser_delayed():
            import time

            time.sleep(1.0)  # Wait for server to start
            webbrowser.open(url)

        threading.Thread(target=open_browser_delayed, daemon=True).start()

    # Run server
    uvicorn.run(
        "gui.server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


def run_cli_mode(args):
    """Run the pipeline in CLI mode (for automation/scripting)."""
    import os
    from dotenv import load_dotenv

    load_dotenv()

    gene = args.cli
    email = args.email or os.getenv("NCBI_EMAIL")
    output = args.output or f"./output/{gene}"

    if not email:
        print("Error: Email is required for CLI mode.")
        print("Either set NCBI_EMAIL environment variable or use --email flag.")
        sys.exit(1)

    # Import and run the automated workflow
    from cli.automated_workflow import automated_variant_extraction_workflow

    print(f"\n{'='*60}")
    print(f"Running CLI mode for gene: {gene}")
    print(f"Output directory: {output}")
    print(f"{'='*60}\n")

    automated_variant_extraction_workflow(
        gene_symbol=gene,
        email=email,
        output_dir=output,
    )


if __name__ == "__main__":
    main()
