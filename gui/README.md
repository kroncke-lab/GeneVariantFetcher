# GeneVariantFetcher GUI

A web-based interface for running variant extraction pipelines with:
- **Background execution** - Jobs run even if you close the browser
- **Checkpoint/resume** - Jobs can be resumed after computer shutdown
- **Real-time progress** - Watch extraction progress via WebSocket

## Quick Start

```bash
# Install GUI dependencies
pip install -r gui/requirements.txt

# Launch the GUI
python main.py
```

This will:
1. Start a local web server at http://localhost:8000
2. Open your browser to the GUI
3. Check for any incomplete jobs that can be resumed

## Usage

### Starting a New Job

1. Go to the **New Job** tab
2. Enter required fields:
   - **Gene Symbol** (e.g., BRCA1, SCN5A, KCNH2)
   - **Email** (for NCBI E-utilities)
   - **Output Directory** (where results will be saved)
3. Optionally adjust:
   - Max PMIDs and downloads
   - Synonym discovery
   - Advanced filtering settings
4. Click **Start Pipeline**

### Monitoring Progress

The **Active** tab shows:
- Current pipeline step
- Real-time statistics (PMIDs found, papers downloaded, etc.)
- Live log output

### Resuming Interrupted Jobs

If your computer shuts down or the process is interrupted:

1. Restart the GUI: `python main.py`
2. You'll see a notification about incomplete jobs
3. Go to **Jobs** tab and click **Resume** on the interrupted job

Jobs are checkpointed after each major step, so you won't lose progress.

### Viewing Results

After completion:
- Results are saved to `{output_dir}/{gene}/{timestamp}/`
- Click **Results** on a completed job to see summary
- Key files:
  - `{GENE}.db` - SQLite database with all extracted data
  - `{GENE}_workflow_summary.json` - Run statistics
  - `extractions/*.json` - Per-paper extraction results

## Architecture

```
gui/
├── checkpoint.py    # Job state persistence
├── worker.py        # Pipeline execution with checkpointing
├── server.py        # FastAPI backend + WebSocket
├── static/
│   └── index.html   # Single-page frontend
└── requirements.txt # GUI-specific dependencies
```

### Checkpoint Storage

Job checkpoints are stored in `~/.gvf_jobs/`:
```
~/.gvf_jobs/
├── brca1_20250112_143022_abc123/
│   └── checkpoint.json
└── scn5a_20250112_150000_def456/
    └── checkpoint.json
```

### API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/jobs` | List all jobs |
| POST | `/api/jobs` | Create new job |
| GET | `/api/jobs/{id}` | Get job details |
| POST | `/api/jobs/{id}/resume` | Resume interrupted job |
| POST | `/api/jobs/{id}/stop` | Request graceful stop |
| DELETE | `/api/jobs/{id}` | Delete job |
| GET | `/api/jobs/{id}/logs` | Get job logs |
| GET | `/api/jobs/{id}/results` | Get completed job results |
| WS | `/ws/{id}` | Real-time progress stream |

## Configuration

### Command Line Options

```bash
python main.py --help

Options:
  --host TEXT      Host to bind to (default: 127.0.0.1)
  --port INTEGER   Port to bind to (default: 8000)
  --no-browser     Don't open browser automatically
  --reload         Enable auto-reload for development
```

### Environment Variables

The GUI uses the same environment variables as the CLI:

```bash
# Required
OPENAI_API_KEY=...
NCBI_EMAIL=...

# Optional
NCBI_API_KEY=...        # Higher rate limits
ANTHROPIC_API_KEY=...   # For certain models
```

## Troubleshooting

### "Missing required packages"

Install GUI dependencies:
```bash
pip install -r gui/requirements.txt
```

### Job stuck at a step

1. Check the logs in the Active tab
2. Try stopping and resuming the job
3. Check that API keys are set correctly

### Port already in use

Use a different port:
```bash
python main.py --port 8080
```

### Can't connect from another machine

Bind to all interfaces:
```bash
python main.py --host 0.0.0.0
```

## Development

Run with auto-reload:
```bash
python main.py --reload
```

The frontend is a single HTML file with embedded CSS/JS for simplicity.
For production use, consider:
- Running behind nginx/Apache
- Adding authentication
- Using a proper job queue (Celery/RQ) for multi-user scenarios
