# Manifest Schema Specification

## Overview

The manifest system provides structured tracking of processing results across pipeline stages. Each stage produces a manifest file that documents what was processed, what succeeded, what failed, and what files were created.

This enables:
- **Resumability**: Identify which items need reprocessing after failures
- **Audit trail**: Track processing history and outcomes
- **Stage handoffs**: Pass structured data between pipeline stages
- **Error analysis**: Categorize and analyze failure patterns

## Schema Version

Current version: **1.0**

## File Location Convention

Manifests are stored alongside processed data:
```
output/
  GENE/
    PMID/
      ...files...
    manifest.json          # Stage manifest for this gene
```

## JSON Schema

```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "title": "GVF Pipeline Manifest",
  "type": "object",
  "required": ["schema_version", "stage", "created_at", "entries"],
  "properties": {
    "schema_version": {
      "type": "string",
      "description": "Schema version for forward compatibility",
      "pattern": "^\\d+\\.\\d+$",
      "examples": ["1.0"]
    },
    "stage": {
      "type": "string",
      "enum": ["download", "scout", "extract", "migrate"],
      "description": "Pipeline stage that produced this manifest"
    },
    "gene": {
      "type": "string",
      "description": "Gene symbol being processed (optional for multi-gene runs)",
      "examples": ["BRCA1", "TP53"]
    },
    "created_at": {
      "type": "string",
      "format": "date-time",
      "description": "ISO 8601 timestamp when manifest was created"
    },
    "entries": {
      "type": "array",
      "description": "Processing results for each item",
      "items": {
        "$ref": "#/definitions/entry"
      }
    }
  },
  "definitions": {
    "entry": {
      "type": "object",
      "required": ["pmid", "status", "timestamp"],
      "properties": {
        "pmid": {
          "type": "string",
          "description": "PubMed ID",
          "pattern": "^\\d+$"
        },
        "status": {
          "type": "string",
          "enum": ["SUCCESS", "FAILED", "SKIPPED", "PAYWALL", "CAPTCHA", "TIMEOUT"],
          "description": "Processing outcome"
        },
        "error_message": {
          "type": "string",
          "description": "Error details when status is not SUCCESS"
        },
        "files_created": {
          "type": "array",
          "items": {"type": "string"},
          "description": "Relative paths to files created during processing"
        },
        "timestamp": {
          "type": "string",
          "format": "date-time",
          "description": "ISO 8601 timestamp when this entry was processed"
        }
      }
    }
  }
}
```

## Status Codes

| Status | Description |
|--------|-------------|
| `SUCCESS` | Item processed successfully |
| `FAILED` | Processing failed (generic error) |
| `SKIPPED` | Item intentionally skipped (e.g., already processed) |
| `PAYWALL` | Content behind paywall, could not access |
| `CAPTCHA` | Blocked by CAPTCHA challenge |
| `TIMEOUT` | Operation timed out |

## Example Manifest

```json
{
  "schema_version": "1.0",
  "stage": "download",
  "gene": "BRCA1",
  "created_at": "2025-01-30T14:30:00Z",
  "entries": [
    {
      "pmid": "12345678",
      "status": "SUCCESS",
      "files_created": ["BRCA1/12345678/full_text.pdf", "BRCA1/12345678/metadata.json"],
      "timestamp": "2025-01-30T14:30:05Z"
    },
    {
      "pmid": "87654321",
      "status": "PAYWALL",
      "error_message": "Requires institutional access",
      "files_created": [],
      "timestamp": "2025-01-30T14:30:10Z"
    },
    {
      "pmid": "11111111",
      "status": "TIMEOUT",
      "error_message": "Request timed out after 30s",
      "files_created": [],
      "timestamp": "2025-01-30T14:30:45Z"
    }
  ]
}
```

## Usage Patterns

### Creating a Manifest

```python
from utils.manifest import Manifest, ManifestEntry, Status

manifest = Manifest(stage="download", gene="BRCA1")

# Add entries as processing completes
manifest.add_entry(ManifestEntry(
    pmid="12345678",
    status=Status.SUCCESS,
    files_created=["BRCA1/12345678/full_text.pdf"]
))

# Save atomically
manifest.save("output/BRCA1/manifest.json")
```

### Loading and Filtering

```python
manifest = Manifest.load("output/BRCA1/manifest.json")

# Get failed entries for retry
failed = [e for e in manifest.entries if e.status != Status.SUCCESS]

# Get entries by status
paywalled = manifest.get_by_status(Status.PAYWALL)
```

### Stage Handoff

```python
# Download stage produces manifest
download_manifest = Manifest.load("output/BRCA1/download_manifest.json")

# Scout stage consumes successful downloads
successful_pmids = [e.pmid for e in download_manifest.entries 
                   if e.status == Status.SUCCESS]
```

## Design Decisions

1. **Atomic writes**: Manifests are written to a temp file first, then renamed. This prevents corruption from interrupted writes.

2. **Append-friendly**: The `add_entry()` method allows incremental building of manifests during processing.

3. **Status enum**: Using an enum for status codes ensures consistency and enables IDE autocompletion.

4. **ISO timestamps**: All timestamps use ISO 8601 format for unambiguous parsing.

5. **Relative paths**: `files_created` uses relative paths for portability across environments.
