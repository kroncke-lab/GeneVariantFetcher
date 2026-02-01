"""Settings management API routes."""

import logging
import os
from pathlib import Path
from typing import Dict

from fastapi import APIRouter, HTTPException

from gui.models import EnvSettingsResponse, EnvSettingsUpdateRequest

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/settings", tags=["settings"])

# Path to .env file
ENV_FILE_PATH = Path(__file__).resolve().parent.parent.parent / ".env"

# Define which settings are configurable via the GUI
CONFIGURABLE_SETTINGS = {
    # API Keys (required)
    "OPENAI_API_KEY": {
        "label": "OpenAI API Key",
        "type": "password",
        "required": True,
        "group": "API Keys",
    },
    "ANTHROPIC_API_KEY": {
        "label": "Anthropic API Key",
        "type": "password",
        "required": False,
        "group": "API Keys",
    },
    "NCBI_EMAIL": {
        "label": "NCBI Email",
        "type": "email",
        "required": True,
        "group": "API Keys",
    },
    "NCBI_API_KEY": {
        "label": "NCBI API Key (optional)",
        "type": "password",
        "required": False,
        "group": "API Keys",
    },
    # Model Configuration
    "TIER2_MODEL": {
        "label": "Tier 2 Model",
        "type": "text",
        "default": "gpt-4o-mini",
        "group": "Models",
    },
    "TIER3_MODELS": {
        "label": "Tier 3 Models (comma-separated)",
        "type": "text",
        "default": "gpt-4o-mini,gpt-4o",
        "group": "Models",
    },
    # Pipeline Defaults
    "ENABLE_TIER1": {
        "label": "Enable Tier 1 (Keyword Filter)",
        "type": "checkbox",
        "default": "true",
        "group": "Pipeline",
    },
    "ENABLE_TIER2": {
        "label": "Enable Tier 2 (LLM Filter)",
        "type": "checkbox",
        "default": "true",
        "group": "Pipeline",
    },
    "TIER2_CONFIDENCE_THRESHOLD": {
        "label": "Filter Confidence Threshold",
        "type": "number",
        "default": "0.5",
        "group": "Pipeline",
    },
    # Literature Sources
    "USE_PUBMIND": {
        "label": "Use PubMind",
        "type": "checkbox",
        "default": "true",
        "group": "Sources",
    },
    "USE_PUBMED": {
        "label": "Use PubMed",
        "type": "checkbox",
        "default": "true",
        "group": "Sources",
    },
    "USE_EUROPEPMC": {
        "label": "Use Europe PMC",
        "type": "checkbox",
        "default": "false",
        "group": "Sources",
    },
    # Scout Configuration
    "SCOUT_ENABLED": {
        "label": "Enable DATA_ZONES (faster extraction)",
        "type": "checkbox",
        "default": "true",
        "group": "Scout",
    },
    "SCOUT_MIN_RELEVANCE": {
        "label": "Scout Min Relevance",
        "type": "number",
        "default": "0.3",
        "group": "Scout",
    },
    # Output Defaults
    "DEFAULT_OUTPUT_DIR": {
        "label": "Default Output Directory",
        "type": "directory",
        "default": "./output",
        "group": "Output",
    },
}


def parse_env_file(path: Path) -> Dict[str, str]:
    """Parse a .env file and return key-value pairs."""
    settings = {}
    if not path.exists():
        return settings

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue
            # Parse key=value
            if "=" in line:
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip()
                # Remove quotes if present
                if (value.startswith('"') and value.endswith('"')) or (
                    value.startswith("'") and value.endswith("'")
                ):
                    value = value[1:-1]
                settings[key] = value
    return settings


def write_env_file(path: Path, settings: Dict[str, str]):
    """Write settings to a .env file, preserving comments and structure."""
    existing_lines = []
    existing_keys = set()

    # Read existing file to preserve structure and comments
    if path.exists():
        with open(path, "r") as f:
            for line in f:
                original_line = line
                line = line.strip()
                if not line or line.startswith("#"):
                    existing_lines.append(original_line)
                elif "=" in line:
                    key, _, _ = line.partition("=")
                    key = key.strip()
                    existing_keys.add(key)
                    if key in settings:
                        # Update existing key
                        value = settings[key]
                        # Quote value if it contains spaces or special chars
                        if " " in value or "#" in value:
                            value = f'"{value}"'
                        existing_lines.append(f"{key}={value}\n")
                    else:
                        # Keep original line
                        existing_lines.append(original_line)

    # Add any new keys that weren't in the file
    new_keys = set(settings.keys()) - existing_keys
    if new_keys:
        if existing_lines and not existing_lines[-1].endswith("\n"):
            existing_lines.append("\n")
        existing_lines.append("\n# Added via GUI\n")
        for key in sorted(new_keys):
            value = settings[key]
            if " " in value or "#" in value:
                value = f'"{value}"'
            existing_lines.append(f"{key}={value}\n")

    # Write file
    with open(path, "w") as f:
        f.writelines(existing_lines)


@router.get("", response_model=EnvSettingsResponse)
async def get_settings():
    """Get current environment settings."""
    # Load from .env file
    file_settings = parse_env_file(ENV_FILE_PATH)

    # Also check current environment (might have been set externally)
    settings = {}
    for key, config in CONFIGURABLE_SETTINGS.items():
        # Priority: .env file > environment variable > default
        if key in file_settings:
            settings[key] = file_settings[key]
        elif key in os.environ:
            settings[key] = os.environ[key]
        else:
            settings[key] = config.get("default", "")

    return EnvSettingsResponse(
        settings=settings,
        env_file_exists=ENV_FILE_PATH.exists(),
        env_file_path=str(ENV_FILE_PATH),
    )


@router.get("/schema")
async def get_settings_schema():
    """Get the schema for configurable settings."""
    return {"settings": CONFIGURABLE_SETTINGS}


@router.post("")
async def update_settings(request: EnvSettingsUpdateRequest):
    """Update environment settings."""
    # Validate that only known settings are being updated
    unknown = set(request.settings.keys()) - set(CONFIGURABLE_SETTINGS.keys())
    if unknown:
        raise HTTPException(
            status_code=400, detail=f"Unknown settings: {', '.join(unknown)}"
        )

    # Load existing settings
    existing = parse_env_file(ENV_FILE_PATH)

    # Merge with new settings (only update non-empty values)
    for key, value in request.settings.items():
        if value or key not in existing:
            existing[key] = value

    # Write back to file
    write_env_file(ENV_FILE_PATH, existing)

    # Reload dotenv to update current process
    from dotenv import load_dotenv

    load_dotenv(ENV_FILE_PATH, override=True)

    logger.info(f"Settings updated: {list(request.settings.keys())}")

    return {"status": "ok", "updated": list(request.settings.keys())}


@router.get("/validate")
async def validate_settings():
    """Validate current settings and return any issues."""
    file_settings = parse_env_file(ENV_FILE_PATH)

    issues = []
    warnings = []

    # Check required settings
    if not file_settings.get("OPENAI_API_KEY") and not os.environ.get("OPENAI_API_KEY"):
        issues.append("OPENAI_API_KEY is required for variant extraction")

    if not file_settings.get("NCBI_EMAIL") and not os.environ.get("NCBI_EMAIL"):
        issues.append("NCBI_EMAIL is required for literature fetching")

    # Check for placeholder values
    placeholders = {"your_api_key", "changeme", "your-openai-api-key"}
    for key, value in file_settings.items():
        if value.lower() in placeholders:
            issues.append(f"{key} appears to be a placeholder value")

    # Warnings
    if not file_settings.get("NCBI_API_KEY") and not os.environ.get("NCBI_API_KEY"):
        warnings.append("NCBI_API_KEY not set - API requests may be rate-limited")

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
    }
