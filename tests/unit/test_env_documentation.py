"""Keep every behavior-changing environment variable on a documented surface.

GVF has typed settings plus direct environment reads in runtime modules and
operator scripts.  An undocumented exported value can silently alter a run, so
every direct read must appear in ``.env.example`` or ``config/settings.py``.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCANNED_PACKAGES = (
    "pipeline",
    "harvesting",
    "config",
    "utils",
    "gene_literature",
    "cli",
    "scripts",
)

# Standard-library/toolchain configuration read but not defined by GVF.
EXTERNAL = frozenset(
    {
        "HOME",
        "PATH",
        "PWD",
        "TMPDIR",
        "USER",
        "USERPROFILE",
        "LOCALAPPDATA",
        "APPDATA",
        "XDG_CONFIG_HOME",
        "XDG_CACHE_HOME",
        "PYTHONPATH",
        "PYTHONHASHSEED",
        "PYTHONDONTWRITEBYTECODE",
        "VIRTUAL_ENV",
        "CI",
        "GITHUB_ACTIONS",
        "TERM",
        "COLUMNS",
        "NO_COLOR",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "NO_PROXY",
        "SSL_CERT_FILE",
        "REQUESTS_CA_BUNDLE",
        "LITELLM_LOG",
    }
)

ENV_NAME = re.compile(r"^[A-Z][A-Z0-9_]{2,}$")
ENV_HELPERS = frozenset({"get_env_bool", "get_env_int", "_env_float"})


def _module_env_constants(tree: ast.Module) -> dict[str, str]:
    """Resolve module-level ``NAME = 'ENV_VAR'`` indirection."""
    constants: dict[str, str] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not isinstance(node.value, ast.Constant) or not isinstance(
            node.value.value, str
        ):
            continue
        if not ENV_NAME.match(node.value.value):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name):
                constants[target.id] = node.value.value
    return constants


def _env_reads(tree: ast.Module) -> set[str]:
    constants = _module_env_constants(tree)
    found: set[str] = set()

    def resolve(node: ast.AST) -> str | None:
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value if ENV_NAME.match(node.value) else None
        if isinstance(node, ast.Name):
            return constants.get(node.id)
        return None

    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and node.args:
            func = node.func
            attr = getattr(func, "attr", None)
            if attr in {"get", "getenv", "pop"}:
                value = ast.unparse(func.value) if hasattr(func, "value") else ""
                if "environ" in value or value == "os":
                    if resolved := resolve(node.args[0]):
                        found.add(resolved)
            callee = attr if attr else getattr(func, "id", None)
            if callee in ENV_HELPERS:
                if resolved := resolve(node.args[0]):
                    found.add(resolved)
        if isinstance(node, ast.Subscript):
            target = ast.unparse(node.value) if hasattr(node, "value") else ""
            if "environ" in target:
                if resolved := resolve(node.slice):
                    found.add(resolved)
    return found


def test_every_env_var_read_by_gvf_is_documented():
    env_example = (REPO / ".env.example").read_text(encoding="utf-8")
    settings = (REPO / "config" / "settings.py").read_text(encoding="utf-8").lower()
    undocumented: dict[str, str] = {}

    for package in SCANNED_PACKAGES:
        for path in sorted((REPO / package).rglob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for name in _env_reads(tree):
                if name in EXTERNAL or name in undocumented:
                    continue
                if name in env_example or name.lower() in settings:
                    continue
                undocumented[name] = str(path.relative_to(REPO))

    assert not undocumented, (
        "Environment variables read by GVF but documented nowhere.\n"
        "Add each to .env.example with its default and effect, or make it a "
        "typed setting:\n"
        + "\n".join(
            f"  {name:38s} first read in {path}"
            for name, path in sorted(undocumented.items())
        )
    )


def test_the_scanner_resolves_indirect_reads():
    tree = ast.parse((REPO / "utils" / "env_utils.py").read_text(encoding="utf-8"))
    assert "GVF_DISABLE_LOCAL_DATA" in _env_reads(tree)
