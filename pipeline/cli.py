"""Command-line interface for the Gene Variant Fetcher toolkit.

This lightweight Typer application exposes a stable entrypoint so the package
can be installed and invoked via the ``gvf`` console script declared in
``pyproject.toml``.
"""

from importlib import metadata
from typing import Optional

import typer

app = typer.Typer(help="Command-line utilities for the Gene Variant Fetcher pipeline.")


def _get_version() -> Optional[str]:
    """Return the installed package version if available."""

    try:
        return metadata.version("genevariantfetcher")
    except metadata.PackageNotFoundError:
        return None


@app.command()
def version() -> None:
    """Print the installed package version."""

    pkg_version = _get_version()
    if pkg_version:
        typer.echo(pkg_version)
    else:
        typer.echo("unknown", err=True)


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context) -> None:
    """Show help when no subcommand is provided."""

    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


if __name__ == "__main__":
    app()
