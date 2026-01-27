"""Directory browser API routes."""

import os
from pathlib import Path

from fastapi import APIRouter, HTTPException

from gui.models import DirectoryEntry, DirectoryListResponse

router = APIRouter(prefix="/api/browse", tags=["browser"])


@router.get("", response_model=DirectoryListResponse)
async def browse_directory(path: str = "~") -> DirectoryListResponse:
    """Browse a directory and return its contents."""
    # Expand user home directory
    browse_path = Path(path).expanduser().resolve()

    # Security: Don't allow browsing outside of reasonable locations
    # Allow home directory and subdirectories, /tmp, and relative paths from cwd
    allowed_roots = [
        Path.home(),
        Path("/tmp"),
        Path.cwd(),
    ]

    is_allowed = any(
        browse_path == root or root in browse_path.parents for root in allowed_roots
    )

    if not is_allowed:
        raise HTTPException(
            status_code=403, detail="Access to this directory is not allowed"
        )

    if not browse_path.exists():
        raise HTTPException(status_code=404, detail="Directory not found")

    if not browse_path.is_dir():
        raise HTTPException(status_code=400, detail="Path is not a directory")

    entries = []
    try:
        for item in sorted(
            browse_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())
        ):
            # Skip hidden files (optional)
            if item.name.startswith("."):
                continue

            try:
                is_readable = os.access(item, os.R_OK)
                entries.append(
                    DirectoryEntry(
                        name=item.name,
                        path=str(item),
                        is_dir=item.is_dir(),
                        is_readable=is_readable,
                    )
                )
            except PermissionError:
                entries.append(
                    DirectoryEntry(
                        name=item.name,
                        path=str(item),
                        is_dir=item.is_dir(),
                        is_readable=False,
                    )
                )
    except PermissionError:
        raise HTTPException(status_code=403, detail="Permission denied")

    # Calculate parent path
    parent_path = None
    if browse_path != Path.home() and browse_path.parent != browse_path:
        parent_path = str(browse_path.parent)

    # Check if we can create directories here
    can_create = os.access(browse_path, os.W_OK)

    return DirectoryListResponse(
        current_path=str(browse_path),
        parent_path=parent_path,
        entries=entries,
        can_create=can_create,
    )


@router.post("/create")
async def create_directory(path: str):
    """Create a new directory."""
    new_path = Path(path).expanduser().resolve()

    # Security check
    allowed_roots = [Path.home(), Path("/tmp"), Path.cwd()]
    is_allowed = any(
        root in new_path.parents or new_path == root for root in allowed_roots
    )

    if not is_allowed:
        raise HTTPException(
            status_code=403, detail="Cannot create directory in this location"
        )

    if new_path.exists():
        raise HTTPException(status_code=400, detail="Path already exists")

    try:
        new_path.mkdir(parents=True, exist_ok=True)
        return {"status": "ok", "path": str(new_path)}
    except PermissionError:
        raise HTTPException(status_code=403, detail="Permission denied")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/shortcuts")
async def get_directory_shortcuts():
    """Get common directory shortcuts."""
    shortcuts = []

    # Home directory
    home = Path.home()
    shortcuts.append({"name": "Home", "path": str(home), "icon": "home"})

    # Desktop (if exists)
    desktop = home / "Desktop"
    if desktop.exists():
        shortcuts.append({"name": "Desktop", "path": str(desktop), "icon": "desktop"})

    # Documents (if exists)
    documents = home / "Documents"
    if documents.exists():
        shortcuts.append(
            {"name": "Documents", "path": str(documents), "icon": "folder"}
        )

    # Current working directory
    cwd = Path.cwd()
    shortcuts.append({"name": "Project Root", "path": str(cwd), "icon": "code"})

    # Default output directory
    default_output = cwd / "output"
    shortcuts.append(
        {"name": "Default Output", "path": str(default_output), "icon": "output"}
    )

    return {"shortcuts": shortcuts}
