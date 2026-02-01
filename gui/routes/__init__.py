"""GUI route modules for FastAPI."""

from gui.routes.settings import router as settings_router
from gui.routes.browser import router as browser_router

__all__ = ["settings_router", "browser_router"]
