"""
Web application module for scLucid.

This module provides a modern web interface for interactive single-cell analysis.

Architecture:
- FastAPI backend for REST API
- Vue 3 frontend for UI
- Celery for async task processing
- Redis for caching and task queue

Usage:
    from scLucid.web import launch_web_app

    launch_web_app(host="0.0.0.0", port=8000)
"""

__version__ = "0.1.0"

from scLucid.web.api.main import app
from scLucid.web.cli import launch_web_app

__all__ = ["app", "launch_web_app"]
