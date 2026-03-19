"""
Command-line interface for launching scLucid web application.
"""

import sys
import logging

log = logging.getLogger(__name__)


def launch_web_app(
    host: str = "0.0.0.0",
    api_port: int = 8000,
    dash_port: int = 8050,
    debug: bool = False,
):
    """
    Launch the scLucid web application.

    This starts both:
    1. FastAPI backend server (port 8000)
    2. Dash frontend server (port 8050)

    Args:
        host: Host to bind to
        api_port: Port for FastAPI backend
        dash_port: Port for Dash frontend
        debug: Enable debug mode

    Examples:
        >>> from scLucid.web import launch_web_app
        >>> launch_web_app()

        Or command line:
        $ sclucid web --host 0.0.0.0 --api-port 8000 --dash-port 8050
    """
    import subprocess
    import time

    # Import after logging setup
    from scLucid.web.api.main import app as api_app
    from scLucid.web.dash_app import QCDashApp
    import uvicorn

    log.info("Starting scLucid web application...")
    log.info(f"Backend API: http://{host}:{api_port}")
    log.info(f"Frontend Dashboard: http://{host}:{dash_port}")

    # Start FastAPI backend in a separate process
    api_process = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "uvicorn",
            "scLucid.web.api.main:app",
            "--host",
            host,
            "--port",
            str(api_port),
        ]
    )

    # Give API server time to start
    time.sleep(2)

    try:
        # Start Dash frontend in main process
        dash_app = QCDashApp(
            api_base_url=f"http://{host}:{api_port}",
            port=dash_port,
            debug=debug,
        )
        dash_app.run(host=host)

    except KeyboardInterrupt:
        log.info("Shutting down...")
    finally:
        # Cleanup API process
        api_process.terminate()
        api_process.wait()
