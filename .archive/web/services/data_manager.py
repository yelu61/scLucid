"""
Data management service for web application.

Handles project data storage and retrieval.
"""

import logging
from pathlib import Path
from typing import Optional

import anndata as ad

log = logging.getLogger(__name__)

# Storage directory
DATA_DIR = Path.home() / ".sclucid" / "web_data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# In-memory cache for active projects
_project_cache: dict = {}


class DataManager:
    """
    Manages project data storage and retrieval.

    Projects are stored as .h5ad files in ~/.sclucid/web_data/
    """

    def __init__(self, data_dir: Path = DATA_DIR):
        """
        Initialize data manager.

        Args:
            data_dir: Directory to store project data
        """
        self.data_dir = data_dir
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def save_project(self, project_id: str, adata: ad.AnnData) -> Path:
        """
        Save project data.

        Args:
            project_id: Project identifier
            adata: AnnData object to save

        Returns:
            Path to saved file
        """
        file_path = self.data_dir / f"{project_id}.h5ad"
        adata.write(file_path)
        log.info(f"Saved project {project_id} to {file_path}")

        # Update cache
        _project_cache[project_id] = adata

        return file_path

    def load_project(self, project_id: str) -> Optional[ad.AnnData]:
        """
        Load project data.

        Args:
            project_id: Project identifier

        Returns:
            AnnData object or None if not found
        """
        # Check cache first
        if project_id in _project_cache:
            return _project_cache[project_id]

        # Load from disk
        file_path = self.data_dir / f"{project_id}.h5ad"

        if not file_path.exists():
            return None

        try:
            adata = ad.read(file_path)
            _project_cache[project_id] = adata
            log.info(f"Loaded project {project_id} from {file_path}")
            return adata
        except Exception as e:
            log.error(f"Error loading project {project_id}: {e}")
            return None

    def delete_project(self, project_id: str) -> bool:
        """
        Delete project data.

        Args:
            project_id: Project identifier

        Returns:
            True if deleted, False if not found
        """
        file_path = self.data_dir / f"{project_id}.h5ad"

        # Remove from cache
        if project_id in _project_cache:
            del _project_cache[project_id]

        if file_path.exists():
            file_path.unlink()
            log.info(f"Deleted project {project_id}")
            return True

        return False

    def list_projects(self) -> list:
        """
        List all projects.

        Returns:
            List of project IDs
        """
        projects = []
        for file_path in self.data_dir.glob("*.h5ad"):
            projects.append(file_path.stem)
        return projects


# Global data manager instance
_data_manager = DataManager()


def get_project_data(project_id: str) -> Optional[ad.AnnData]:
    """
    Get project data.

    Args:
        project_id: Project identifier

    Returns:
        AnnData object or None if not found
    """
    return _data_manager.load_project(project_id)


def update_project_data(project_id: str, adata: ad.AnnData):
    """
    Update project data.

    Args:
        project_id: Project identifier
        adata: Updated AnnData object
    """
    _data_manager.save_project(project_id, adata)


def create_project(project_id: str, adata: ad.AnnData) -> bool:
    """
    Create a new project.

    Args:
        project_id: Project identifier
        adata: AnnData object

    Returns:
        True if successful
    """
    try:
        _data_manager.save_project(project_id, adata)
        return True
    except Exception as e:
        log.error(f"Error creating project {project_id}: {e}")
        return False
