"""Test bootstrap helpers for Sage scripts."""

from pathlib import Path
import sys


# Sage sets sys.path[0] to the script directory when running a .sage file by path.
# Add the repository root so local package imports work without external PYTHONPATH.
REPO_ROOT = Path(__file__).resolve().parents[3]
repo_root = str(REPO_ROOT)

if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
