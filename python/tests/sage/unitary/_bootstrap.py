"""Test bootstrap for Sage scripts executed from tests/sage/unitary."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
repo_root = str(REPO_ROOT)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
