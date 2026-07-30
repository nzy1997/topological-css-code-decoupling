"""Runtime metadata helpers for generated benchmark and validation outputs."""

from pathlib import Path
import platform
import subprocess
import sys

from sage.env import SAGE_VERSION


def source_revision(source_root):
    """Return the Git revision containing the source tree."""
    resolved_root = str(Path(source_root).resolve())
    completed = subprocess.run(
        ["git", "-C", resolved_root, "rev-parse", "HEAD"],
        capture_output=True,
        check=True,
        text=True,
    )
    return completed.stdout.strip()


def source_tree_dirty(source_root):
    """Return whether tracked files under the source root have local changes."""
    resolved_root = str(Path(source_root).resolve())
    completed = subprocess.run(
        [
            "git",
            "-C",
            resolved_root,
            "status",
            "--porcelain",
            "--untracked-files=no",
            "--",
            ".",
        ],
        capture_output=True,
        check=True,
        text=True,
    )
    return bool(completed.stdout.strip())


def runtime_metadata(*, adapter, suite, source_root):
    """Return metadata needed to interpret one result file.

    Args:
        adapter: Decoder adapter name recorded in benchmark metadata.
        suite: Benchmark suite name recorded in metadata.
        source_root: Repository root used when collecting benchmark metadata.

    Returns:
        dict: Environment and source metadata for one result file.
    """
    resolved_root = str(Path(source_root).resolve())
    return {
        "adapter": adapter,
        "suite": suite,
        "source_root": resolved_root,
        "source_revision": source_revision(resolved_root),
        "source_tree_dirty": source_tree_dirty(resolved_root),
        "sage_version": SAGE_VERSION,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
    }
