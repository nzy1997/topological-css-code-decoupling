"""Runtime metadata helpers for generated benchmark and validation outputs."""

from pathlib import Path
import platform
import subprocess
import sys

from sage.env import SAGE_VERSION


def source_revision(source_root):
    """Return the Git revision for one source tree.

    Args:
        source_root: Repository root used when collecting benchmark metadata.

    Returns:
        object: The Git revision for one source tree.
    """
    # Keep syndrome normalization separate from the matching solve.
    resolved_root = str(Path(source_root).resolve())
    completed = subprocess.run(
        ["git", "-C", resolved_root, "rev-parse", "HEAD"],
        capture_output=True,
        check=True,
        text=True,
    )
    return completed.stdout.strip()


def runtime_metadata(*, adapter, suite, source_root):
    """Return metadata needed to interpret one result file.

    Args:
        adapter: Decoder adapter name recorded in benchmark metadata.
        suite: Benchmark suite name recorded in metadata.
        source_root: Repository root used when collecting benchmark metadata.

    Returns:
        float: Metadata needed to interpret one result file.
    """
    # Keep syndrome normalization separate from the matching solve.
    resolved_root = str(Path(source_root).resolve())
    return {
        "adapter": adapter,
        "suite": suite,
        "source_root": resolved_root,
        "source_revision": source_revision(resolved_root),
        "sage_version": SAGE_VERSION,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
    }
