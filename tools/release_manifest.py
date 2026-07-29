#!/usr/bin/env python3
"""Write or verify the release SHA-256 manifest."""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path


MANIFEST = "MANIFEST.sha256"


def repo_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(result.stdout.strip())


def tracked_paths(root: Path) -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    paths = [p.decode("utf-8") for p in result.stdout.split(b"\0") if p]
    return sorted(path for path in paths if path != MANIFEST)


def manifest_lines(root: Path) -> list[str]:
    lines: list[str] = []
    for path in tracked_paths(root):
        digest = hashlib.sha256((root / path).read_bytes()).hexdigest()
        lines.append(f"{digest}  {path}\n")
    return lines


def write_manifest(root: Path) -> None:
    (root / MANIFEST).write_text("".join(manifest_lines(root)), encoding="utf-8")


def check_manifest(root: Path) -> int:
    manifest = root / MANIFEST
    if not manifest.is_file():
        print(f"{MANIFEST} is missing", file=sys.stderr)
        return 1
    expected = "".join(manifest_lines(root))
    actual = manifest.read_text(encoding="utf-8")
    if actual != expected:
        print(f"{MANIFEST} is stale; run tools/release_manifest.py --write", file=sys.stderr)
        return 1
    print(f"{MANIFEST}: ok")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="write MANIFEST.sha256")
    mode.add_argument("--check", action="store_true", help="verify MANIFEST.sha256")
    args = parser.parse_args()

    root = repo_root()
    if args.write:
        write_manifest(root)
        print(f"{MANIFEST}: wrote {len(tracked_paths(root))} entries")
        return 0
    return check_manifest(root)


if __name__ == "__main__":
    raise SystemExit(main())
