from __future__ import annotations

import hashlib
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def tracked_paths() -> set[str]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return {
        path.decode("utf-8")
        for path in result.stdout.split(b"\0")
        if path
    }


class ReleaseManifestTests(unittest.TestCase):
    def test_exported_source_layout_is_present(self) -> None:
        for package in (
            "decoder_core",
            "decoupling",
            "unitary_decouple_based_decoder",
        ):
            self.assertTrue((ROOT / "python" / package).is_dir())
        self.assertFalse((ROOT / "python" / "isomorphism").exists())
        self.assertFalse((ROOT / "python" / "unitary_decouple_decoder").exists())
        self.assertTrue((ROOT / "julia" / "ToricBuilder" / "src").is_dir())

    def test_exported_tree_excludes_private_artifacts(self) -> None:
        self.assertFalse((ROOT / "python" / "RELEASE_BLOCKERS.md").exists())
        self.assertFalse((ROOT / "julia" / "ToricBuilder" / "docs" / "plans").exists())
        non_build_caches = [
            path
            for path in ROOT.rglob("*.jls")
            if path.relative_to(ROOT).parts[:1] != ("build",)
        ]
        self.assertFalse(non_build_caches)
        self.assertFalse([path for path in tracked_paths() if path.endswith(".jls")])

    def test_manifest_matches_every_tracked_file_except_itself(self) -> None:
        manifest = ROOT / "MANIFEST.sha256"
        if not manifest.is_file():
            self.skipTest("MANIFEST.sha256 is frozen in Task 7")
        manifest_paths: set[str] = set()
        for line in manifest.read_text(encoding="utf-8").splitlines():
            digest, path = line.split("  ", 1)
            self.assertEqual(64, len(digest))
            self.assertEqual(digest, hashlib.sha256((ROOT / path).read_bytes()).hexdigest())
            manifest_paths.add(path)
        self.assertEqual(manifest_paths, tracked_paths() - {"MANIFEST.sha256"})

    def test_manifest_checker_accepts_the_tracked_tree(self) -> None:
        if not (ROOT / "MANIFEST.sha256").is_file():
            self.skipTest("MANIFEST.sha256 is frozen in Task 7")
        result = subprocess.run(
            [sys.executable, "tools/release_manifest.py", "--check"],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(0, result.returncode, result.stdout + result.stderr)
