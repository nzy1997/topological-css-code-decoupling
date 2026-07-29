"""Check that the published package roots stay independently importable."""

import ast
from importlib.util import find_spec
from pathlib import Path
import subprocess
import sys

import _bootstrap

import decoder_core


PUBLIC_PACKAGES = (
    "decoder_core",
    "isomorphism",
    "unitary_decouple_decoder",
)
assert find_spec("decoder_core.shot_streams") is None
assert find_spec("decoder_core.thresholds") is None
for name in (
    "estimate_pairwise_crossings",
    "find_or_create_color_666_shots",
    "matrix_hash",
    "read_shot_stream",
    "write_shot_stream",
):
    assert name not in decoder_core.__all__
    assert not hasattr(decoder_core, name)

for package_name in PUBLIC_PACKAGES:
    package_root = _bootstrap.REPO_ROOT / package_name
    for source_path in package_root.rglob("*.py"):
        ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))

probe = """
import sys
for package_name in ('decoder_core', 'isomorphism', 'unitary_decouple_decoder'):
    __import__(package_name)
"""
subprocess.run(
    [sys.executable, "-c", probe],
    cwd=Path(_bootstrap.REPO_ROOT),
    check=True,
)

print("test_public_package_boundaries: ok")
