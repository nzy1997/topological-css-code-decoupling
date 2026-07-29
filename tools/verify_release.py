#!/usr/bin/env python3
"""Static release-boundary checks for the public snapshot."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


FORBIDDEN_PATH_PARTS = (
    "docs/release-checklist.md",
    "docs/superpowers",
    "skills",
    "python/RELEASE_BLOCKERS.md",
    "python/.sage-decoder-public-release",
    "python/CITATION.cff",
    "python/MANIFEST.sha256",
    "julia/ToricBuilder/docs/plans",
    "julia/ToricBuilder/example/testing_case_cache",
    "julia/ToricBuilder/example/reference/results_macmini.md",
    "julia/ToricBuilder/bench",
    "julia/ToricBuilder/src/matching",
    "julia/ToricBuilder/src/export/decode_bundle_export.jl",
    "julia/ToricBuilder/src/export/decoder_bundle_types.jl",
    "julia/ToricBuilder/docs/expert/Matching_decoder_note.md",
    "julia/ToricBuilder/docs/expert/local_solve_explanation.md",
    "julia/ToricBuilder/docs/expert/matching_maps_phi_psi_K_explanation.md",
    "julia/ToricBuilder/docs/formats/css_decode_bundle_format.md",
    "julia/ToricBuilder/docs/workflows/matching_decoder_workflow.md",
    "julia/ToricBuilder/example/generated",
    "julia/ToricBuilder/example/scripts/color_code_decoder.jl",
    "julia/ToricBuilder/example/scripts/export_css_decode_bundle.jl",
    "julia/ToricBuilder/example/scripts/interface_demo.jl",
    "julia/ToricBuilder/example/scripts/save_code.jl",
    "python/results/unitary_decouple_decoder/color_666_unitary_vs_bposd.csv",
    "python/results/unitary_decouple_decoder/color_666_unitary_vs_bposd.md",
    "python/results/unitary_decouple_decoder/color_666_unitary_vs_bposd.png",
)

FORBIDDEN_FILENAMES = {
    "cleanup_incomplete_decoupled_caches.jl",
    "read_testing_case_cache.jl",
    "save_debug_matrices.jl",
    "test_plot_area_comparison_annotations.jl",
    "bench_harness.jl",
}

SECRET_LIKE_NAMES = (
    ".env",
    ".netrc",
    "id_rsa",
    "id_ed25519",
)

FORBIDDEN_BYTES = (
    b"/Use" + b"rs/",
    b"git@github.com:nzy1997/topological-code" + b"-decoupling.git",
    b"https://github.com/nzy1997/topological-code" + b"-decoupling",
    b"link-will" + b"-be-inserted",
    b"to be confirmed " + b"before release",
    b"Software authors to be confirmed " + b"before release",
    b"local bb_code_instances notebook",
    b"julia_code/",
    b"build_matching_maps",
    b"solve_locally",
    b"export_css_decode_bundle",
    b"TToricDecoder",
    b"d26b037",
    b"14291af",
    b"5f90c32",
    b'"ttoric"',
)

REQUIRED_DIRS = (
    "python/isomorphism",
    "python/decoder_core",
    "python/unitary_decouple_decoder",
    "julia/ToricBuilder/src",
    "julia/ToricBuilder/results",
)

REQUIRED_FILES = (
    "python/README.md",
    "python/pyproject.toml",
    "julia/ToricBuilder/Project.toml",
    "julia/ToricBuilder/Manifest.toml",
    "julia/ToricBuilder/README.md",
    "julia/ToricBuilder/results/bb_code_decoupling_results.md",
    "julia/ToricBuilder/results/additional_bb_code_decoupling_results.md",
    "julia/ToricBuilder/results/decoding_benchmark.json",
    "julia/ToricBuilder/results/transported_cnot_support.json",
    "julia/ToricBuilder/example/scripts/plot_area_comparison.jl",
    "julia/ToricBuilder/example/scripts/plot_decoding_result_from_data.jl",
    "julia/ToricBuilder/example/scripts/transported_cnot_support.jl",
)

CONTENT_SCAN_EXEMPTIONS = {
    "tools/verify_release.py",
    "tests/release/test_release_metadata.py",
}


def repo_root() -> Path:
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(result.stdout.strip())


def candidate_paths(root: Path) -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    paths = (p.decode("utf-8") for p in result.stdout.split(b"\0") if p)
    return sorted(path for path in paths if (root / path).exists())


def is_binary(payload: bytes) -> bool:
    return b"\0" in payload[:4096]


def check_layout(root: Path) -> list[str]:
    errors: list[str] = []
    for path in REQUIRED_DIRS:
        if not (root / path).is_dir():
            errors.append(f"missing required directory: {path}")
    for path in REQUIRED_FILES:
        if not (root / path).is_file():
            errors.append(f"missing required file: {path}")
    return errors


def check_paths(paths: list[str]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        name = Path(path).name
        if path.endswith(".jls"):
            errors.append(f"serialized cache must not be tracked: {path}")
        if name in FORBIDDEN_FILENAMES:
            errors.append(f"forbidden development file: {path}")
        if name in SECRET_LIKE_NAMES or name.endswith((".pem", ".key", ".p12")):
            errors.append(f"secret-like filename: {path}")
        for forbidden in FORBIDDEN_PATH_PARTS:
            if path == forbidden or path.startswith(f"{forbidden}/"):
                errors.append(f"forbidden release path: {path}")
    return errors


def check_contents(root: Path, paths: list[str]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        if path in CONTENT_SCAN_EXEMPTIONS:
            continue
        full_path = root / path
        if not full_path.is_file():
            continue
        try:
            payload = full_path.read_bytes()
        except OSError as exc:
            errors.append(f"cannot read {path}: {exc}")
            continue
        if is_binary(payload):
            continue
        for marker in FORBIDDEN_BYTES:
            if marker in payload:
                errors.append(f"forbidden marker {marker.decode('utf-8')!r} in {path}")
    return errors


def main() -> int:
    root = repo_root()
    paths = candidate_paths(root)
    errors = check_layout(root)
    errors.extend(check_paths(paths))
    errors.extend(check_contents(root, paths))
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("release boundary: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
