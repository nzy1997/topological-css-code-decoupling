from __future__ import annotations

import subprocess
import sys
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def citation_authors(section: str) -> list[str]:
    citation = ROOT / "CITATION.cff"
    if not citation.is_file():
        return []
    lines = citation.read_text(encoding="utf-8").splitlines()
    authors: list[str] = []
    in_section = False
    in_authors = section == "authors"
    for line in lines:
        if line == f"{section}:":
            in_section = True
            in_authors = section == "authors"
            continue
        if in_section and section != "authors" and line.strip() == "authors:":
            in_authors = True
            continue
        if in_section and line and not line.startswith(" "):
            break
        if in_authors and line.startswith("  ") and not line.startswith("    ") and section != "authors":
            break
        stripped = line.strip()
        if in_authors and (
            stripped.startswith("family-names:") or stripped.startswith("- family-names:")
        ):
            family = line.split(":", 1)[1].strip().strip('"')
        elif in_authors and stripped.startswith("given-names:"):
            given = line.split(":", 1)[1].strip().strip('"')
            authors.append(f"{given} {family}")
    return authors


class ReleaseMetadataTests(unittest.TestCase):
    def test_public_release_boundary_is_minimal_and_reproducible(self) -> None:
        required_paths = (
            "julia/ToricBuilder/results/bb_code_decoupling_results.md",
            "julia/ToricBuilder/results/additional_bb_code_decoupling_results.md",
            "julia/ToricBuilder/results/decoding_benchmark.json",
            "julia/ToricBuilder/results/transported_cnot_support.json",
            "julia/ToricBuilder/example/scripts/transported_cnot_support.jl",
            "python/tests/sage/decoupling/test_paper_color_code_matrices.sage",
            "python/tests/sage/unitary/test_decoder_logical_failure.sage",
            "python/tests/sage/unitary/test_decoder_subdistance_error.sage",
        )
        forbidden_paths = (
            "docs/release-checklist.md",
            "docs/superpowers",
            "skills",
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
            "python/isomorphism",
            "python/unitary_decouple_decoder",
        )

        for path in required_paths:
            self.assertTrue((ROOT / path).is_file(), f"missing required release file: {path}")
        for path in forbidden_paths:
            self.assertFalse((ROOT / path).exists(), f"forbidden release path: {path}")

    def test_root_license_and_repository_identity_are_stable(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")

        self.assertIn("topological-css-code-decoupling", readme)
        self.assertIn("Decoupling 2D translation-invariant topological CSS codes", readme)
        self.assertIn("MIT License", license_text)
        self.assertIn("Copyright (c) 2026 The Authors", license_text)

    def test_software_author_order_is_stable(self) -> None:
        self.assertEqual(
            ["Zhongyi Ni", "Mingxin He"],
            citation_authors("authors"),
        )

    def test_preferred_paper_citation_author_order_is_stable(self) -> None:
        self.assertEqual(
            ["Yifei Wang", "Zhongyi Ni", "Mingxin He", "Jinguo Liu", "Yingfei Gu"],
            citation_authors("preferred-citation"),
        )

    def test_package_metadata_names_code_authors(self) -> None:
        python_metadata = (ROOT / "python" / "pyproject.toml").read_text(encoding="utf-8")
        julia_metadata = (ROOT / "julia" / "ToricBuilder" / "Project.toml").read_text(encoding="utf-8")
        python_author_placeholder = (
            "Software authors to be confirmed " + "before release"
        )

        self.assertIn('{ name = "Mingxin He" }', python_metadata)
        self.assertNotIn(python_author_placeholder, python_metadata)
        self.assertIn('authors = ["Zhongyi Ni', julia_metadata)

    def test_python_test_extra_contains_verification_tools(self) -> None:
        python_metadata = (ROOT / "python" / "pyproject.toml").read_text(encoding="utf-8")

        for dependency in ('"coverage>=7"', '"jupyter>=1"', '"nbconvert>=7"'):
            self.assertIn(dependency, python_metadata)

    def test_python_ldpc_version_is_reproducible(self) -> None:
        with (ROOT / "python" / "pyproject.toml").open("rb") as stream:
            metadata = tomllib.load(stream)
        extras = metadata["project"]["optional-dependencies"]

        for group in ("benchmark", "test"):
            ldpc_requirements = [
                dependency
                for dependency in extras[group]
                if dependency.startswith("ldpc")
            ]
            self.assertEqual(["ldpc==2.4.1"], ldpc_requirements)

    def test_release_docs_do_not_keep_blocker_language(self) -> None:
        for path in (
            ROOT / "README.md",
            ROOT / "python" / "README.md",
            ROOT / "julia" / "ToricBuilder" / "README.md",
        ):
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("release candidate", text.lower())
            self.assertNotIn("RELEASE_BLOCKERS.md", text)
            self.assertNotIn("license and software-author list are " + "confirmed", text)

    def test_claude_md_points_to_agents_md(self) -> None:
        agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
        claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")

        self.assertIn("AGENTS.md", claude)
        self.assertLess(len(claude.splitlines()), 12)
        self.assertLess(len(claude), len(agents) / 2)

    def test_makefile_declares_stable_command_surface(self) -> None:
        makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

        for target in (
            "test-python",
            "test-julia",
            "verify-manifest",
            "verify",
            "reproduce",
            "reproduce-full",
        ):
            self.assertIn(f"{target}:", makefile)

        for command in (
            "scripts/decoupling/reproduce_666_color_code.sage",
            "scripts/decoupling/reproduce_488_color_code.sage",
            "reproduce_bb_codes.sage --rows benchmark-01",
            "--output-json ../python-bb-benchmark-01.json",
            "tests/sage/decoupling/test_*.sage",
            "tests/sage/unitary/test_*.sage",
            "benchmark_paper_bb_family.sage --smoke",
            "julia/ToricBuilder/example/scripts/color_code.jl",
            "main(list=ab_list1[1:1]",
            "plot_area_comparison.jl",
            "plot_decoding_result_from_data.jl",
            "transported_cnot_support.jl",
        ):
            self.assertIn(command, makefile)
        self.assertNotIn("--inplace", makefile)

    def test_julia_environment_check_instantiates_and_loads_oscar(self) -> None:
        result = subprocess.run(
            ["make", "-n", "check-julia-env"],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )

        self.assertEqual(0, result.returncode, result.stderr)
        self.assertIn("--startup-file=no", result.stdout)
        self.assertIn("--project=julia/ToricBuilder", result.stdout)
        self.assertIn("Pkg.instantiate(; allow_autoprecomp=false)", result.stdout)
        self.assertIn("using Oscar", result.stdout)
        self.assertIn("Julia/Oscar environment: ok", result.stdout)

    def test_reproduction_docs_reference_existing_entrypoints(self) -> None:
        reproduction = (ROOT / "docs" / "reproduction.md").read_text(encoding="utf-8")

        required_paths = (
            "python/results/decoupling/666_color_code.ipynb",
            "python/results/decoupling/488_color_code.ipynb",
            "python/scripts/decoupling/reproduce_bb_codes.sage",
            "python/scripts/unitary_decouple_based_decoder/benchmark_paper_bb_family.sage",
            "julia/ToricBuilder/example/scripts/color_code.jl",
            "julia/ToricBuilder/example/scripts/decouple_bbcodes.jl",
            "julia/ToricBuilder/example/scripts/decouple_additional_bb_codes.jl",
            "julia/ToricBuilder/example/scripts/plot_area_comparison.jl",
            "julia/ToricBuilder/example/scripts/plot_decoding_result_from_data.jl",
            "julia/ToricBuilder/example/scripts/transported_cnot_support.jl",
            "julia/ToricBuilder/results/bb_code_decoupling_results.md",
            "julia/ToricBuilder/results/additional_bb_code_decoupling_results.md",
            "julia/ToricBuilder/results/decoding_benchmark.json",
            "julia/ToricBuilder/results/transported_cnot_support.json",
        )
        for path in required_paths:
            self.assertTrue((ROOT / path).exists(), path)
            self.assertIn(path, reproduction)

        for phrase in (
            "6.6.6",
            "4.8.8",
            "BB",
            "make reproduce",
            "make reproduce-full",
        ):
            self.assertIn(phrase, reproduction)

    def test_julia_degree_one_map_names_match_the_paper(self) -> None:
        legacy_degree_one = "ph" + "i_1"
        legacy_compact = "ph" + "i1"
        paths = (
            ROOT / "docs" / "reproduction.md",
            ROOT / "julia" / "ToricBuilder" / "docs" / "workflows" / "toric_form_workflow.md",
            ROOT / "julia" / "ToricBuilder" / "docs" / "expert" / "laurent_gaussian_degree_growth.md",
            ROOT / "julia" / "ToricBuilder" / "example" / "scripts" / "decouple_bbcodes.jl",
            ROOT / "julia" / "ToricBuilder" / "example" / "scripts" / "decouple_additional_bb_codes.jl",
            ROOT / "julia" / "ToricBuilder" / "example" / "scripts" / "transported_cnot_support.jl",
            ROOT / "julia" / "ToricBuilder" / "example" / "scripts" / "laurent_gaussian_degree_growth.jl",
            ROOT / "julia" / "ToricBuilder" / "results" / "bb_code_decoupling_results.md",
            ROOT / "julia" / "ToricBuilder" / "results" / "additional_bb_code_decoupling_results.md",
            ROOT / "julia" / "ToricBuilder" / "results" / "transported_cnot_support.json",
        )
        for path in paths:
            text = path.read_text(encoding="utf-8")
            self.assertNotIn(legacy_degree_one, text, path)
            self.assertNotIn(legacy_compact, text, path)

        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("### Julia/Oscar map notation", readme)
        self.assertIn("psi_1_inverse", readme)
        self.assertIn("result.psi_inverse.psi_1_inverse", readme)

    def test_ci_workflow_declares_release_sage_and_julia_jobs(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        codecov = (ROOT / ".codecov.yml").read_text(encoding="utf-8")

        for phrase in (
            "push:",
            "pull_request:",
            "release-checks:",
            "sagemath:",
            "julia:",
            "make verify-manifest",
            "sagemath/sagemath:10.8",
            "version: '1.12.5'",
            "Instantiate and test Julia project",
            "Pkg.instantiate(; allow_autoprecomp=false)",
            'include("julia/ToricBuilder/test/runtests.jl")',
            "import decoder_core, decoupling, unitary_decouple_based_decoder",
            "tests/sage/decoupling/test_*.sage",
            "tests/sage/unitary/test_*.sage",
            "benchmark_paper_bb_family.sage --smoke",
        ):
            self.assertIn(phrase, workflow)

        self.assertIn("python:", codecov)
        self.assertNotIn("julia:", codecov)
        self.assertEqual(2, workflow.count("timeout-minutes: 60"))

    def test_release_verifier_accepts_complete_release_metadata(self) -> None:
        result = subprocess.run(
            [sys.executable, "tools/verify_release.py"],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(0, result.returncode, result.stdout + result.stderr)
