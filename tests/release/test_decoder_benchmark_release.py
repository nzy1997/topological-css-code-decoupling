from __future__ import annotations

import hashlib
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "julia" / "ToricBuilder" / "results"
DATASET = RESULTS / "decoding_benchmark.json"
SOURCES = RESULTS / "decoder_benchmark_sources.json"


class DecoderBenchmarkReleaseTests(unittest.TestCase):
    def test_decoder_checkout_directory_cannot_redirect_recursive_delete(self) -> None:
        result = subprocess.run(
            [
                "make",
                "--dry-run",
                "prepare-decoder-benchmark",
                "DECODER_BENCHMARK_DIR=/",
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertNotIn("rm -rf /\n", result.stdout)
        self.assertIn(
            "rm -rf build/reproduction/decoder-benchmark",
            result.stdout,
        )

    def test_decoder_checkout_verifies_the_tag_commit(self) -> None:
        result = subprocess.run(
            ["make", "--dry-run", "prepare-decoder-benchmark"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertIn(
            "d3690ec0ab1018dcc4f49166715fecc5d3e6d5ab",
            result.stdout,
        )
        self.assertIn("git -C build/reproduction/decoder-benchmark rev-parse HEAD", result.stdout)

    def test_external_sources_are_public_and_immutable(self) -> None:
        payload = json.loads(SOURCES.read_text(encoding="utf-8"))
        repositories = payload["repositories"]

        self.assertEqual(1, payload["schema_version"])
        self.assertEqual("decoding_benchmark.json", payload["dataset"])
        self.assertEqual(
            "60ef2ab0d7c1d3c5aa9e214cc096872ef8af3614d7ec5a895ff27f1f15675905",
            payload["dataset_sha256"],
        )
        self.assertEqual(
            "9306524aada307e9310a6a981834d983ce987503",
            repositories["matching_wrapper"]["commit"],
        )
        self.assertEqual(
            ("paper-decoder-2026.1", "1eb67a8712d66ad5a6b16895c52b5a0e9089a4c5"),
            (
                repositories["tensorqec"]["tag"],
                repositories["tensorqec"]["commit"],
            ),
        )
        self.assertEqual(
            (
                "paper-decoder-benchmark-2026.2",
                "d3690ec0ab1018dcc4f49166715fecc5d3e6d5ab",
            ),
            (
                repositories["decoder_benchmarks"]["tag"],
                repositories["decoder_benchmarks"]["commit"],
            ),
        )
        for repository in repositories.values():
            self.assertTrue(repository["repository"].startswith("https://github.com/"))

    def test_normalized_dataset_matches_the_exported_bytes(self) -> None:
        sources = json.loads(SOURCES.read_text(encoding="utf-8"))
        digest = hashlib.sha256(DATASET.read_bytes()).hexdigest()
        self.assertEqual(sources["dataset_sha256"], digest)

        dataset = json.loads(DATASET.read_text(encoding="utf-8"))
        self.assertEqual("archived", dataset["dataset_class"])
        self.assertEqual("not_recorded", dataset["seed_status"])
        self.assertEqual(
            "not_recorded",
            dataset["provenance"]["environment"]["status"],
        )
        self.assertEqual(
            "profile likelihood",
            dataset["statistics"]["interval"],
        )
        self.assertEqual(
            "independent_depolarizing_pauli",
            dataset["provenance"]["noise_model"]["name"],
        )
        self.assertEqual(
            {
                "I": "1-p",
                "X": "p/3",
                "Y": "p/3",
                "Z": "p/3",
            },
            dataset["provenance"]["noise_model"]["single_qubit_probabilities"],
        )
        self.assertEqual(
            "ceil_remaining_per_worker_then_sum",
            dataset["provenance"]["stopping_rule"]["parallel_budget"],
        )
        self.assertEqual(
            "maximum_failures + workers - 1",
            dataset["provenance"]["stopping_rule"]["maximum_aggregate_failures"],
        )
        self.assertEqual(
            2040,
            dataset["logical_error"]["4"]["unitary"]["error_count"][0],
        )

    def test_code_parameters_and_unknown_distance_are_explicit(self) -> None:
        payload = json.loads(SOURCES.read_text(encoding="utf-8"))
        expected = {
            "4": (224, 6),
            "6": (504, 6),
            "8": (896, 6),
            "10": (1400, 6),
        }
        for toric_distance, (n, k) in expected.items():
            instance = payload["code_instances"][toric_distance]
            self.assertEqual(toric_distance, str(instance["toric_code_distance"]))
            self.assertEqual((n, k), (instance["n"], instance["k"]))
            self.assertEqual("not_recorded", instance["code_distance_status"])

    def test_make_and_docs_expose_archive_smoke_and_replot_workflows(self) -> None:
        makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
        reproduction = (ROOT / "docs" / "reproduction.md").read_text(encoding="utf-8")

        for target in (
            "verify-decoder-archive",
            "reproduce-decoder-smoke",
            "replot-decoder-benchmark",
        ):
            self.assertIn(f"{target}:", makefile)
            self.assertIn(f"make {target}", reproduction)
        for phrase in (
            "paper-decoder-2026.1",
            "paper-decoder-benchmark-2026.2",
            "profile likelihood",
            "unweighted",
            "matched prior",
            "not recorded",
        ):
            self.assertIn(phrase, reproduction)


if __name__ == "__main__":
    unittest.main()
