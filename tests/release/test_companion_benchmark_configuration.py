from __future__ import annotations

import importlib.machinery
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]


class RecordingBpOsdDecoder:
    def __init__(self, parity_check_matrix, **kwargs):
        self.parity_check_matrix = parity_check_matrix
        self.options = kwargs


def load_ldpc_adapter():
    decoder_core = types.ModuleType("decoder_core")
    decoder_core.__path__ = [str(ROOT / "python" / "decoder_core")]

    arrays = types.ModuleType("decoder_core.arrays")
    arrays.to_numpy_uint8 = lambda value: value

    sage = types.ModuleType("sage")
    sage_all = types.ModuleType("sage.all")
    sage_all.GF = lambda order: order
    sage_all.vector = lambda field, values: (field, values)
    sage.all = sage_all

    ldpc_dependency = types.ModuleType("ldpc")
    ldpc_dependency.BpOsdDecoder = RecordingBpOsdDecoder

    modules = {
        "decoder_core": decoder_core,
        "decoder_core.arrays": arrays,
        "sage": sage,
        "sage.all": sage_all,
        "ldpc": ldpc_dependency,
    }
    path = ROOT / "python" / "decoder_core" / "ldpc.py"
    spec = importlib.util.spec_from_file_location("decoder_core.ldpc", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load the BP-OSD adapter")
    module = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, modules):
        spec.loader.exec_module(module)
    module.require_bp_osd = lambda: RecordingBpOsdDecoder
    return module


class LaurentSymbol:
    def __add__(self, other):
        return self

    __radd__ = __add__

    def __mul__(self, other):
        return self

    __rmul__ = __mul__

    def __pow__(self, exponent):
        return self


def load_companion_benchmark():
    decoder_core = types.ModuleType("decoder_core")
    for name in (
        "DecoderStats",
        "bp_osd_decoder",
        "classify_attempt",
        "decode_with_bp_osd",
        "numpy_rng",
        "plot_logical_error_rates",
        "sample_data_error",
        "to_numpy_uint8",
    ):
        setattr(decoder_core, name, object())

    benchmark_metadata = types.ModuleType("decoder_core.benchmark_metadata")
    benchmark_metadata.runtime_metadata = object()

    decoupling = types.ModuleType("decoupling")
    decoupling.x = LaurentSymbol()
    decoupling.y = LaurentSymbol()

    unitary = types.ModuleType("unitary_decouple_based_decoder")
    unitary.UnitaryDecoupleBasedDecoder = object()

    modules = {
        "decoder_core": decoder_core,
        "decoder_core.benchmark_metadata": benchmark_metadata,
        "decoupling": decoupling,
        "unitary_decouple_based_decoder": unitary,
    }
    path = (
        ROOT
        / "python"
        / "scripts"
        / "unitary_decouple_based_decoder"
        / "benchmark_paper_bb_family.sage"
    )
    loader = importlib.machinery.SourceFileLoader(
        "companion_benchmark_paper_bb_family",
        str(path),
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    if spec is None:
        raise RuntimeError("failed to load the companion benchmark script")
    module = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, modules):
        loader.exec_module(module)
    return module


class CompanionBenchmarkConfigurationTests(unittest.TestCase):
    def test_bp_osd_decoder_defaults_to_serial_schedule(self) -> None:
        adapter = load_ldpc_adapter()

        decoder = adapter.bp_osd_decoder(
            "H_X",
            0.01,
            max_iter=100,
            osd_order=0,
        )

        self.assertEqual("serial", decoder.options["schedule"])

    def test_bp_osd_decoder_forwards_explicit_schedule(self) -> None:
        adapter = load_ldpc_adapter()

        decoder = adapter.bp_osd_decoder(
            "H_X",
            0.01,
            max_iter=100,
            osd_order=0,
            schedule="parallel",
        )

        self.assertEqual("parallel", decoder.options["schedule"])

    def test_companion_benchmark_cli_defaults_to_serial_schedule(self) -> None:
        benchmark = load_companion_benchmark()

        with patch.object(sys, "argv", ["benchmark_paper_bb_family.sage"]):
            args = benchmark.parse_args()

        self.assertEqual("serial", args.bp_schedule)

    def test_companion_benchmark_records_selected_schedule(self) -> None:
        benchmark = load_companion_benchmark()
        args = SimpleNamespace(
            seed=20260607,
            bp_method="minimum_sum",
            bp_schedule="parallel",
            ms_scaling_factor=0.625,
            osd_order=2,
            max_iter=288,
        )
        stats = {
            label: SimpleNamespace(
                shots=5,
                logical_failures=0,
                decode_failures=0,
                logical_error_rate=0.0,
                standard_error=0.0,
                elapsed_seconds=0.1,
            )
            for label in (benchmark.DECODER_LABEL, benchmark.BPOSD_LABEL)
        }

        rows = benchmark.result_rows(
            SimpleNamespace(correction_size=224),
            4,
            0.01,
            stats,
            args,
        )

        self.assertEqual(
            ["parallel", "parallel"],
            [row["bp_schedule"] for row in rows],
        )


if __name__ == "__main__":
    unittest.main()
