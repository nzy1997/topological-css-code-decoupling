"""Check additional unitary-vs-BPOSD benchmark entry points."""

from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
import tempfile

import _bootstrap


REPO_ROOT = _bootstrap.REPO_ROOT
SCRIPT_488 = REPO_ROOT / "scripts" / "unitary" / "benchmark_color_488_unitary_vs_bposd.sage"
SCRIPT_BB = REPO_ROOT / "scripts" / "unitary" / "benchmark_bb_instances_unitary_vs_bposd.sage"

benchmark_488 = SourceFileLoader("benchmark_color_488_unitary_vs_bposd", str(SCRIPT_488)).load_module()
benchmark_bb = SourceFileLoader("benchmark_bb_instances_unitary_vs_bposd", str(SCRIPT_BB)).load_module()


def parsed_args(module, script_path, *args):
    """Parse script args without leaking argv between tests."""
    old_argv = sys.argv
    try:
        sys.argv = [str(script_path), *args]
        return module.parse_args()
    finally:
        sys.argv = old_argv


args_488 = parsed_args(benchmark_488, SCRIPT_488)
assert benchmark_488.parse_distances(args_488.distances) == [4, 6, 8]
assert benchmark_488.parse_probabilities(args_488.p_list) == [0.01, 0.02, 0.03, 0.05, 0.08, 0.1]
assert args_488.shots == 100000
assert args_488.output == "results/unitary_decouple_decoder/color_488_unitary_vs_bposd.csv"
assert args_488.plot_output == "results/unitary_decouple_decoder/color_488_unitary_vs_bposd.png"

decoder_488 = benchmark_488.color_488_decoder(2)
assert decoder_488.check_chain_relation()
assert decoder_488.num_x_checks == 4
assert decoder_488.num_qubits == 8

args_bb = parsed_args(benchmark_bb, SCRIPT_BB)
assert [row.row for row in benchmark_bb.parse_rows(args_bb.rows)] == [2, 3, 5]
assert benchmark_bb.parse_distances(args_bb.distances) == [3, 5]
assert args_bb.shots == 100000
assert args_bb.output == "results/unitary_decouple_decoder/bb_instances_unitary_vs_bposd.csv"

row_2 = benchmark_bb.parse_rows("2")[0]
decoder_bb = benchmark_bb.bb_decoder(row_2, 2)
assert decoder_bb.check_chain_relation()

rows = [
    {
        "code": "bb_row_2",
        "row": 2,
        "distance": 2,
        "n": decoder_bb.correction_size,
        "p": 0.05,
        "decoder": "unitary_decouple",
        "shots": 5,
        "logical_failures": 1,
        "decode_failures": 0,
        "logical_error_rate": 0.2,
        "standard_error": 0.1788854382,
        "elapsed_seconds": 0.5,
    },
    {
        "code": "bb_row_2",
        "row": 2,
        "distance": 2,
        "n": decoder_bb.correction_size,
        "p": 0.05,
        "decoder": "bposd",
        "shots": 5,
        "logical_failures": 0,
        "decode_failures": 0,
        "logical_error_rate": 0.0,
        "standard_error": 0.0,
        "elapsed_seconds": 0.7,
    },
]
with tempfile.TemporaryDirectory() as tmpdir:
    output = Path(tmpdir) / "bb.csv"
    md_path, plot_paths = benchmark_bb.write_outputs(rows, output, "BB Test")
    assert output.exists()
    assert md_path.exists()
    assert len(plot_paths) == 1
    assert plot_paths[0].exists()
    assert "| bb_row_2 | 2 | 2 |" in md_path.read_text()

benchmark_488.run_smoke(parsed_args(benchmark_488, SCRIPT_488, "--smoke"))
benchmark_bb.run_smoke(parsed_args(benchmark_bb, SCRIPT_BB, "--smoke"))

print("test_additional_benchmarks: ok")
