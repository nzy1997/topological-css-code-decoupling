"""Check the 6.6.6 unitary-vs-BPOSD benchmark configuration."""

from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
import tempfile

import _bootstrap

from decoder_core import (
    LogicalFailureClassifier,
    bposd_decoder,
    decode_with_bposd,
    to_numpy_uint8,
)
from isomorphism import x, y
from unitary_decouple_decoder import DecouplingDecoder
from unitary_decouple_decoder.benchmarking import build_hz_dagger_matrix


REPO_ROOT = _bootstrap.REPO_ROOT
SCRIPT_PATH = REPO_ROOT / "scripts" / "unitary" / "benchmark_color_666_unitary_vs_bposd.sage"
benchmark = SourceFileLoader("benchmark_color_666_unitary_vs_bposd", str(SCRIPT_PATH)).load_module()

old_argv = sys.argv
try:
    sys.argv = [str(SCRIPT_PATH)]
    default_args = benchmark.parse_args()
finally:
    sys.argv = old_argv

assert benchmark.parse_distances(default_args.distances) == [4, 6, 8]
assert benchmark.parse_probabilities(default_args.p_list) == [0.01, 0.02, 0.03, 0.05, 0.08, 0.1]
assert default_args.shots == 100000
assert default_args.progress_interval == 10000
assert default_args.output == "results/unitary_decouple_decoder/color_666_unitary_vs_bposd.csv"
assert default_args.plot_output == "results/unitary_decouple_decoder/color_666_unitary_vs_bposd.png"

decoder = DecouplingDecoder.from_bb(1 + x + x * y, 1 + y + x * y, 2, 2)
hz_dagger = build_hz_dagger_matrix(decoder)
classifier = LogicalFailureClassifier(decoder.hx_source_finite, hz_dagger)

assert hz_dagger.nrows() == decoder.correction_size
assert decoder.hx_source_finite.nrows() == decoder.syndrome_size
assert not classifier.is_logical_failure(decoder.zero_correction())

bposd = bposd_decoder(
    to_numpy_uint8(decoder.hx_source_finite),
    0.05,
    int(20),
    int(1),
)
decoded = decode_with_bposd(bposd, decoder.zero_syndrome())
assert len(decoded) == decoder.correction_size
assert decoder.syndrome(decoded).is_zero()

rows = [
    {
        "distance": 2,
        "n": decoder.correction_size,
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
        "distance": 2,
        "n": decoder.correction_size,
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
    output = Path(tmpdir) / "benchmark.csv"
    plot = Path(tmpdir) / "benchmark.png"
    md_path, plot_path = benchmark.write_outputs(rows, output, "Color 666 Test", plot)
    assert output.exists()
    assert md_path.exists()
    assert plot_path == plot
    assert plot.exists()
    assert "| 2 | 24 | 0.05 | unitary_decouple | 5 | 1 | 0 | 0.2000 |" in md_path.read_text()

try:
    sys.argv = [
        str(SCRIPT_PATH),
        "--distances",
        "2",
        "--p-list",
        "0.05",
        "--shots",
        "3",
        "--progress-interval",
        "3",
    ]
    args = benchmark.parse_args()
finally:
    sys.argv = old_argv

assert benchmark.parse_distances(args.distances) == [2]
assert benchmark.parse_probabilities(args.p_list) == [0.05]
assert args.shots == 3

try:
    sys.argv = [str(SCRIPT_PATH), "--smoke"]
    smoke_args = benchmark.parse_args()
finally:
    sys.argv = old_argv
benchmark.run_smoke(smoke_args)

print("test_color_666_benchmark: ok")
