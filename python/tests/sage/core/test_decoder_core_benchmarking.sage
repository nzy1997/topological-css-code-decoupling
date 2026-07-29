"""Check shared benchmark statistics and plotting helpers."""

import _bootstrap  # noqa: F401

from math import sqrt
from pathlib import Path
import tempfile

from decoder_core import (
    DecoderStats,
    binomial_rate_error_bar,
    binomial_wilson_interval,
    plot_logical_error_rates,
    plot_rate_points,
)


stats = DecoderStats("example")
for logical_failure in (True, False, False, False):
    stats.record(logical_failure, False, 0.25)

assert stats.shots == 4
assert stats.logical_failures == 1
assert stats.decode_failures == 0
assert abs(stats.logical_error_rate - 0.25) < 1e-12
assert abs(stats.standard_error - sqrt(0.25 * 0.75 / 4)) < 1e-12
assert abs(stats.elapsed_seconds - 1.0) < 1e-12

rows = [
    {
        "distance": 4,
        "p": 0.01,
        "decoder": "unitary_decouple",
        "shots": 100,
        "logical_failures": 0,
        "logical_error_rate": 0.0,
        "standard_error": 0.0,
    },
    {
        "distance": 4,
        "p": 0.02,
        "decoder": "unitary_decouple",
        "shots": 100,
        "logical_failures": 10,
        "logical_error_rate": 0.1,
        "standard_error": 0.03,
    },
    {
        "distance": 4,
        "p": 0.01,
        "decoder": "bposd",
        "shots": 100,
        "logical_failures": 0,
        "logical_error_rate": 0.0,
        "standard_error": 0.0,
    },
]

probabilities, rates, zero_rate, yerr = plot_rate_points(rows[:2])
assert probabilities == [0.01, 0.02]
assert rates == [0.005, 0.1]
assert zero_rate == [True, False]
assert yerr[0][0] == 0.0
assert yerr[1][0] > 0.0
low, high = binomial_wilson_interval(10, 100)
error_bar = binomial_rate_error_bar(10, 100)
assert 0.0 < low < 0.1 < high < 1.0
assert abs(error_bar.interval_low - low) < 1e-12
assert abs(error_bar.interval_high - high) < 1e-12
assert abs(yerr[0][1] - (0.1 - low)) < 1e-12
assert abs(yerr[1][1] - (high - 0.1)) < 1e-12

with tempfile.TemporaryDirectory() as tmpdir:
    plot_path = Path(tmpdir) / "plot.png"
    written = plot_logical_error_rates(rows, plot_path, "Benchmark Core Test")
    assert written == plot_path
    assert plot_path.exists()

print("test_decoder_core_benchmarking: ok")
