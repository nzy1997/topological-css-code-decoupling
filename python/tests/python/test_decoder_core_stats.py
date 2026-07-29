"""Check pure-Python benchmark statistics helpers."""

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
repo_root = str(REPO_ROOT)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from decoder_core.stats import (  # noqa: E402
    DecoderStats,
    binomial_rate_error_bar,
    binomial_standard_error,
    binomial_wilson_interval,
)


stats = DecoderStats("example")
for logical_failure in (True, False, False, False):
    stats.record(logical_failure, False, 0.25)

assert stats.shots == 4
assert stats.logical_failures == 1
assert abs(stats.logical_error_rate - 0.25) < 1e-12
assert abs(stats.standard_error - binomial_standard_error(1, 4)) < 1e-12

zero_bar = binomial_rate_error_bar(0, 100)
assert zero_bar.rate == 0.0
assert zero_bar.display_rate == 0.005
assert zero_bar.zero_rate
assert zero_bar.lower_error == 0.0
assert zero_bar.upper_error > 0.0

nonzero_bar = binomial_rate_error_bar(10, 100)
low, high = binomial_wilson_interval(10, 100)
assert nonzero_bar.rate == 0.1
assert nonzero_bar.display_rate == 0.1
assert not nonzero_bar.zero_rate
assert abs(nonzero_bar.interval_low - low) < 1e-12
assert abs(nonzero_bar.interval_high - high) < 1e-12
assert abs(nonzero_bar.lower_error - (0.1 - max(low, 0.005))) < 1e-12
assert abs(nonzero_bar.upper_error - (high - 0.1)) < 1e-12

print("test_decoder_core_stats: ok")
