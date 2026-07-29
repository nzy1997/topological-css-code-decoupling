"""Shared benchmark statistics containers and binomial error bars."""

from dataclasses import dataclass
from math import sqrt


@dataclass(frozen=True)
class RateErrorBar:
    """Error-bar data for one empirical binomial rate."""

    failures: int
    shots: int
    rate: float
    display_rate: float
    zero_rate: bool
    interval_low: float
    interval_high: float
    lower_error: float
    upper_error: float


@dataclass
class DecoderStats:
    """Aggregate logical-error and runtime statistics for one decoder."""

    decoder: str
    shots: int = 0
    logical_failures: int = 0
    decode_failures: int = 0
    elapsed_seconds: float = 0.0

    @property
    def logical_error_rate(self):
        """Return the empirical logical error rate.

        Args:
            None.

        Returns:
            object: Statistical estimate or interval data for benchmark reporting.
        """
        # Keep syndrome normalization separate from the matching solve.
        if self.shots == 0:
            return 0.0
        return self.logical_failures / self.shots

    @property
    def standard_error(self):
        """Return the binomial standard error for the logical error rate.

        Args:
            None.

        Returns:
            float: The binomial standard error for the logical error rate.
        """
        # Keep syndrome normalization separate from the matching solve.
        return binomial_standard_error(self.logical_failures, self.shots)

    def record(self, logical_failure, decode_failure, elapsed_seconds):
        """Record one decoding shot.

        Args:
            logical_failure: Whether the residual is a nontrivial logical error.
            decode_failure: Whether a decoder failed to reproduce the syndrome.
            elapsed_seconds: Wall-clock decode time to accumulate in benchmark stats.

        Returns:
            None: This function mutates local state or performs validation only.
        """
        # Keep syndrome normalization separate from the matching solve.
        self.shots += 1
        self.elapsed_seconds += float(elapsed_seconds)
        if logical_failure:
            self.logical_failures += 1
        if decode_failure:
            self.decode_failures += 1


def binomial_standard_error(failures, shots):
    """Return the binomial standard error for ``failures / shots``.

    Args:
        failures: Number of logical failures observed.
        shots: Number of Monte Carlo shots sampled.

    Returns:
        vector: The binomial standard error for ``failures / shots``.
    """
    # Keep syndrome normalization separate from the matching solve.
    shots = int(shots)
    if shots <= 0:
        return 0.0
    failures = int(failures)
    _validate_binomial_counts(failures, shots)
    rate = failures / shots
    return sqrt(rate * (1.0 - rate) / shots)


def binomial_wilson_interval(failures, shots, z=1.0):
    """Return a Wilson binomial interval for one failure count.

    Args:
        failures: Number of logical failures observed.
        shots: Number of Monte Carlo shots sampled.
        z: Normal-approximation z value for confidence intervals.

    Returns:
        object: Statistical estimate or interval data for benchmark reporting.
    """
    # Keep syndrome normalization separate from the matching solve.
    shots = int(shots)
    if shots <= 0:
        return 0.0, 0.0
    failures = int(failures)
    _validate_binomial_counts(failures, shots)
    phat = failures / shots
    z = float(z)
    denom = 1.0 + z * z / shots
    center = (phat + z * z / (2.0 * shots)) / denom
    half_width = (
        z
        * (phat * (1.0 - phat) / shots + z * z / (4.0 * shots * shots)) ** 0.5
        / denom
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def binomial_rate_error_bar(failures, shots, *, z=1.0, display_floor=None):
    """Return log-plot-safe rate and asymmetric Wilson error bars.

    Args:
        failures: Number of logical failures observed.
        shots: Number of Monte Carlo shots sampled.
        z: Normal-approximation z value for confidence intervals.
        display_floor: Small positive plotting floor for zero-rate estimates.

    Returns:
        object: Statistical estimate or interval data for benchmark reporting.
    """
    # Keep syndrome normalization separate from the matching solve.
    shots = int(shots)
    if shots <= 0:
        return RateErrorBar(
            failures=0,
            shots=0,
            rate=0.0,
            display_rate=0.0,
            zero_rate=True,
            interval_low=0.0,
            interval_high=0.0,
            lower_error=0.0,
            upper_error=0.0,
        )
    failures = int(failures)
    _validate_binomial_counts(failures, shots)
    rate = failures / shots
    interval_low, interval_high = binomial_wilson_interval(failures, shots, z=z)
    floor = 0.5 / shots if display_floor is None else float(display_floor)
    if rate == 0.0:
        # Log-scale plots cannot display an exact zero.  Display zero-failure
        # points at half a shot while keeping the statistical rate itself zero.
        display_rate = floor
        lower_error = 0.0
        upper_error = max(0.0, interval_high - display_rate)
        zero_rate = True
    else:
        display_rate = rate
        lower = max(interval_low, floor)
        lower_error = max(0.0, rate - lower)
        upper_error = max(0.0, interval_high - rate)
        zero_rate = False
    return RateErrorBar(
        failures=failures,
        shots=shots,
        rate=rate,
        display_rate=display_rate,
        zero_rate=zero_rate,
        interval_low=interval_low,
        interval_high=interval_high,
        lower_error=lower_error,
        upper_error=upper_error,
    )


def _validate_binomial_counts(failures, shots):
    """Validate binomial counts.

    Args:
        failures: Number of logical failures observed.
        shots: Number of Monte Carlo shots sampled.

    Returns:
        None: This function mutates local state or performs validation only.
    """
    # Reject invalid inputs early so downstream algebra sees canonical data.
    if failures < 0 or failures > shots:
        raise ValueError("failures must be between 0 and shots.")
