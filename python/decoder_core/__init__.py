"""Shared finite decoder primitives.

Exports are loaded on demand so pure-Python helpers such as statistics and
plotting can be imported without requiring a SageMath runtime.
"""

from importlib import import_module


_EXPORT_MODULES = {
    "BBRow": "bb_cases",
    "BB_ROWS": "bb_cases",
    "BB_ROWS_LTE_63": "bb_cases",
    "DecoderStats": "stats",
    "LogicalFailureClassifier": "logical",
    "RateErrorBar": "stats",
    "ToricMatchingDecoder": "toric",
    "TorusShape": "finite_torus",
    "binomial_rate_error_bar": "stats",
    "binomial_standard_error": "stats",
    "binomial_wilson_interval": "stats",
    "bposd_decoder": "ldpc",
    "classify_attempt": "logical",
    "classify_verified_attempt": "logical",
    "decode_with_bposd": "ldpc",
    "numpy_rng": "sampling",
    "plot_logical_error_rates": "plotting",
    "plot_rate_points": "plotting",
    "require_bposd": "ldpc",
    "runtime_metadata": "benchmark_metadata",
    "sample_data_error": "sampling",
    "source_revision": "benchmark_metadata",
    "to_numpy_uint8": "arrays",
}

__all__ = sorted(_EXPORT_MODULES)


def __getattr__(name):
    """Load exported helpers only when requested.

    Args:
        name: Human-readable field name used in validation errors.

    Returns:
        object: Load exported helpers only when requested.
    """
    # Keep syndrome normalization separate from the matching solve.
    module_name = _EXPORT_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, name)
    globals()[name] = value
    return value
