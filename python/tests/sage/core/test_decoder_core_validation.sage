"""Check shared validation and NumPy sampling."""

import builtins

import _bootstrap  # noqa: F401

from sage.all import GF

from decoder_core.validation import sample_bsc_error


original_import = builtins.__import__


def reject_ldpc_import(name, *args, **kwargs):
    if name == "ldpc" or name.startswith("ldpc."):
        raise AssertionError("sample_bsc_error must not import ldpc")
    return original_import(name, *args, **kwargs)


builtins.__import__ = reject_ldpc_import
try:
    all_zero = sample_bsc_error(8, 0.0)
    all_one = sample_bsc_error(5, 1.0)
finally:
    builtins.__import__ = original_import

assert len(all_zero) == 8
assert all_zero.base_ring() == GF(2)
assert list(all_zero) == [0] * 8
assert len(all_one) == 5
assert all_one.base_ring() == GF(2)
assert list(all_one) == [1] * 5

print("test_decoder_core_validation: ok")
