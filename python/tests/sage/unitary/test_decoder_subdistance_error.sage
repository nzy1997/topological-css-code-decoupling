"""A correctable weight-one error must leave a stabilizer residual."""

import _bootstrap  # noqa: F401

from sage.all import GF, vector

from decoupling import x, y
from unitary_decouple_based_decoder import UnitaryDecoupleBasedDecoder
from unitary_decouple_based_decoder.benchmarking import build_h_z_dagger_finite


decoder = UnitaryDecoupleBasedDecoder.from_bb(
    1 + x + x**-1 * y,
    1 + y + x * y,
    4,
    4,
    max_period=7,
)
error = vector(GF(2), decoder.correction_size)
error[0] = 1

assert hasattr(decoder, "decode_error"), (
    "The decoder must classify a supplied error independently of sampling."
)
record = decoder.decode_error(error, verify=True)
stabilizer_space = build_h_z_dagger_finite(decoder).column_space()

assert record["residual"] == error + record["correction"]
assert record["residual_syndrome"].is_zero()
assert record["residual"] in stabilizer_space
assert not record["logical_failure"]
assert not record["decode_failure"]
assert record["success"]

print("test_decoder_subdistance_error: ok")
