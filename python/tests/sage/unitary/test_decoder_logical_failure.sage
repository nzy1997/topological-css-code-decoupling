"""A syndrome-consistent logical residual must not be reported as success."""

import _bootstrap  # noqa: F401

from sage.all import GF, vector

from decoupling import x, y
from unitary_decouple_based_decoder import UnitaryDecoupleBasedDecoder
from unitary_decouple_based_decoder.benchmarking import build_h_z_dagger_finite


decoder = UnitaryDecoupleBasedDecoder.from_bb(
    1 + x + x * y,
    1 + y + x * y,
    2,
    2,
)
error = vector(GF(2), decoder.correction_size)
error[0] = 1

assert decoder._h_z_dagger_input_finite is None
decoder.decode(decoder.syndrome(error), verify=True)
assert decoder._h_z_dagger_input_finite is None

assert hasattr(decoder, "decode_error"), (
    "The decoder must classify a supplied error independently of sampling."
)
record = decoder.decode_error(error, verify=True)
stabilizer_space = build_h_z_dagger_finite(decoder).column_space()

assert decoder._h_z_dagger_input_finite is not None
assert decoder.syndrome(record["correction"]) == record["syndrome"]
assert record["residual"] == error + record["correction"]
assert record["residual_syndrome"].is_zero()
assert record["residual"] not in stabilizer_space
assert record["logical_failure"]
assert not record["decode_failure"]
assert not record["success"]

sample_error = decoder.sample_error
decoder.sample_error = lambda _p: error
try:
    sampled_record = decoder.decode_sample(0.25, verify=True)
finally:
    decoder.sample_error = sample_error

assert sampled_record["residual"] == record["residual"]
assert sampled_record["logical_failure"]
assert not sampled_record["success"]

decode = decoder.decode
decoder.decode = lambda _syndrome, *, verify=False: decoder.zero_correction()
try:
    failed_record = decoder.decode_error(error, verify=False)
finally:
    decoder.decode = decode

assert failed_record["residual_syndrome"] == failed_record["syndrome"]
assert failed_record["logical_failure"]
assert failed_record["decode_failure"]
assert not failed_record["success"]

print("test_decoder_logical_failure: ok")
