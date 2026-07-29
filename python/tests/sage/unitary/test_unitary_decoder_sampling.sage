"""Check unitary decoder sampling and validation helpers."""

import _bootstrap  # noqa: F401

from isomorphism import x, y
from unitary_decouple_decoder import DecouplingDecoder


decoder = DecouplingDecoder.from_bb(1 + x + x * y, 1 + y + x * y, 2, 2)
sample = decoder.sample_error(0.0)
record = decoder.decode_sample(0.0)

assert len(sample) == decoder.correction_size
assert all(int(value) == 0 for value in sample)
assert record["success"]
assert all(int(value) == 0 for value in record["residual"])

print("test_unitary_decoder_sampling: ok")
