"""Check standard product-state plus toric-sector decoding."""

import _bootstrap  # noqa: F401

from sage.all import Matrix

from decoder_core import TorusShape
from isomorphism import R, x, y
from unitary_decouple_decoder.standard_code import StandardCodeDecoder

shape = TorusShape(3, 3)
hx_standard = Matrix(R, [[1, 0, 0], [0, 1 + x, 1 + y]])
decoder = StandardCodeDecoder(hx_standard, shape)

syndrome = decoder.zero_syndrome()
syndrome[shape.index(0, 0)] = 1
syndrome[shape.size + shape.index(0, 0)] = 1
syndrome[shape.size + shape.index(1, 0)] = 1

correction = decoder.decode(syndrome)
assert decoder.hx_standard_finite * correction == syndrome
assert correction[shape.index(0, 0)] == 1

bad = decoder.zero_syndrome()
bad[shape.size + shape.index(0, 0)] = 1
try:
    decoder.decode(bad)
except ValueError as exc:
    assert "even parity" in str(exc)
else:
    raise AssertionError("Expected odd toric syndrome to fail.")

print("test_standard_code: ok")
