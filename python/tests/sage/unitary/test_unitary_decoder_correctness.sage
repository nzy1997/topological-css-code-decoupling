"""Check the full unitary-decouple decode pipeline on a BB code."""

import _bootstrap  # noqa: F401

from sage.all import GF, vector

from isomorphism import x, y
from unitary_decouple_decoder import DecouplingDecoder


decoder = DecouplingDecoder.from_bb(1 + x + x * y, 1 + y + x * y, 2, 2)
assert decoder.check_chain_relation()

error = vector(GF(2), decoder.correction_size)
error[0] = 1
error[decoder.shape.size + decoder.shape.index(1, 0)] = 1

syndrome = decoder.syndrome(error)
correction = decoder.decode(syndrome, verify=True)

assert len(correction) == decoder.correction_size
assert decoder.syndrome(correction) == syndrome

print("test_unitary_decoder_correctness: ok")
