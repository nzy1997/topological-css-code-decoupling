"""Check shared toric decoder backend."""

import _bootstrap  # noqa: F401

from decoder_core import ToricMatchingDecoder, TorusShape

shape = TorusShape(3, 3)
decoder = ToricMatchingDecoder(shape, num_stacks=1)
syndrome = decoder.zero_syndrome()
syndrome[decoder.vertex_index(0, 0, 0)] = 1
syndrome[decoder.vertex_index(0, 1, 0)] = 1
correction = decoder.decode(syndrome)

assert decoder.boundary(correction) == syndrome
assert correction[decoder.edge_index(0, "h", 0, 0)] == 1

print("test_decoder_core_toric: ok")
