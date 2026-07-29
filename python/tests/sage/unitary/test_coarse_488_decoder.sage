"""Check unitary-decouple decoding from a pre-coarse-grained 4.8.8 color code."""

import _bootstrap  # noqa: F401

from sage.all import GF, Matrix, block_matrix, vector, zero_matrix

from isomorphism import R, oblique_coarse_grain, x, y
from unitary_decouple_decoder import DecouplingDecoder
from unitary_decouple_decoder.benchmarking import build_hz_dagger_matrix


check_488 = Matrix(
    R,
    [
        [1 + y, 1 + y, 1 + x, 1 + x],
        [x * y, y, x * y, x],
    ],
)
zero = zero_matrix(R, 2, 4)
epsilon_488 = block_matrix(R, [[check_488, zero], [zero, check_488]])
coarse_488 = oblique_coarse_grain(epsilon_488, (2, 0), (1, 1))

decoder = DecouplingDecoder.from_coarse_css(
    coarse_488,
    num_x_checks=4,
    num_qubits=8,
    Lx=2,
    Ly=2,
)
assert decoder.check_chain_relation()
assert decoder.num_x_checks == 4
assert decoder.num_qubits == 8

hz_dagger = build_hz_dagger_matrix(decoder)
assert hz_dagger.nrows() == decoder.correction_size

error = vector(GF(2), decoder.correction_size)
error[0] = 1
error[decoder.shape.size + decoder.shape.index(1, 0)] = 1

syndrome = decoder.syndrome(error)
correction = decoder.decode(syndrome, verify=True)

assert len(correction) == decoder.correction_size
assert decoder.syndrome(correction) == syndrome
assert decoder.decode_sample(0.0)["success"]

print("test_coarse_488_decoder: ok")
