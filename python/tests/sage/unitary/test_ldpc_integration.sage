"""Check ldpc integration points for sampling and BP-OSD comparison."""

import _bootstrap  # noqa: F401

from decoder_core.validation import matrix_to_numpy_uint8, sample_bsc_error
from isomorphism import x, y
from unitary_decouple_decoder import DecouplingDecoder


decoder = DecouplingDecoder.from_bb(1 + x + x * y, 1 + y + x * y, 2, 2)
sample = sample_bsc_error(decoder.correction_size, 0.0)
syndrome = decoder.syndrome(sample)

from ldpc import BpOsdDecoder

bp_osd = BpOsdDecoder(
    matrix_to_numpy_uint8(decoder.hx_source_finite),
    error_rate=float(0.05),
    max_iter=int(5),
    osd_order=int(0),
)
decoded = bp_osd.decode(matrix_to_numpy_uint8(syndrome))

assert len(sample) == decoder.correction_size
assert len(decoded) == decoder.correction_size

print("test_ldpc_integration: ok")
