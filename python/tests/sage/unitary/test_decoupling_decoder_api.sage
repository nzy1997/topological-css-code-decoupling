"""Check the supplement-oriented public decoupling decoder API."""

import _bootstrap  # noqa: F401

from dataclasses import fields

import unitary_decouple_decoder
from isomorphism import DecouplingMaps, DecouplingResult, InverseDecouplingMaps, x, y
from unitary_decouple_decoder import DecouplingDecoder
from unitary_decouple_decoder.standard_code import StandardCodeDecoder


assert [item.name for item in fields(DecouplingMaps)] == ["phi2", "phi1", "phi0"]
assert [item.name for item in fields(InverseDecouplingMaps)] == [
    "phi2_inverse",
    "phi1_inverse",
    "phi0_inverse",
]
assert [item.name for item in fields(DecouplingResult)] == [
    "source_coarse_matrix",
    "hx_standard",
    "hz_standard",
    "inverse_maps",
    "maps",
    "diagnostics",
]
assert unitary_decouple_decoder.__all__ == ["DecouplingDecoder"]
assert not hasattr(unitary_decouple_decoder, "Unitary" + "DecoupleDecoder")
assert not hasattr(unitary_decouple_decoder, "CoarseUnitary" + "DecoupleDecoder")
assert not hasattr(unitary_decouple_decoder, "StandardCodeDecoder")

bb_decoder = DecouplingDecoder.from_bb(
    1 + x + x * y,
    1 + y + x * y,
    2,
    2,
)
assert bb_decoder.decoupling_result.maps is None
assert bb_decoder.inverse_maps is bb_decoder.decoupling_result.inverse_maps
assert bb_decoder.hx_standard == bb_decoder.decoupling_result.hx_standard
assert not hasattr(bb_decoder.inverse_maps, "hx_standard")
assert not hasattr(bb_decoder.inverse_maps, "hz_standard")
assert not hasattr(bb_decoder.inverse_maps, "diagnostics")
assert isinstance(bb_decoder.standard_decoder, StandardCodeDecoder)
assert bb_decoder.check_chain_relation()
assert (
    bb_decoder.hx_source_finite * bb_decoder.phi1_inverse_finite
    == bb_decoder.phi0_inverse_finite * bb_decoder.hx_standard_finite
)
assert (
    bb_decoder.phi0_finite * bb_decoder.phi0_inverse_finite
).is_one()

for removed_name in (
    "maps",
    "q" + "_finite",
    "sig" + "ma_finite",
    "sig" + "ma_inverse_finite",
    "hx_" + "target",
    "hx_" + "target_finite",
):
    assert not hasattr(bb_decoder, removed_name)

coarse_decoder = DecouplingDecoder.from_coarse_css(
    bb_decoder.source_epsilon,
    num_x_checks=bb_decoder.num_x_checks,
    num_qubits=bb_decoder.num_qubits,
    Lx=2,
    Ly=2,
)
assert coarse_decoder.source_epsilon == bb_decoder.source_epsilon
assert coarse_decoder.num_x_checks == bb_decoder.num_x_checks
assert coarse_decoder.num_qubits == bb_decoder.num_qubits
assert coarse_decoder.check_chain_relation()

print("test_decoupling_decoder_api: ok")
