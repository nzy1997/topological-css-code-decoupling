"""ldpc package adapters used by decoder benchmarks."""

from sage.all import GF, vector

from .arrays import to_numpy_uint8


def require_bp_osd():
    """Import and return the BP-OSD decoder class with a clear error."""
    try:
        from ldpc import BpOsdDecoder
    except ImportError as exc:
        raise ValueError("ldpc is required for BPOSD benchmarking.") from exc
    return BpOsdDecoder


def bp_osd_decoder(
    h_x_numpy,
    p,
    max_iter,
    osd_order,
    *,
    bp_method="minimum_sum",
    ms_scaling_factor=0.625,
):
    """Construct one ``ldpc.BpOsdDecoder`` with Python-native numeric args.

    Args:
        h_x_numpy: NumPy representation of the X-check matrix.
        p: Physical error probability.
        max_iter: Maximum number of BP iterations.
        osd_order: Ordered-statistics decoding order for BPOSD.
        bp_method: Belief-propagation update rule requested from the LDPC package.
        ms_scaling_factor: Minimum-sum BP scaling factor.

    Returns:
        object: Configured ``ldpc.BpOsdDecoder`` instance.

    Example:
        decoder = bp_osd_decoder(h_x_numpy, 0.01, max_iter=50, osd_order=2)
    """
    BpOsdDecoder = require_bp_osd()
    return BpOsdDecoder(
        h_x_numpy,
        error_rate=float(p),
        max_iter=int(max_iter),
        bp_method=bp_method,
        ms_scaling_factor=float(ms_scaling_factor),
        osd_method="OSD_CS",
        osd_order=int(osd_order),
    )


def decode_with_bp_osd(decoder, syndrome):
    """Decode a syndrome with BPOSD and return a Sage GF(2) vector.

    Args:
        decoder: Configured ``ldpc.BpOsdDecoder``.
        syndrome: GF(2) syndrome vector to decode.

    Returns:
        vector: GF(2) correction vector in the relevant finite layout.

    """
    decoded = decoder.decode(to_numpy_uint8(syndrome))
    return vector(GF(2), [int(value) for value in decoded])
