"""ldpc package adapters used by decoder benchmarks."""

from sage.all import GF, vector

from .arrays import to_numpy_uint8


def require_bposd():
    """Import and return the BPOSD decoder class with a clear error.

    Args:
        None.

    Returns:
        object: Constructed object or imported class described by the function summary.
    """
    # Reject invalid inputs early so downstream algebra sees canonical data.
    try:
        from ldpc import BpOsdDecoder
    except ImportError as exc:
        raise ValueError("ldpc is required for BPOSD benchmarking.") from exc
    return BpOsdDecoder


def bposd_decoder(
    hx_numpy,
    p,
    max_iter,
    osd_order,
    *,
    bp_method="minimum_sum",
    ms_scaling_factor=0.625,
):
    """Construct one ``ldpc.BpOsdDecoder`` with Python-native numeric args.

    Args:
        hx_numpy: NumPy representation of the X-check matrix.
        p: Physical error probability.
        max_iter: Maximum number of BP iterations.
        osd_order: Ordered-statistics decoding order for BPOSD.
        bp_method: Belief-propagation update rule requested from the LDPC package.
        ms_scaling_factor: Minimum-sum BP scaling factor.

    Returns:
        object: Construct one ``ldpc.BpOsdDecoder`` with Python-native numeric args.

    Example:
        decoder = bposd_decoder(hx_numpy, 0.01, max_iter=50, osd_order=2)
    """
    # Keep syndrome normalization separate from the matching solve.
    BpOsdDecoder = require_bposd()
    return BpOsdDecoder(
        hx_numpy,
        error_rate=float(p),
        max_iter=int(max_iter),
        bp_method=bp_method,
        ms_scaling_factor=float(ms_scaling_factor),
        osd_method="OSD_CS",
        osd_order=int(osd_order),
    )


def decode_with_bposd(decoder, syndrome):
    """Decode a syndrome with BPOSD and return a Sage GF(2) vector.

    Args:
        decoder: Decoder instance that provides finite maps, target metadata, or LDPC state.
        syndrome: GF(2) syndrome vector to decode or classify.

    Returns:
        vector: GF(2) correction vector in the relevant finite layout.

    """
    # Normalize inputs before running the decoding pipeline.
    decoded = decoder.decode(to_numpy_uint8(syndrome))
    return vector(GF(2), [int(value) for value in decoded])
