"""Shared validation and result writer for the two paper color-code examples."""

from pathlib import Path
import csv
import io
import json
import os
import tempfile
import time

from sage.all import identity_matrix

from decoder_core.benchmark_metadata import runtime_metadata
from decoupling import (
    coarse_grain_to_superlattice,
    compute_translation_representation,
    dagger_matrix,
    find_anyon_preserving_superlattice_from_generators,
    find_anyon_preserving_superlattice_from_translation_representation,
    solve_decoupling_unitary,
)
from decoupling.chain_maps.decoupling import verify_qca_decoupling
from decoupling.css import check_commutation, split_css_blocks


def _matrix_degree(matrix):
    """Return the maximal Laurent L1 degree among matrix monomials."""
    degree = 0
    for row, column in matrix.nonzero_positions():
        for x_power, y_power in matrix[row, column].dict():
            degree = max(degree, abs(x_power) + abs(y_power))
    return int(degree)


def _matrix_rows(matrix):
    """Return a JSON-safe row representation of a Laurent matrix."""
    return [
        [str(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ]


def _basis_index(basis):
    """Return the absolute determinant of a two-vector basis."""
    first, second = basis
    return abs(first[0] * second[1] - first[1] * second[0])


def _period_data(case, excitation_map):
    """Compute anyon-preserving periods from the paper input."""
    if case.label == "666":
        h_x = excitation_map[:1, :2]
        return find_anyon_preserving_superlattice_from_generators(
            h_x[0, 0],
            h_x[0, 1],
            max_period=case.reference_square_period,
        )
    representation = compute_translation_representation(excitation_map)
    return find_anyon_preserving_superlattice_from_translation_representation(
        representation.tx,
        representation.ty,
        max_period=case.reference_square_period,
    )


def _verify_superlattice(case, periods):
    """Verify membership and index for the paper-specified basis."""
    if int(periods.square_period) != case.reference_square_period:
        raise AssertionError("The computed square period does not match the paper.")
    period_vectors = set(periods.period_vectors)
    if any(vector not in period_vectors for vector in case.superlattice_basis):
        raise AssertionError("A paper basis vector is not an anyon-preserving period.")
    paper_index = _basis_index(case.superlattice_basis)
    computed_index = min(
        abs(u[0] * v[1] - u[1] * v[0])
        for position, u in enumerate(periods.period_vectors)
        for v in periods.period_vectors[position + 1 :]
        if u[0] * v[1] - u[1] * v[0]
    )
    if paper_index != computed_index:
        raise AssertionError("The paper basis does not span the computed superlattice.")
    return paper_index


def _verify_chain_maps(result):
    """Verify both forward equations, both inverse equations, and invertibility."""
    num_x_checks = result.h_x_tilde.nrows()
    num_qubits = result.h_x_tilde.ncols()
    h_z, h_x = split_css_blocks(
        result.input_excitation_map,
        num_x_checks,
        num_qubits,
    )
    psi = result.psi
    inverse = result.psi_inverse
    if psi is None:
        raise AssertionError("Color-code reproduction requires explicit forward maps.")
    if result.h_x_tilde * psi.psi_1 != psi.psi_0 * h_x:
        raise AssertionError("The forward degree-one equation failed.")
    if result.h_z_tilde_dagger * psi.psi_2 != psi.psi_1 * dagger_matrix(h_z):
        raise AssertionError("The forward degree-two equation failed.")
    if h_x * inverse.psi_1_inverse != inverse.psi_0_inverse * result.h_x_tilde:
        raise AssertionError("The inverse degree-one equation failed.")
    if (
        dagger_matrix(h_z) * inverse.psi_2_inverse
        != inverse.psi_1_inverse * result.h_z_tilde_dagger
    ):
        raise AssertionError("The inverse degree-two equation failed.")
    for forward, backward in (
        (psi.psi_2, inverse.psi_2_inverse),
        (psi.psi_1, inverse.psi_1_inverse),
        (psi.psi_0, inverse.psi_0_inverse),
    ):
        identity = identity_matrix(forward.base_ring(), forward.nrows())
        if forward * backward != identity or backward * forward != identity:
            raise AssertionError("A chain-isomorphism invertibility check failed.")


def _atomic_write(path, text):
    """Atomically replace one UTF-8 result file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        dir=path.parent,
        text=True,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _csv_text(metadata, record):
    """Render one complete result row, including runtime metadata, as CSV."""
    columns = [
        "code",
        "input_excitation_map",
        "actual_square_period",
        "superlattice_basis",
        "actual_q",
        "actual_p_x",
        "actual_p_z",
        "actual_t",
        "actual_psi_1_inverse_degree",
        "decoupling_runtime_seconds",
        "actual_runtime_seconds",
        "verification_status",
        "adapter",
        "suite",
        "source_revision",
        "source_tree_dirty",
        "sage_version",
        "python_version",
        "platform",
        "paper_runtime_included",
    ]
    stream = io.StringIO()
    writer = csv.DictWriter(
        stream,
        fieldnames=columns,
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    flat = dict(record)
    flat.update(metadata)
    flat["input_excitation_map"] = json.dumps(
        record["input_excitation_map"],
        separators=(",", ":"),
    )
    flat["superlattice_basis"] = json.dumps(
        record["superlattice_basis"],
        separators=(",", ":"),
    )
    writer.writerow(flat)
    return stream.getvalue()


def _markdown_text(metadata, record):
    """Render one color-code result summary as Markdown."""
    return "\n".join(
        [
            f"# {record['code']} Color-Code Decoupling Reproduction",
            "",
            f"- Input excitation map: `{record['input_excitation_map']}`",
            f"- Verification: `{record['verification_status']}`",
            f"- Square period L: `{record['actual_square_period']}`",
            f"- Paper superlattice basis: `{record['superlattice_basis']}`",
            f"- Standard sectors (q, p_X, p_Z, t): "
            f"`({record['actual_q']}, {record['actual_p_x']}, "
            f"{record['actual_p_z']}, {record['actual_t']})`",
            f"- Actual deg(psi_1^-1): `{record['actual_psi_1_inverse_degree']}`",
            f"- Actual decoupling runtime: `{record['decoupling_runtime_seconds']:.6f}` s",
            f"- Actual total runtime: `{record['actual_runtime_seconds']:.6f}` s",
            f"- Sage: `{metadata['sage_version']}`",
            f"- Python: `{metadata['python_version']}`",
            f"- Platform: `{metadata['platform']}`",
            f"- Adapter: `{metadata['adapter']}`",
            f"- Suite: `{metadata['suite']}`",
            f"- Source revision: `{metadata['source_revision']}`",
            f"- Source tree dirty: `{metadata['source_tree_dirty']}`",
            f"- Paper runtime included: `{metadata['paper_runtime_included']}`",
            "",
            "The chain isomorphisms are not unique. A valid Gröbner basis and "
            "basis-selection path can produce matrices that differ entry by entry "
            "from those printed in the paper while satisfying every chain, inverse, "
            "and symplectic identity.",
            "",
        ]
    )


def run_color_code(case, source_root):
    """Validate a color code and publish canonical JSON, CSV, and Markdown."""
    source_root = Path(source_root)
    started = time.perf_counter()
    excitation_map = case.excitation_map()
    original_num_x_checks = excitation_map.nrows() // 2
    original_num_qubits = excitation_map.ncols() // 2
    if not check_commutation(excitation_map, original_num_qubits):
        raise AssertionError("The paper input violates CSS commutation.")

    periods = _period_data(case, excitation_map)
    index = _verify_superlattice(case, periods)
    coarse = coarse_grain_to_superlattice(
        excitation_map,
        *case.superlattice_basis,
    )
    actual_q = coarse.ncols() // 2
    if actual_q != case.reference_q or actual_q != original_num_qubits * index:
        raise AssertionError("The coarse qubit count does not match the paper.")

    decoupling_started = time.perf_counter()
    result = solve_decoupling_unitary(
        coarse,
        num_x_checks=coarse.nrows() // 2,
        num_qubits=actual_q,
        compute_psi=True,
    )
    decoupling_seconds = time.perf_counter() - decoupling_started
    _verify_chain_maps(result)
    if not verify_qca_decoupling(result):
        raise AssertionError("The symplectic QCA decoupling certificate failed.")

    actual_p_x = int(result.diagnostics["product_x_rank"])
    actual_p_z = int(result.diagnostics["product_z_rank"])
    actual_t_x = result.h_x_tilde.nrows() - actual_p_x
    actual_t_z = result.h_z_tilde_dagger.ncols() - actual_p_z
    if actual_t_x != actual_t_z:
        raise AssertionError("The standard complexes have unequal toric counts.")
    actual_t = int(actual_t_x)
    actual = (actual_q, actual_p_x, actual_p_z, actual_t)
    expected = (
        case.reference_q,
        case.reference_p_x,
        case.reference_p_z,
        case.reference_t,
    )
    if actual != expected:
        raise AssertionError(f"Computed (q,p_X,p_Z,t)={actual}, expected {expected}.")

    record = {
        "code": case.label,
        "input_excitation_map": _matrix_rows(excitation_map),
        "coarse_excitation_map": _matrix_rows(coarse),
        "actual_square_period": int(periods.square_period),
        "superlattice_basis": [list(vector) for vector in case.superlattice_basis],
        "actual_q": actual_q,
        "actual_p_x": actual_p_x,
        "actual_p_z": actual_p_z,
        "actual_t": actual_t,
        "h_x_tilde": _matrix_rows(result.h_x_tilde),
        "h_z_tilde_dagger": _matrix_rows(result.h_z_tilde_dagger),
        "psi_2": _matrix_rows(result.psi.psi_2),
        "psi_1": _matrix_rows(result.psi.psi_1),
        "psi_0": _matrix_rows(result.psi.psi_0),
        "psi_2_inverse": _matrix_rows(result.psi_inverse.psi_2_inverse),
        "psi_1_inverse": _matrix_rows(result.psi_inverse.psi_1_inverse),
        "psi_0_inverse": _matrix_rows(result.psi_inverse.psi_0_inverse),
        "actual_psi_1_inverse_degree": _matrix_degree(
            result.psi_inverse.psi_1_inverse
        ),
        "decoupling_runtime_seconds": decoupling_seconds,
        "actual_runtime_seconds": time.perf_counter() - started,
        "verification_status": "passed",
    }
    metadata = runtime_metadata(
        adapter="decoupling",
        suite=f"paper_color_code_{case.label}",
        source_root=source_root,
    )
    metadata.pop("source_root", None)
    metadata["paper_runtime_included"] = False

    result_directory = source_root / "results" / "decoupling"
    stem = f"color_code_{case.label}"
    _atomic_write(
        result_directory / f"{stem}.json",
        json.dumps(
            {"metadata": metadata, "result": record},
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    _atomic_write(result_directory / f"{stem}.csv", _csv_text(metadata, record))
    _atomic_write(result_directory / f"{stem}.md", _markdown_text(metadata, record))
    return record
