"""Reproduce the BB rows at or below the benchmark-19 runtime cutoff."""

from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import math
import multiprocessing
import os
import platform
import sys
import tempfile
import time
import traceback

from sage.all import identity_matrix, zero_matrix
from sage.env import SAGE_VERSION


SOURCE_ROOT = Path(__file__).resolve().parents[2]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from decoder_core.benchmark_metadata import runtime_metadata  # noqa: E402
from decoupling import (  # noqa: E402
    build_two_generator_css_excitation_map,
    coarse_grain_to_superlattice,
    dagger_matrix,
    find_anyon_preserving_superlattice_from_generators,
    solve_decoupling_unitary,
)
from decoupling.css import check_commutation, split_css_blocks  # noqa: E402
from decoupling.paper_examples import (  # noqa: E402
    BB_CODE_INSTANCES,
    select_bb_code_instances,
)


RESULT_DIRECTORY = SOURCE_ROOT / "results" / "decoupling"
DEFAULT_JSON = RESULT_DIRECTORY / "bb_codes.json"
DEFAULT_CSV = RESULT_DIRECTORY / "bb_codes.csv"
DEFAULT_MARKDOWN = RESULT_DIRECTORY / "bb_codes.md"
DEFAULT_CHECKPOINT = RESULT_DIRECTORY / ".bb_codes.checkpoint.json"
EXPLICIT_FORWARD_Q_LIMIT = 18
RESULT_SCHEMA_VERSION = 1
RUNTIME_CUTOFF_REFERENCE_ID = "benchmark-19"
EXCLUDED_LONGER_RUNTIME_ROW_IDS = frozenset({"benchmark-29", "table-17"})
DEFAULT_VALIDATION_CASES = tuple(
    case
    for case in BB_CODE_INSTANCES
    if case.identifier not in EXCLUDED_LONGER_RUNTIME_ROW_IDS
)


def compute_implementation_fingerprint():
    """Hash the mathematical implementation used to validate resume records."""
    paths = sorted((SOURCE_ROOT / "decoupling").rglob("*.py"))
    paths.append(
        SOURCE_ROOT / "scripts" / "decoupling" / "reproduce_bb_codes.sage"
    )
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(SOURCE_ROOT)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


RUN_IMPLEMENTATION_FINGERPRINT = compute_implementation_fingerprint()


def validation_environment():
    """Return the execution environment attached to each durable row."""
    return {
        "validation_sage_version": SAGE_VERSION,
        "validation_python_version": sys.version.split()[0],
        "validation_platform": platform.platform(),
    }


def positive_float(value):
    """Parse an optional timeout that must be strictly positive."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("timeout must be finite and positive")
    return parsed


def nonnegative_int(value):
    """Parse a nonnegative retry count."""
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("retry count must be nonnegative")
    return parsed


def parse_args():
    """Parse row selection, checkpointing, and retry options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rows",
        help="Comma-separated stable identifiers such as benchmark-01,table-03.",
    )
    parser.add_argument(
        "--row-timeout-seconds",
        type=positive_float,
        default=None,
        help="Optional positive process-level timeout for each row.",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Rerun failed or timed-out rows loaded by --resume.",
    )
    parser.add_argument("--retries", type=nonnegative_int, default=0)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output-markdown", type=Path, default=DEFAULT_MARKDOWN)
    parser.add_argument(
        "--checkpoint-json",
        type=Path,
        default=DEFAULT_CHECKPOINT,
        help="Private atomic checkpoint used for resume; never a published result.",
    )
    return parser.parse_args()


def json_ready(value):
    """Convert Sage values recursively into canonical JSON values."""
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_ready(item) for item in value]
    try:
        return int(value)
    except (TypeError, ValueError):
        return str(value)


def matrix_degree(matrix):
    """Return the maximal Laurent L1 degree among matrix monomials."""
    degree = 0
    for row, column in matrix.nonzero_positions():
        exponents = matrix[row, column].dict()
        if exponents:
            degree = max(
                degree,
                max(abs(x_power) + abs(y_power) for x_power, y_power in exponents),
            )
    return int(degree)


def determinant_index(basis):
    """Return the index of a two-vector superlattice basis."""
    first, second = basis
    return abs(first[0] * second[1] - first[1] * second[0])


def verify_reference_superlattice(case, period_data):
    """Prove the paper basis belongs to and spans the computed superlattice."""
    if int(period_data.square_period) != case.reference_square_period:
        raise AssertionError("The computed square period does not match the paper.")
    period_vectors = set(period_data.period_vectors)
    missing = [vector for vector in case.superlattice_basis if vector not in period_vectors]
    if missing:
        raise AssertionError(f"Paper basis vectors are not computed periods: {missing}")
    paper_index = determinant_index(case.superlattice_basis)
    if paper_index != case.reference_q // 2:
        raise AssertionError("The paper basis index does not equal q/2.")
    computed_index = min(
        abs(u[0] * v[1] - u[1] * v[0])
        for position, u in enumerate(period_data.period_vectors)
        for v in period_data.period_vectors[position + 1 :]
        if u[0] * v[1] - u[1] * v[0]
    )
    if paper_index != computed_index:
        raise AssertionError("The paper basis does not span the full computed superlattice.")
    return paper_index


def sparse_product(left, right):
    """Multiply Laurent matrices while visiting only nonzero entries."""
    if left.ncols() != right.nrows():
        raise ValueError("Incompatible matrix dimensions for multiplication.")
    ring = left.base_ring()
    zero = ring.zero()
    result = zero_matrix(ring, left.nrows(), right.ncols())
    for row in range(left.nrows()):
        values = {}
        for middle in left.nonzero_positions_in_row(row):
            coefficient = left[row, middle]
            for column in right.nonzero_positions_in_row(middle):
                values[column] = (
                    values.get(column, zero) + coefficient * right[middle, column]
                )
        for column, value in values.items():
            if value != zero:
                result[row, column] = value
    return result


def sparse_product_equals(left, right, expected):
    """Return whether a sparse Laurent product equals an expected matrix."""
    if (
        left.ncols() != right.nrows()
        or left.nrows() != expected.nrows()
        or right.ncols() != expected.ncols()
    ):
        return False
    ring = left.base_ring()
    zero = ring.zero()
    for row in range(left.nrows()):
        values = {}
        for middle in left.nonzero_positions_in_row(row):
            coefficient = left[row, middle]
            for column in right.nonzero_positions_in_row(middle):
                values[column] = (
                    values.get(column, zero) + coefficient * right[middle, column]
                )
        actual_columns = {
            column for column, value in values.items() if value != zero
        }
        expected_columns = set(expected.nonzero_positions_in_row(row))
        if actual_columns != expected_columns:
            return False
        if any(values[column] != expected[row, column] for column in actual_columns):
            return False
    return True


def verify_chain_equations(result):
    """Verify both inverse equations and any explicitly requested forward maps."""
    num_x_checks = result.h_x_tilde.nrows()
    num_qubits = result.h_x_tilde.ncols()
    h_z, h_x = split_css_blocks(
        result.input_excitation_map,
        num_x_checks,
        num_qubits,
    )
    inverse = result.psi_inverse
    if not sparse_product_equals(
        h_x,
        inverse.psi_1_inverse,
        sparse_product(inverse.psi_0_inverse, result.h_x_tilde),
    ):
        raise AssertionError("The inverse degree-one chain equation failed.")
    if not sparse_product_equals(
        dagger_matrix(h_z),
        inverse.psi_2_inverse,
        sparse_product(inverse.psi_1_inverse, result.h_z_tilde_dagger),
    ):
        raise AssertionError("The inverse degree-two chain equation failed.")
    if result.psi is None:
        return
    forward = result.psi
    if not sparse_product_equals(
        result.h_x_tilde,
        forward.psi_1,
        sparse_product(forward.psi_0, h_x),
    ):
        raise AssertionError("The forward degree-one chain equation failed.")
    if not sparse_product_equals(
        result.h_z_tilde_dagger,
        forward.psi_2,
        sparse_product(forward.psi_1, dagger_matrix(h_z)),
    ):
        raise AssertionError("The forward degree-two chain equation failed.")
    for forward_map, inverse_map in (
        (forward.psi_2, inverse.psi_2_inverse),
        (forward.psi_1, inverse.psi_1_inverse),
        (forward.psi_0, inverse.psi_0_inverse),
    ):
        identity = identity_matrix(forward_map.base_ring(), forward_map.nrows())
        if not sparse_product_equals(
            forward_map,
            inverse_map,
            identity,
        ) or not sparse_product_equals(inverse_map, forward_map, identity):
            raise AssertionError("Forward and inverse chain maps are not mutual inverses.")


def qca_z_product_from_degree_two_map(h_z_tilde, degree_two_map, product_z_rank):
    """Construct the Z-side QCA product from the triangular certificate."""
    product_z_rank = int(product_z_rank)
    toric_rank = h_z_tilde.nrows() - product_z_rank
    if toric_rank == 0:
        return h_z_tilde[:, :]
    leading = degree_two_map[:product_z_rank, :product_z_rank]
    lower_left = degree_two_map[product_z_rank:, :product_z_rank]
    if leading != identity_matrix(degree_two_map.base_ring(), product_z_rank):
        raise AssertionError("The leading degree-two block is not identity.")
    if lower_left != zero_matrix(
        degree_two_map.base_ring(),
        toric_rank,
        product_z_rank,
    ):
        raise AssertionError("The lower-left degree-two block is not zero.")

    top = h_z_tilde[:product_z_rank, :]
    bottom = h_z_tilde[product_z_rank:, :]
    upper_right = degree_two_map[:product_z_rank, product_z_rank:]
    lower_right = degree_two_map[product_z_rank:, product_z_rank:]
    lower_inverse = dagger_matrix(lower_right).inverse()
    coupling = dagger_matrix(upper_right)
    bottom_product = sparse_product(
        lower_inverse,
        bottom + sparse_product(coupling, top),
    )
    result = zero_matrix(
        h_z_tilde.base_ring(),
        h_z_tilde.nrows(),
        h_z_tilde.ncols(),
    )
    result[:product_z_rank, :] = top
    result[product_z_rank:, :] = bottom_product
    return result


def verify_qca_certificate(result):
    """Verify symplectic decoupling using explicit or inverse-only maps."""
    num_x_checks = result.h_x_tilde.nrows()
    num_qubits = result.h_x_tilde.ncols()
    h_z, h_x = split_css_blocks(
        result.input_excitation_map,
        num_x_checks,
        num_qubits,
    )
    inverse = result.psi_inverse
    raw_x = sparse_product(h_x, inverse.psi_1_inverse)
    h_z_tilde = dagger_matrix(result.h_z_tilde_dagger)
    if result.psi is not None:
        raw_z = sparse_product(h_z, dagger_matrix(result.psi.psi_1))
    else:
        diagnostics = result.diagnostics
        post_product = qca_z_product_from_degree_two_map(
            h_z_tilde,
            diagnostics["post_phi_2_prime"],
            diagnostics["product_z_rank"],
        )
        xi_2_dagger_inverse = dagger_matrix(
            diagnostics["clearing_xi_2"]
        ).inverse()
        raw_z = sparse_product(xi_2_dagger_inverse, post_product)
    psi_0 = (
        result.psi.psi_0
        if result.psi is not None
        else inverse.psi_0_inverse.inverse()
    )
    if not sparse_product_equals(psi_0, raw_x, result.h_x_tilde):
        raise AssertionError("The stabilizer-redefined QCA X product failed.")
    if not sparse_product_equals(
        dagger_matrix(inverse.psi_2_inverse),
        raw_z,
        h_z_tilde,
    ):
        raise AssertionError("The stabilizer-redefined QCA Z product failed.")


def validate_case(case):
    """Run one complete paper-conformance validation."""
    started = time.perf_counter()
    record = case.as_reference_record()
    record.update(
        {
            "result_schema_version": RESULT_SCHEMA_VERSION,
            "implementation_fingerprint": RUN_IMPLEMENTATION_FINGERPRINT,
            **validation_environment(),
            "status": "running",
            "actual_square_period": None,
            "actual_q": None,
            "actual_p_x": None,
            "actual_p_z": None,
            "actual_t": None,
            "actual_psi_1_inverse_degree": None,
            "actual_runtime_seconds": None,
            "decoupling_runtime_seconds": None,
            "forward_maps_computed": False,
            "verification_status": "pending",
            "error": None,
        }
    )
    try:
        f, g = case.polynomials()
        excitation_map = build_two_generator_css_excitation_map(f, g)
        if excitation_map.dimensions() != (2, 4):
            raise AssertionError("The BB excitation map must have dimensions 2 x 4.")
        if not check_commutation(excitation_map, 2):
            raise AssertionError("The BB excitation map violates CSS commutation.")

        periods = find_anyon_preserving_superlattice_from_generators(
            f,
            g,
            max_period=case.reference_square_period,
        )
        index = verify_reference_superlattice(case, periods)
        coarse = coarse_grain_to_superlattice(
            excitation_map,
            *case.superlattice_basis,
        )
        actual_q = coarse.ncols() // 2
        if actual_q != 2 * index or actual_q != case.reference_q:
            raise AssertionError("The coarse qubit count does not match the paper.")

        compute_psi = actual_q <= EXPLICIT_FORWARD_Q_LIMIT
        decoupling_started = time.perf_counter()
        result = solve_decoupling_unitary(
            coarse,
            num_x_checks=coarse.nrows() // 2,
            num_qubits=actual_q,
            compute_psi=compute_psi,
        )
        decoupling_seconds = time.perf_counter() - decoupling_started
        verify_chain_equations(result)
        verify_qca_certificate(result)

        actual_p_x = int(result.diagnostics["product_x_rank"])
        actual_p_z = int(result.diagnostics["product_z_rank"])
        actual_t_x = result.h_x_tilde.nrows() - actual_p_x
        actual_t_z = result.h_z_tilde_dagger.ncols() - actual_p_z
        if actual_t_x != actual_t_z:
            raise AssertionError("The X and Z standard complexes have unequal toric counts.")
        actual_t = int(actual_t_x)
        actual = (actual_q, actual_p_x, actual_p_z, actual_t)
        reference = (
            case.reference_q,
            case.reference_p_x,
            case.reference_p_z,
            case.reference_t,
        )
        if actual != reference:
            raise AssertionError(f"Computed (q,p_X,p_Z,t)={actual}, expected {reference}.")

        record.update(
            {
                "status": "ok",
                "actual_square_period": int(periods.square_period),
                "actual_q": actual_q,
                "actual_p_x": actual_p_x,
                "actual_p_z": actual_p_z,
                "actual_t": actual_t,
                "actual_psi_1_inverse_degree": matrix_degree(
                    result.psi_inverse.psi_1_inverse
                ),
                "decoupling_runtime_seconds": decoupling_seconds,
                "forward_maps_computed": result.psi is not None,
                "verification_status": "passed",
            }
        )
    except Exception as error:  # Persist reproducible row failures.
        record.update(
            {
                "status": "failed",
                "verification_status": "failed",
                "error": f"{type(error).__name__}: {error}",
                "traceback": traceback.format_exc(limit=40),
            }
        )
    record["actual_runtime_seconds"] = time.perf_counter() - started
    return json_ready(record)


def worker(case, connection):
    """Validate a row in a child process."""
    try:
        connection.send(validate_case(case))
    finally:
        connection.close()


def timeout_record(case, elapsed, timeout_seconds):
    """Return a durable timeout record."""
    record = case.as_reference_record()
    record.update(
        {
            "result_schema_version": RESULT_SCHEMA_VERSION,
            "implementation_fingerprint": RUN_IMPLEMENTATION_FINGERPRINT,
            **validation_environment(),
            "status": "timeout",
            "actual_square_period": None,
            "actual_q": None,
            "actual_p_x": None,
            "actual_p_z": None,
            "actual_t": None,
            "actual_psi_1_inverse_degree": None,
            "actual_runtime_seconds": elapsed,
            "decoupling_runtime_seconds": None,
            "forward_maps_computed": False,
            "verification_status": "failed",
            "error": f"TimeoutError: exceeded {timeout_seconds:.3f} seconds",
        }
    )
    return record


def run_with_optional_timeout(case, timeout_seconds):
    """Run a row directly or enforce a process-level timeout."""
    if timeout_seconds is None:
        return validate_case(case)
    context = multiprocessing.get_context("fork")
    parent_connection, child_connection = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(case, child_connection))
    started = time.perf_counter()
    try:
        process.start()
        child_connection.close()
        process.join(timeout_seconds)
        elapsed = time.perf_counter() - started
        if process.is_alive():
            return timeout_record(case, elapsed, timeout_seconds)
        if process.exitcode != 0 or not parent_connection.poll(1.0):
            record = timeout_record(case, elapsed, timeout_seconds)
            record["status"] = "failed"
            record["error"] = (
                f"ProcessExitError: child exited with code {process.exitcode} "
                "without returning a record"
            )
            return record
        return parent_connection.recv()
    finally:
        child_connection.close()
        if process.pid is not None and process.is_alive():
            process.terminate()
            process.join()
        parent_connection.close()


def atomic_write(path, text):
    """Atomically replace one UTF-8 result file."""
    path = Path(path)
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


def build_metadata(cases, args, records):
    """Return runtime and completeness metadata for one checkpoint."""
    metadata = runtime_metadata(
        adapter="decoupling",
        suite="paper_bb_codes",
        source_root=SOURCE_ROOT,
    )
    metadata.pop("source_root", None)
    expected = [case.identifier for case in cases]
    completed = {record["id"] for record in records}
    metadata.update(
        {
            "result_schema_version": RESULT_SCHEMA_VERSION,
            "implementation_fingerprint": RUN_IMPLEMENTATION_FINGERPRINT,
            "row_count": len(cases),
            "catalog_row_count": len(BB_CODE_INSTANCES),
            "row_ids": expected,
            "runtime_cutoff_reference_id": RUNTIME_CUTOFF_REFERENCE_ID,
            "excluded_longer_runtime_row_ids": sorted(
                EXCLUDED_LONGER_RUNTIME_ROW_IDS
            ),
            "row_timeout_seconds": args.row_timeout_seconds,
            "explicit_forward_q_limit": EXPLICIT_FORWARD_Q_LIMIT,
            "complete": completed == set(expected),
            "source_unchanged_during_run": (
                compute_implementation_fingerprint()
                == RUN_IMPLEMENTATION_FINGERPRINT
            ),
            "all_records_current": completed == set(expected)
            and all(
                resume_record_is_current(
                    next(case for case in cases if case.identifier == record["id"]),
                    record,
                )
                for record in records
            ),
            "paper_runtime_included": False,
        }
    )
    metadata["all_passed"] = (
        metadata["complete"]
        and metadata["source_unchanged_during_run"]
        and metadata["all_records_current"]
        and all(record["status"] == "ok" for record in records)
    )
    return metadata


def csv_text(metadata, records):
    """Render canonical row records and runtime metadata as CSV."""
    columns = [
        "id",
        "collection",
        "row",
        "f",
        "g",
        "reference_square_period",
        "actual_square_period",
        "superlattice_basis",
        "reference_q",
        "actual_q",
        "reference_p_x",
        "actual_p_x",
        "reference_p_z",
        "actual_p_z",
        "reference_t",
        "actual_t",
        "actual_psi_1_inverse_degree",
        "decoupling_runtime_seconds",
        "actual_runtime_seconds",
        "forward_maps_computed",
        "verification_status",
        "status",
        "error",
        "result_schema_version",
        "implementation_fingerprint",
        "validation_sage_version",
        "validation_python_version",
        "validation_platform",
        "adapter",
        "suite",
        "source_revision",
        "source_tree_dirty",
        "sage_version",
        "python_version",
        "platform",
        "source_unchanged_during_run",
        "all_records_current",
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
    for record in records:
        flattened = dict(record)
        for key, value in metadata.items():
            flattened.setdefault(key, value)
        flattened["superlattice_basis"] = json.dumps(
            record["superlattice_basis"],
            separators=(",", ":"),
        )
        writer.writerow(flattened)
    return stream.getvalue()


def markdown_text(metadata, records):
    """Render a human-readable result table."""
    lines = [
        "# Paper BB-Code Decoupling Reproduction",
        "",
        f"- Rows: {len(records)} / {metadata['row_count']}",
        f"- Full paper catalog rows: {metadata['catalog_row_count']}",
        f"- Runtime cutoff reference: `{metadata['runtime_cutoff_reference_id']}`",
        "- Rows with longer paper runtimes excluded: "
        f"`{metadata['excluded_longer_runtime_row_ids']}`",
        f"- Complete: `{metadata['complete']}`",
        f"- All passed: `{metadata['all_passed']}`",
        f"- Sage: `{metadata['sage_version']}`",
        f"- Python: `{metadata['python_version']}`",
        f"- Platform: `{metadata['platform']}`",
        f"- Adapter: `{metadata['adapter']}`",
        f"- Suite: `{metadata['suite']}`",
        f"- Source revision: `{metadata['source_revision']}`",
        f"- Source tree dirty: `{metadata['source_tree_dirty']}`",
        f"- Result schema: `{metadata['result_schema_version']}`",
        f"- Implementation fingerprint: `{metadata['implementation_fingerprint']}`",
        f"- Source unchanged during run: `{metadata['source_unchanged_during_run']}`",
        f"- All records current: `{metadata['all_records_current']}`",
        f"- Paper runtime included: `{metadata['paper_runtime_included']}`",
        "- Paper timing data are intentionally not stored.",
        "",
        "| ID | f | g | L | Basis | q | p_X | p_Z | t | deg(psi_1^-1) | Decoupling s | Total s | Forward maps | Verification | Status |",
        "| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for record in records:
        lines.append(
            "| "
            + " | ".join(
                [
                    record["id"],
                    f"`{record['f']}`",
                    f"`{record['g']}`",
                    str(record["actual_square_period"] or ""),
                    str(record["superlattice_basis"]),
                    str(record["actual_q"] or ""),
                    str(record["actual_p_x"] if record["actual_p_x"] is not None else ""),
                    str(record["actual_p_z"] if record["actual_p_z"] is not None else ""),
                    str(record["actual_t"] if record["actual_t"] is not None else ""),
                    str(
                        record["actual_psi_1_inverse_degree"]
                        if record["actual_psi_1_inverse_degree"] is not None
                        else ""
                    ),
                    (
                        f"{record['decoupling_runtime_seconds']:.6f}"
                        if record["decoupling_runtime_seconds"] is not None
                        else ""
                    ),
                    f"{record['actual_runtime_seconds']:.6f}",
                    str(record["forward_maps_computed"]),
                    record["verification_status"],
                    record["status"],
                ]
            )
            + " |"
        )
    return "\n".join(lines) + "\n"


def result_payload(cases, args, records):
    """Return sorted rows, metadata, and canonical JSON text."""
    records = sorted(records, key=lambda record: (record["collection"], record["row"]))
    metadata = build_metadata(cases, args, records)
    payload = json.dumps(
        {"metadata": json_ready(metadata), "rows": records},
        indent=2,
        sort_keys=True,
    ) + "\n"
    return records, metadata, payload


def write_resume_checkpoint(args, cases, records):
    """Atomically persist private resume state without publishing partial results."""
    _, metadata, payload = result_payload(cases, args, records)
    atomic_write(args.checkpoint_json, payload)
    return metadata


def publish_results(args, cases, records):
    """Atomically publish JSON, CSV, and Markdown after every row has passed."""
    records, metadata, payload = result_payload(cases, args, records)
    atomic_write(args.output_json, payload)
    atomic_write(args.output_csv, csv_text(metadata, records))
    atomic_write(args.output_markdown, markdown_text(metadata, records))
    return metadata


def load_resume_records(path):
    """Load prior row records from a canonical JSON checkpoint."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = payload.get("rows")
    if not isinstance(rows, list):
        raise ValueError("Resume JSON does not contain a row list.")
    return {record["id"]: record for record in rows}


def resume_record_is_current(case, record):
    """Return whether a resume row matches this source and paper input exactly."""
    reference = case.as_reference_record()
    if record.get("result_schema_version") != RESULT_SCHEMA_VERSION:
        return False
    if record.get("implementation_fingerprint") != RUN_IMPLEMENTATION_FINGERPRINT:
        return False
    if any(
        record.get(key) != value
        for key, value in validation_environment().items()
    ):
        return False
    return all(record.get(key) == value for key, value in reference.items())


def configure_selected_output_paths(args, cases, selectors):
    """Keep selected-row outputs separate from the canonical 59-row files."""
    if selectors is None:
        return
    identifiers = "_".join(case.identifier for case in cases)
    if len(identifiers) > 120:
        identifiers = hashlib.sha256(identifiers.encode("utf-8")).hexdigest()[:16]
    stem = f"bb_codes.selection.{identifiers}"
    if args.output_json == DEFAULT_JSON:
        args.output_json = RESULT_DIRECTORY / f"{stem}.json"
    if args.output_csv == DEFAULT_CSV:
        args.output_csv = RESULT_DIRECTORY / f"{stem}.csv"
    if args.output_markdown == DEFAULT_MARKDOWN:
        args.output_markdown = RESULT_DIRECTORY / f"{stem}.md"
    if args.checkpoint_json == DEFAULT_CHECKPOINT:
        args.checkpoint_json = RESULT_DIRECTORY / f".{stem}.checkpoint.json"


def main():
    """Run selected rows, checkpoint atomically, and enforce completeness."""
    args = parse_args()
    selectors = (
        [item.strip() for item in args.rows.split(",") if item.strip()]
        if args.rows
        else None
    )
    cases = (
        select_bb_code_instances(selectors)
        if selectors is not None
        else DEFAULT_VALIDATION_CASES
    )
    configure_selected_output_paths(args, cases, selectors)
    full_publication_run = selectors is None
    records_by_id = {}
    if args.resume:
        resume_path = (
            args.checkpoint_json
            if args.checkpoint_json.exists()
            else args.output_json
        )
        if resume_path.exists():
            loaded = load_resume_records(resume_path)
            records_by_id.update(
                {
                    case.identifier: loaded[case.identifier]
                    for case in cases
                    if case.identifier in loaded
                    and resume_record_is_current(case, loaded[case.identifier])
                }
            )

    selected_ids = {case.identifier for case in cases}
    records_by_id = {
        identifier: record
        for identifier, record in records_by_id.items()
        if identifier in selected_ids
    }
    write_resume_checkpoint(args, cases, list(records_by_id.values()))

    for case in cases:
        previous = records_by_id.get(case.identifier)
        if previous is not None:
            if previous.get("status") == "ok":
                continue
            if not args.retry_failed:
                continue
        attempts = 1 + args.retries
        for attempt in range(1, attempts + 1):
            print(
                f"[bb] {case.identifier} attempt={attempt}/{attempts}",
                file=sys.stderr,
                flush=True,
            )
            record = run_with_optional_timeout(case, args.row_timeout_seconds)
            records_by_id[case.identifier] = record
            write_resume_checkpoint(args, cases, list(records_by_id.values()))
            if record["status"] == "ok":
                break

    records = list(records_by_id.values())
    metadata = build_metadata(cases, args, records)
    failures = [record for record in records if record["status"] != "ok"]
    missing = selected_ids - set(records_by_id)
    if failures or missing or not metadata["all_passed"]:
        raise SystemExit(1)
    if full_publication_run and (
        len(records) != len(DEFAULT_VALIDATION_CASES)
        or not metadata["complete"]
        or not metadata["all_passed"]
    ):
        raise SystemExit(1)
    metadata = publish_results(args, cases, records)
    if args.checkpoint_json.exists():
        args.checkpoint_json.unlink()


if __name__ == "__main__":
    main()
