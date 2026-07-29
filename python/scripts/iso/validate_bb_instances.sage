"""Validate the full isomorphism pipeline on every BB Table I instance."""

from pathlib import Path
import argparse
import json
import multiprocessing
import sys
import time
import traceback

from sage.all import identity_matrix, zero_matrix


REPO_ROOT = Path(__file__).resolve().parents[2]
repo_root = str(REPO_ROOT)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from decoder_core.bb_cases import BB_ROWS  # noqa: E402
from decoder_core.benchmark_metadata import runtime_metadata  # noqa: E402
from isomorphism import (  # noqa: E402
    choose_smallest_oblique_cell,
    construct_excitation_map,
    decouple_coarse_matrix,
    oblique_coarse_grain,
    periods_from_generators,
    x,
    y,
)
from isomorphism.css import antipode_matrix, split_css_blocks  # noqa: E402


DEFAULT_OUTPUT_JSON = REPO_ROOT / "results" / "isomorphism" / "isomorphism_bb_full_validation.json"
DEFAULT_OUTPUT_MD = REPO_ROOT / "results" / "isomorphism" / "isomorphism_bb_full_validation.md"


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", help="Comma-separated BB row numbers. Defaults to all 28.")
    parser.add_argument("--row-timeout-seconds", type=float, default=1800.0)
    parser.add_argument(
        "--qca-check",
        choices=("product", "inverse", "skip"),
        default="product",
        help=(
            "'product' verifies the QCA epsilon product from the inverse degree-two "
            "certificate without constructing forward maps; 'inverse' also constructs "
            "and checks the forward maps; 'skip' validates only the supplement "
            "period/coarse-graining and inverse chain-map algorithms."
        ),
    )
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT_JSON)
    parser.add_argument("--output-md", type=Path, default=DEFAULT_OUTPUT_MD)
    parser.add_argument("--stop-on-failure", action="store_true")
    return parser.parse_args()


def selected_rows(row_text):
    """Return the selected BB rows in table order."""
    if not row_text:
        return BB_ROWS
    wanted = {int(part.strip()) for part in row_text.split(",") if part.strip()}
    if not wanted:
        raise ValueError("--rows must contain at least one row number.")
    rows = tuple(row for row in BB_ROWS if row.row in wanted)
    missing = sorted(wanted - {row.row for row in rows})
    if missing:
        raise ValueError(f"Unknown BB row numbers: {missing}")
    return rows


def bb_generators(row):
    """Return the two BB generator polynomials for one table row."""
    a = x ** row.a_exponents[0] * y ** row.a_exponents[1]
    b = x ** row.b_exponents[0] * y ** row.b_exponents[1]
    return 1 + x + a, 1 + y + b


def json_ready(value):
    """Convert Sage/Python values used by the runner into JSON-safe data."""
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, tuple):
        return [json_ready(item) for item in value]
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    try:
        return int(value)
    except (TypeError, ValueError):
        return str(value)


def emit_stage(progress_queue, event, stage):
    """Emit one stage progress event to the parent process if requested."""
    if progress_queue is not None:
        progress_queue.put((event, json_ready(stage)))


def timed(stages, name, thunk, *, progress_queue=None):
    """Run one validation stage and append a timed stage record."""
    running = {"stage": name, "status": "running", "seconds": None, "error": None}
    emit_stage(progress_queue, "stage_start", running)
    started = time.perf_counter()
    try:
        value = thunk()
    except BaseException as exc:  # noqa: BLE001 - validation must persist failures
        stage = {
            "stage": name,
            "status": "failed",
            "seconds": time.perf_counter() - started,
            "error": f"{type(exc).__name__}: {exc}",
        }
        stages.append(stage)
        emit_stage(progress_queue, "stage_finish", stage)
        raise
    stage = {
        "stage": name,
        "status": "ok",
        "seconds": time.perf_counter() - started,
        "error": None,
    }
    stages.append(stage)
    emit_stage(progress_queue, "stage_finish", stage)
    return value


def verify_period(row, periods):
    """Validate the discovered period data against the table period."""
    if int(periods.square_period) != int(row.paper_period):
        raise AssertionError(
            f"computed square period {periods.square_period} != paper period {row.paper_period}"
        )
    return True


def verify_chain_maps(decoupled):
    """Validate inverse and optional forward decoupling identities."""
    coarse_epsilon = decoupled.source_coarse_matrix
    inverse_maps = decoupled.inverse_maps
    maps = decoupled.maps
    num_x_checks = coarse_epsilon.nrows() // 2
    num_qubits = coarse_epsilon.ncols() // 2
    HZ, HX = split_css_blocks(coarse_epsilon, num_x_checks, num_qubits)
    if not sparse_product_equals(
        HX,
        inverse_maps.phi1_inverse,
        sparse_product(
            inverse_maps.phi0_inverse,
            decoupled.hx_standard,
        ),
    ):
        raise AssertionError("The inverse degree-one chain equation failed.")
    if not sparse_product_equals(
        antipode_matrix(HZ),
        inverse_maps.phi2_inverse,
        sparse_product(
            inverse_maps.phi1_inverse,
            antipode_matrix(decoupled.hz_standard),
        ),
    ):
        raise AssertionError("The inverse degree-two chain equation failed.")
    if maps is not None:
        hz_source_dagger = antipode_matrix(HZ)
        hz_standard_dagger = antipode_matrix(decoupled.hz_standard)
        if not sparse_product_equals(
            decoupled.hx_standard,
            maps.phi1,
            sparse_product(maps.phi0, HX),
        ):
            raise AssertionError("The forward degree-zero chain equation failed.")
        if not sparse_product_equals(
            hz_standard_dagger,
            maps.phi2,
            sparse_product(maps.phi1, hz_source_dagger),
        ):
            raise AssertionError("The forward degree-two chain equation failed.")

        inverse_pairs = (
            ("degree two", maps.phi2, inverse_maps.phi2_inverse),
            ("degree one", maps.phi1, inverse_maps.phi1_inverse),
            ("degree zero", maps.phi0, inverse_maps.phi0_inverse),
        )
        for label, forward, inverse in inverse_pairs:
            identity = identity_matrix(forward.base_ring(), forward.nrows())
            if not sparse_product_equals(forward, inverse, identity):
                raise AssertionError(f"The {label} maps are not mutual inverses.")
            if not sparse_product_equals(inverse, forward, identity):
                raise AssertionError(f"The {label} maps are not mutual inverses.")
    return True


def sparse_product(A, B):
    """Return ``A * B`` using only nonzero matrix entries."""
    if A.ncols() != B.nrows():
        raise ValueError("Incompatible matrix dimensions for multiplication.")
    ring = A.base_ring()
    result = zero_matrix(ring, A.nrows(), B.ncols())
    zero = ring.zero()
    for row in range(A.nrows()):
        values = {}
        for middle in A.nonzero_positions_in_row(row):
            left = A[row, middle]
            for col in B.nonzero_positions_in_row(middle):
                values[col] = values.get(col, zero) + left * B[middle, col]
        for col, value in values.items():
            if value != zero:
                result[row, col] = value
    return result


def sparse_product_equals(A, B, C):
    """Return whether ``A * B == C`` using only nonzero matrix entries."""
    if A.ncols() != B.nrows() or A.nrows() != C.nrows() or B.ncols() != C.ncols():
        return False
    ring = A.base_ring()
    zero = ring.zero()
    for row in range(A.nrows()):
        values = {}
        for middle in A.nonzero_positions_in_row(row):
            left = A[row, middle]
            for col in B.nonzero_positions_in_row(middle):
                values[col] = values.get(col, zero) + left * B[middle, col]
        expected_cols = set(C.nonzero_positions_in_row(row))
        actual_cols = {col for col, value in values.items() if value != zero}
        if actual_cols != expected_cols:
            return False
        for col in actual_cols:
            if values[col] != C[row, col]:
                return False
    return True


def qca_z_product_from_phi2_inverse(
    hz_standard,
    phi2_inverse,
    product_z_rank,
):
    """Return the QCA Z block from the inverse degree-two certificate."""
    product_z_rank = int(product_z_rank)
    toric_rank = hz_standard.nrows() - product_z_rank
    if toric_rank == 0:
        return hz_standard[:, :]
    leading = phi2_inverse[:product_z_rank, :product_z_rank]
    lower_left = phi2_inverse[product_z_rank:, :product_z_rank]
    if leading != identity_matrix(phi2_inverse.base_ring(), product_z_rank):
        raise AssertionError("The leading inverse degree-two block is not identity.")
    if lower_left != zero_matrix(
        phi2_inverse.base_ring(),
        toric_rank,
        product_z_rank,
    ):
        raise AssertionError("The lower-left inverse degree-two block is not zero.")

    top_standard = hz_standard[:product_z_rank, :]
    bottom_standard = hz_standard[product_z_rank:, :]
    upper_right = phi2_inverse[:product_z_rank, product_z_rank:]
    lower_right = phi2_inverse[product_z_rank:, product_z_rank:]
    lower_right_dagger_inverse = antipode_matrix(lower_right).inverse()
    coupling = antipode_matrix(upper_right)
    bottom_product = sparse_product(
        lower_right_dagger_inverse,
        bottom_standard + sparse_product(coupling, top_standard),
    )
    result = zero_matrix(
        hz_standard.base_ring(),
        hz_standard.nrows(),
        hz_standard.ncols(),
    )
    result[:product_z_rank, :] = top_standard
    result[product_z_rank:, :] = bottom_product
    return result


def qca_summary(decoupled, *, qca_check):
    """Return QCA epsilon-product checks for one BB row."""
    inverse_maps = decoupled.inverse_maps
    maps = decoupled.maps
    diagnostics = decoupled.diagnostics
    num_x_checks = decoupled.hx_standard.nrows()
    num_qubits = decoupled.hx_standard.ncols()
    HZ, HX = split_css_blocks(
        decoupled.source_coarse_matrix,
        num_x_checks,
        num_qubits,
    )
    raw_x = sparse_product(HX, inverse_maps.phi1_inverse)
    if qca_check == "inverse":
        raw_z = sparse_product(HZ, antipode_matrix(maps.phi1))
        verification_method = "explicit_forward_degree_one_map"
    else:
        post_z_product = qca_z_product_from_phi2_inverse(
            decoupled.hz_standard,
            diagnostics["post_inverse_maps"].phi2_inverse,
            diagnostics["product_z_rank"],
        )
        xi2_dagger_inverse = antipode_matrix(
            diagnostics["clearing_xi2"],
        ).inverse()
        raw_z = sparse_product(xi2_dagger_inverse, post_z_product)
        verification_method = "inverse_degree_two_triangular_certificate"
    raw_matches = (
        raw_x == decoupled.hx_standard
        and raw_z == decoupled.hz_standard
    )
    phi0 = maps.phi0 if maps is not None else inverse_maps.phi0_inverse.inverse()
    redefined_matches = sparse_product_equals(
        phi0,
        raw_x,
        decoupled.hx_standard,
    ) and sparse_product_equals(
        antipode_matrix(inverse_maps.phi2_inverse),
        raw_z,
        decoupled.hz_standard,
    )
    if not redefined_matches:
        raise AssertionError("The stabilizer-redefined QCA product is not standard.")
    return {
        "qca_matrix_shape": (2 * num_qubits, 2 * num_qubits),
        "standard_epsilon_shape": (
            decoupled.hx_standard.nrows() + decoupled.hz_standard.nrows(),
            2 * num_qubits,
        ),
        "forward_maps_computed": maps is not None,
        "raw_qca_product_matches_standard": raw_matches,
        "after_stabilizer_redefinition_matches_standard": redefined_matches,
        "qca_verification_method": verification_method,
    }


def validate_row(row, qca_check="product", progress_queue=None):
    """Run the full validation pipeline for one BB row."""
    started = time.perf_counter()
    stages = []
    result = {
        "row": row.row,
        "a_exponents": row.a_exponents,
        "b_exponents": row.b_exponents,
        "paper_period": row.paper_period,
        "status": "ok",
        "total_seconds": None,
        "stages": stages,
        "error": None,
    }
    try:
        f, g = bb_generators(row)
        periods = timed(
            stages,
            "period_search",
            lambda: periods_from_generators(f, g, max_period=row.paper_period),
            progress_queue=progress_queue,
        )
        result["computed_square_period"] = int(periods.square_period)
        result["period_count"] = len(periods.vectors)
        timed(
            stages,
            "verify_period",
            lambda: verify_period(row, periods),
            progress_queue=progress_queue,
        )
        cell = timed(
            stages,
            "choose_cell",
            lambda: choose_smallest_oblique_cell(periods.vectors),
            progress_queue=progress_queue,
        )
        result["cell"] = cell
        epsilon = timed(
            stages,
            "construct_excitation_map",
            lambda: construct_excitation_map(f, g),
            progress_queue=progress_queue,
        )
        result["epsilon_shape"] = (epsilon.nrows(), epsilon.ncols())
        coarse = timed(
            stages,
            "coarse_grain",
            lambda: oblique_coarse_grain(epsilon, *cell),
            progress_queue=progress_queue,
        )
        result["coarse_epsilon_shape"] = (coarse.nrows(), coarse.ncols())
        decoupled = timed(
            stages,
            "decouple_total",
            lambda: decouple_coarse_matrix(
                coarse,
                num_x_checks=coarse.nrows() // 2,
                num_qubits=coarse.ncols() // 2,
                compute_forward_maps=qca_check == "inverse",
            ),
            progress_queue=progress_queue,
        )
        result["product_ranks"] = (
            decoupled.diagnostics["product_x_rank"],
            decoupled.diagnostics["product_z_rank"],
        )
        result["phi1_inverse_shape"] = (
            decoupled.inverse_maps.phi1_inverse.nrows(),
            decoupled.inverse_maps.phi1_inverse.ncols(),
        )
        timed(
            stages,
            "verify_chain_maps",
            lambda: verify_chain_maps(decoupled),
            progress_queue=progress_queue,
        )
        result["qca_check"] = qca_check
        if qca_check in ("product", "inverse"):
            result.update(
                timed(
                    stages,
                    "verify_qca_decoupling",
                    lambda: qca_summary(
                        decoupled,
                        qca_check=qca_check,
                    ),
                    progress_queue=progress_queue,
                )
            )
        else:
            result.update(
                {
                    "qca_matrix_shape": None,
                    "standard_epsilon_shape": (
                        decoupled.hx_standard.nrows()
                        + decoupled.hz_standard.nrows(),
                        2 * decoupled.hx_standard.ncols(),
                    ),
                    "forward_maps_computed": decoupled.maps is not None,
                    "raw_qca_product_matches_standard": None,
                    "after_stabilizer_redefinition_matches_standard": None,
                    "qca_verification_method": "skipped",
                }
            )
    except BaseException as exc:  # noqa: BLE001 - caller needs the captured bug
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc(limit=40)
    result["total_seconds"] = time.perf_counter() - started
    result = json_ready(result)
    if progress_queue is not None:
        progress_queue.put(("final", result))
    return result


def row_worker(row, qca_check, queue):
    """Run one row in a child process and return JSON-safe data."""
    validate_row(row, qca_check=qca_check, progress_queue=queue)


def base_row_result(row):
    """Return the parent-side partial result for one row."""
    return {
        "row": row.row,
        "a_exponents": list(row.a_exponents),
        "b_exponents": list(row.b_exponents),
        "paper_period": row.paper_period,
        "status": "running",
        "total_seconds": None,
        "stages": [],
        "error": None,
    }


def update_partial_stage(partial, stage):
    """Update one parent-side partial result with a stage progress record."""
    stages = partial.setdefault("stages", [])
    for index in range(len(stages) - 1, -1, -1):
        if stages[index].get("stage") == stage.get("stage"):
            if stages[index].get("status") == "running" or stage.get("status") != "running":
                stages[index] = stage
                return
    stages.append(stage)


def drain_queue(queue, partial):
    """Drain child progress events and return the final result if available."""
    final = None
    while not queue.empty():
        event, payload = queue.get()
        if event in ("stage_start", "stage_finish"):
            update_partial_stage(partial, payload)
        elif event == "final":
            final = payload
        else:
            raise ValueError(f"Unknown worker event: {event}")
    return final


def timeout_result(row, elapsed, timeout_seconds):
    """Return a timeout record for one row."""
    return {
        "row": row.row,
        "a_exponents": list(row.a_exponents),
        "b_exponents": list(row.b_exponents),
        "paper_period": row.paper_period,
        "status": "timeout",
        "total_seconds": elapsed,
        "stages": [],
        "error": f"TimeoutError: exceeded {float(timeout_seconds):.1f} seconds",
    }


def run_row_with_timeout(row, timeout_seconds, *, qca_check):
    """Run one BB row with a process-level timeout."""
    started = time.perf_counter()
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    partial = base_row_result(row)
    process = context.Process(target=row_worker, args=(row, qca_check, queue))
    process.start()
    final = None
    while process.is_alive():
        process.join(0.25)
        drained = drain_queue(queue, partial)
        if drained is not None:
            final = drained
        elapsed = time.perf_counter() - started
        if elapsed > float(timeout_seconds):
            process.terminate()
            process.join()
            drained = drain_queue(queue, partial)
            if drained is not None:
                final = drained
            if final is not None:
                return final
            result = timeout_result(row, elapsed, timeout_seconds)
            result["stages"] = partial.get("stages", []) + result["stages"]
            result["stages"].append(
                {
                    "stage": "row_timeout",
                    "status": "failed",
                    "seconds": elapsed,
                    "error": f"TimeoutError: exceeded {float(timeout_seconds):.1f} seconds",
                }
            )
            return result
    elapsed = time.perf_counter() - started
    drained = drain_queue(queue, partial)
    if drained is not None:
        final = drained
    if final is not None:
        return final
    return {
        "row": row.row,
        "a_exponents": list(row.a_exponents),
        "b_exponents": list(row.b_exponents),
        "paper_period": row.paper_period,
        "status": "failed",
        "total_seconds": elapsed,
        "stages": partial.get("stages", [])
        + [
            {
                "stage": "row_process",
                "status": "failed",
                "seconds": elapsed,
                "error": f"ProcessExitError: child exited with code {process.exitcode}",
            }
        ],
        "error": f"ProcessExitError: child exited with code {process.exitcode}",
    }


def validation_metadata(rows, args):
    """Return metadata for one full BB validation run."""
    metadata = runtime_metadata(
        adapter="isomorphism",
        suite="bb_full_validation",
        source_root=REPO_ROOT,
    )
    # Durable public artifacts identify the revision and runtime, not the
    # machine-specific checkout path used to produce them.
    metadata.pop("source_root", None)
    metadata.update(
        {
            "row_count": len(rows),
            "rows": [row.row for row in rows],
            "row_timeout_seconds": float(args.row_timeout_seconds),
            "qca_check": args.qca_check,
            "checks": [
                "period_matches_supplement_table",
                "supplement_oblique_coarse_grain",
                "supplement_inverse_chain_equations",
            ],
        }
    )
    if args.qca_check == "inverse":
        metadata["checks"].extend(
            [
                "forward_inverse_map_pairs",
                "qca_after_stabilizer_redefinition",
            ]
        )
    elif args.qca_check == "product":
        metadata["checks"].append(
            "inverse_degree_two_triangular_qca_certificate"
        )
    return metadata


def write_json(path, metadata, rows):
    """Write validation rows to JSON."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": json_ready(metadata), "rows": json_ready(rows)}
    target.write_text(json.dumps(payload, indent=2) + "\n")


def stage_seconds(row, stage_name):
    """Return the seconds value for one stage, if present."""
    for stage in row.get("stages", []):
        if stage.get("stage") == stage_name:
            return stage.get("seconds")
    return None


def format_seconds(value):
    """Return a compact seconds string for Markdown output."""
    if value is None:
        return ""
    return f"{float(value):.3f}"


def write_markdown(path, metadata, rows):
    """Write a Markdown summary table for the validation run."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    ok_count = sum(1 for row in rows if row.get("status") == "ok")
    lines = [
        "# Isomorphism BB Full Validation",
        "",
        f"- Source revision: `{metadata['source_revision']}`",
        f"- Sage version: `{metadata['sage_version']}`",
        f"- Rows checked: {len(rows)}",
        f"- Rows passed: {ok_count}",
        f"- Row timeout seconds: {metadata['row_timeout_seconds']}",
        f"- QCA check: `{metadata.get('qca_check', 'explicit')}`",
        "",
        "| Row | L | Status | Total s | Period s | Coarse s | Decouple s | QCA s | Cell | phi1_inverse shape | Forward maps | Raw QCA=standard | Redefined QCA=standard | Verification method | Error |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        error = (row.get("error") or "").replace("|", "\\|").replace("\n", " ")
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row.get("row")),
                    str(row.get("paper_period")),
                    str(row.get("status")),
                    format_seconds(row.get("total_seconds")),
                    format_seconds(stage_seconds(row, "period_search")),
                    format_seconds(stage_seconds(row, "coarse_grain")),
                    format_seconds(stage_seconds(row, "decouple_total")),
                    format_seconds(stage_seconds(row, "verify_qca_decoupling")),
                    str(row.get("cell", "")),
                    str(row.get("phi1_inverse_shape", "")),
                    str(row.get("forward_maps_computed", "")),
                    str(row.get("raw_qca_product_matches_standard", "")),
                    str(row.get("after_stabilizer_redefinition_matches_standard", "")),
                    str(row.get("qca_verification_method", "")),
                    error,
                ]
            )
            + " |"
        )
    target.write_text("\n".join(lines) + "\n")


def write_outputs(args, metadata, rows):
    """Write both JSON and Markdown outputs."""
    write_json(args.output_json, metadata, rows)
    write_markdown(args.output_md, metadata, rows)


def main():
    """Run selected BB validation rows and checkpoint the result files."""
    args = parse_args()
    rows = selected_rows(args.rows)
    metadata = validation_metadata(rows, args)
    results = []
    write_outputs(args, metadata, results)
    for row in rows:
        print(
            f"[bb-validation] start row={row.row} L={row.paper_period}",
            file=sys.stderr,
            flush=True,
        )
        result = run_row_with_timeout(
            row,
            args.row_timeout_seconds,
            qca_check=args.qca_check,
        )
        results.append(result)
        write_outputs(args, metadata, results)
        print(
            f"[bb-validation] done row={row.row} status={result['status']} "
            f"seconds={float(result['total_seconds']):.3f}",
            file=sys.stderr,
            flush=True,
        )
        if args.stop_on_failure and result["status"] != "ok":
            break
    failed = [row for row in results if row.get("status") != "ok"]
    print(
        f"[bb-validation] completed rows={len(results)} failed={len(failed)} "
        f"json={args.output_json} md={args.output_md}",
        file=sys.stderr,
        flush=True,
    )
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
