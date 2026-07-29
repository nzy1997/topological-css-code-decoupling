"""Benchmark 6.6.6 color-code unitary-decouple decoding against BPOSD."""

from pathlib import Path
import argparse
import csv
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[2]
repo_root = str(REPO_ROOT)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from decoder_core import (  # noqa: E402
    DecoderStats,
    LogicalFailureClassifier,
    binomial_rate_error_bar,
    bposd_decoder,
    classify_attempt,
    classify_verified_attempt,
    decode_with_bposd,
    numpy_rng,
    plot_logical_error_rates,
    sample_data_error,
    to_numpy_uint8,
)
from isomorphism import x, y  # noqa: E402
from unitary_decouple_decoder import DecouplingDecoder  # noqa: E402
from unitary_decouple_decoder.benchmarking import build_hz_dagger_matrix  # noqa: E402


DEFAULT_DISTANCES = "4,6,8"
DEFAULT_P_LIST = "0.01,0.02,0.03,0.05,0.08,0.1"
DEFAULT_OUTPUT = "results/unitary_decouple_decoder/color_666_unitary_vs_bposd.csv"
DEFAULT_PLOT_OUTPUT = "results/unitary_decouple_decoder/color_666_unitary_vs_bposd.png"


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distances", default=DEFAULT_DISTANCES)
    parser.add_argument("--p-list", default=DEFAULT_P_LIST)
    parser.add_argument("--shots", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--bp-method", choices=("minimum_sum", "product_sum"), default="minimum_sum")
    parser.add_argument("--ms-scaling-factor", type=float, default=0.625)
    parser.add_argument("--osd-order", type=int, default=2)
    parser.add_argument("--max-iter", type=int, default=288)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--plot-output", default=DEFAULT_PLOT_OUTPUT)
    parser.add_argument("--title", default="6.6.6 Color Code: Unitary Decouple vs BPOSD")
    parser.add_argument("--progress-interval", type=int, default=10000)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args):
    """Validate benchmark arguments."""
    if args.shots < 1:
        raise ValueError("--shots must be positive.")
    if args.max_iter < 1:
        raise ValueError("--max-iter must be positive.")
    if args.osd_order < 0:
        raise ValueError("--osd-order must be nonnegative.")
    if args.progress_interval < 1:
        raise ValueError("--progress-interval must be positive.")


def parse_distances(text):
    """Parse comma-separated coarse toric periods."""
    distances = []
    for item in text.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        value = int(stripped)
        if value < 1:
            raise ValueError("Distances must be positive.")
        distances.append(value)
    if not distances:
        raise ValueError("At least one distance is required.")
    return distances


def parse_probabilities(text):
    """Parse comma-separated physical error rates."""
    probabilities = []
    for item in text.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        value = float(stripped)
        if value < 0.0 or value > 1.0:
            raise ValueError("Physical error rates must be between 0 and 1.")
        probabilities.append(value)
    if not probabilities:
        raise ValueError("At least one physical error rate is required.")
    return probabilities


def color_666_decoder(distance):
    """Construct the compatible finite 6.6.6 color-code unitary decoder."""
    return DecouplingDecoder.from_bb(
        1 + x + x * y,
        1 + y + x * y,
        int(distance),
        int(distance),
    )


def decode_with_unitary(decoder, syndrome):
    """Decode with the unitary-decouple decoder."""
    return decoder.decode(syndrome, verify=False)


def run_one_probability(distance, decoder, classifier, bposd, p, args):
    """Run both decoders for one distance and physical error rate."""
    p = float(p)
    distance = int(distance)
    seed_offset = int(round(p * float(1_000_000)))
    rng = numpy_rng(int(args.seed) + 1009 * distance + seed_offset)
    stats = {
        "unitary_decouple": DecoderStats("unitary_decouple"),
        "bposd": DecoderStats("bposd"),
    }

    while stats["unitary_decouple"].shots < int(args.shots):
        true_error = sample_data_error(decoder.correction_size, p, rng)
        syndrome = decoder.syndrome(true_error)

        start = time.perf_counter()
        try:
            unitary_correction = decode_with_unitary(decoder, syndrome)
            elapsed = time.perf_counter() - start
            logical_failure, decode_failure = classify_verified_attempt(
                decoder.hx_source_finite,
                classifier,
                true_error,
                syndrome,
                unitary_correction,
                decoder_name="unitary_decouple",
            )
        except ValueError as exc:
            elapsed = time.perf_counter() - start
            shot = stats["unitary_decouple"].shots + 1
            raise ValueError(
                f"unitary_decouple invariant failed at d={distance} p={p:.6g} shot={shot}"
            ) from exc
        stats["unitary_decouple"].record(logical_failure, decode_failure, elapsed)

        start = time.perf_counter()
        try:
            bposd_correction = decode_with_bposd(bposd, syndrome)
            elapsed = time.perf_counter() - start
            logical_failure, decode_failure = classify_attempt(
                decoder.hx_source_finite,
                classifier,
                true_error,
                syndrome,
                bposd_correction,
            )
        except ValueError:
            elapsed = time.perf_counter() - start
            logical_failure, decode_failure = True, True
        stats["bposd"].record(logical_failure, decode_failure, elapsed)

        if stats["unitary_decouple"].shots % int(args.progress_interval) == 0:
            print(
                f"d={distance} p={p:.4g} shots={stats['unitary_decouple'].shots} "
                f"unitary={stats['unitary_decouple'].logical_failures} "
                f"bposd={stats['bposd'].logical_failures}",
                flush=True,
            )
    return stats


def stats_rows(distance, decoder, p, stats):
    """Convert per-decoder stats into CSV rows."""
    rows = []
    for decoder_name in ("unitary_decouple", "bposd"):
        item = stats[decoder_name]
        error_bar = binomial_rate_error_bar(item.logical_failures, item.shots)
        rows.append(
            {
                "distance": int(distance),
                "n": int(decoder.correction_size),
                "p": float(p),
                "decoder": item.decoder,
                "shots": item.shots,
                "logical_failures": item.logical_failures,
                "decode_failures": item.decode_failures,
                "logical_error_rate": item.logical_error_rate,
                "standard_error": item.standard_error,
                "logical_error_bar_low": error_bar.interval_low,
                "logical_error_bar_high": error_bar.interval_high,
                "elapsed_seconds": item.elapsed_seconds,
            }
        )
    return rows


def rows_with_error_bars(rows):
    """Return rows with Wilson error-bar bounds populated."""
    out = []
    for row in rows:
        item = dict(row)
        if "logical_error_bar_low" not in item or "logical_error_bar_high" not in item:
            error_bar = binomial_rate_error_bar(
                item.get("logical_failures", 0),
                item.get("shots", 0),
            )
            item.setdefault("logical_error_bar_low", error_bar.interval_low)
            item.setdefault("logical_error_bar_high", error_bar.interval_high)
        out.append(item)
    return out


def write_outputs(rows, output_path, title, plot_path=None):
    """Write CSV, Markdown summary, and PNG plot files."""
    rows = rows_with_error_bars(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plot_path is None:
        plot_path = output_path.with_suffix(".png")
    else:
        plot_path = Path(plot_path)
    fieldnames = [
        "distance",
        "n",
        "p",
        "decoder",
        "shots",
        "logical_failures",
        "decode_failures",
        "logical_error_rate",
        "standard_error",
        "logical_error_bar_low",
        "logical_error_bar_high",
        "elapsed_seconds",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    written_plot_path = plot_logical_error_rates(
        rows,
        plot_path,
        title,
        decoder_markers={"unitary_decouple": "o", "bposd": "s"},
    )
    md_path = output_path.with_suffix(".md")
    with md_path.open("w") as handle:
        handle.write(f"# {title}\n\n")
        handle.write(f"![Logical error rate comparison]({written_plot_path.name})\n\n")
        handle.write("| distance | n | p | decoder | shots | failures | decode failures | rate | stderr | err low | err high | seconds |\n")
        handle.write("|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            handle.write(
                f"| {row['distance']} | {row['n']} | {float(row['p']):.5g} | "
                f"{row['decoder']} | {row['shots']} | {row['logical_failures']} | "
                f"{row['decode_failures']} | {row['logical_error_rate']:.4f} | "
                f"{row['standard_error']:.4f} | {row['logical_error_bar_low']:.4f} | "
                f"{row['logical_error_bar_high']:.4f} | {row['elapsed_seconds']:.3f} |\n"
            )
    return md_path, written_plot_path


def run_smoke(args):
    """Run a short deterministic smoke check."""
    smoke_args = argparse.Namespace(
        shots=5,
        seed=args.seed,
        max_iter=min(int(args.max_iter), 20),
        osd_order=min(int(args.osd_order), 1),
        bp_method=args.bp_method,
        ms_scaling_factor=args.ms_scaling_factor,
        progress_interval=args.progress_interval,
    )
    decoder = color_666_decoder(2)
    classifier = LogicalFailureClassifier(
        decoder.hx_source_finite,
        build_hz_dagger_matrix(decoder),
    )
    bposd = bposd_decoder(
        to_numpy_uint8(decoder.hx_source_finite),
        0.05,
        smoke_args.max_iter,
        smoke_args.osd_order,
        bp_method=smoke_args.bp_method,
        ms_scaling_factor=smoke_args.ms_scaling_factor,
    )
    stats = run_one_probability(2, decoder, classifier, bposd, 0.05, smoke_args)
    for item in stats.values():
        if item.shots != 5:
            raise AssertionError("Smoke benchmark did not run the requested shots.")
    if stats["unitary_decouple"].decode_failures != 0:
        raise AssertionError("Unitary decoder smoke run recorded decode failures.")
    print("smoke: ok")


def run_benchmark(args):
    """Run the configured benchmark and return output rows."""
    rows = []
    distances = parse_distances(args.distances)
    probabilities = parse_probabilities(args.p_list)
    for distance in distances:
        print(f"constructing d={distance} unitary decoder", flush=True)
        decoder = color_666_decoder(distance)
        classifier = LogicalFailureClassifier(
            decoder.hx_source_finite,
            build_hz_dagger_matrix(decoder),
        )
        hx_numpy = to_numpy_uint8(decoder.hx_source_finite)
        print(f"d={distance} n={decoder.correction_size}", flush=True)
        for p in probabilities:
            print(f"running d={distance} p={p}", flush=True)
            bposd = bposd_decoder(
                hx_numpy,
                p,
                args.max_iter,
                args.osd_order,
                bp_method=args.bp_method,
                ms_scaling_factor=args.ms_scaling_factor,
            )
            stats = run_one_probability(distance, decoder, classifier, bposd, p, args)
            rows.extend(stats_rows(distance, decoder, p, stats))
    return rows


def main():
    """Run the benchmark."""
    args = parse_args()
    if args.smoke:
        run_smoke(args)
        return
    rows = run_benchmark(args)
    output_path = Path(args.output)
    md_path, plot_path = write_outputs(rows, output_path, args.title, args.plot_output)
    print(f"wrote {output_path}")
    print(f"wrote {md_path}")
    print(f"wrote {plot_path}")


if __name__ == "__main__":
    main()
