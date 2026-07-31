"""Compare the paper's BB-family decoder with BP-OSD."""

from pathlib import Path
import argparse
import csv
import sys
import time


SOURCE_ROOT = Path(__file__).resolve().parents[2]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from decoder_core import (  # noqa: E402
    DecoderStats,
    bp_osd_decoder,
    classify_attempt,
    decode_with_bp_osd,
    numpy_rng,
    plot_logical_error_rates,
    sample_data_error,
    to_numpy_uint8,
)
from decoder_core.benchmark_metadata import runtime_metadata  # noqa: E402
from decoupling import x, y  # noqa: E402
from unitary_decouple_based_decoder import UnitaryDecoupleBasedDecoder  # noqa: E402


DEFAULT_DISTANCES = "4,6"
DEFAULT_PROBABILITIES = "0.01,0.03,0.05,0.07,0.10"
DEFAULT_OUTPUT = (
    SOURCE_ROOT
    / "results"
    / "unitary_decouple_based_decoder"
    / "paper_bb_family_unitary_decouple_based_vs_bp_osd.csv"
)
DECODER_LABEL = "unitary_decouple_based"
BPOSD_LABEL = "bp_osd"
F = 1 + x + x**-1 * y
G = 1 + y + x * y


def parse_args():
    """Parse publication defaults and optional smoke overrides."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distances", default=DEFAULT_DISTANCES)
    parser.add_argument("--probabilities", default=DEFAULT_PROBABILITIES)
    parser.add_argument("--shots", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--bp-method", default="minimum_sum")
    parser.add_argument(
        "--bp-schedule",
        choices=("serial", "parallel"),
        default="serial",
    )
    parser.add_argument("--ms-scaling-factor", type=float, default=0.625)
    parser.add_argument("--osd-order", type=int, default=2)
    parser.add_argument("--max-iter", type=int, default=288)
    parser.add_argument("--progress-interval", type=int, default=10000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.shots <= 0:
        raise ValueError("--shots must be positive.")
    if args.max_iter <= 0:
        raise ValueError("--max-iter must be positive.")
    if args.osd_order < 0:
        raise ValueError("--osd-order must be nonnegative.")
    if args.progress_interval <= 0:
        raise ValueError("--progress-interval must be positive.")
    return args


def parse_positive_integers(text):
    """Parse comma-separated positive integers."""
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values or any(value <= 0 for value in values):
        raise ValueError("At least one positive distance is required.")
    return values


def parse_probabilities(text):
    """Parse comma-separated probabilities in the closed unit interval."""
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values or any(value < 0 or value > 1 for value in values):
        raise ValueError("Probabilities must lie between zero and one.")
    return values


def build_decoder(distance):
    """Construct the unitary-decouple-based decoder for one torus distance."""
    return UnitaryDecoupleBasedDecoder.from_bb(
        F,
        G,
        distance,
        distance,
        max_period=7,
    )


def run_probability(
    decoder,
    classifier,
    bp_osd,
    distance,
    probability,
    rng,
    args,
):
    """Run paired Monte Carlo trials using an identical sampled error stream."""
    stats = {
        DECODER_LABEL: DecoderStats(DECODER_LABEL),
        BPOSD_LABEL: DecoderStats(BPOSD_LABEL),
    }
    while stats[DECODER_LABEL].shots < args.shots:
        error = sample_data_error(decoder.correction_size, probability, rng)
        syndrome = decoder.syndrome(error)

        started = time.perf_counter()
        correction = decoder.decode(syndrome, verify=False)
        elapsed = time.perf_counter() - started
        if decoder.syndrome(correction) != syndrome:
            raise AssertionError(
                "The unitary-decouple-based correction has the wrong syndrome."
            )
        logical_failure, decode_failure = classify_attempt(
            decoder.h_x_input_finite,
            classifier,
            error,
            syndrome,
            correction,
        )
        stats[DECODER_LABEL].record(
            logical_failure,
            decode_failure,
            elapsed,
        )

        started = time.perf_counter()
        try:
            correction = decode_with_bp_osd(bp_osd, syndrome)
            elapsed = time.perf_counter() - started
        except ValueError:
            elapsed = time.perf_counter() - started
            logical_failure, decode_failure = True, True
        else:
            logical_failure, decode_failure = classify_attempt(
                decoder.h_x_input_finite,
                classifier,
                error,
                syndrome,
                correction,
            )
        stats[BPOSD_LABEL].record(
            logical_failure,
            decode_failure,
            elapsed,
        )
        if stats[DECODER_LABEL].shots % args.progress_interval == 0:
            print(
                f"d={distance} p={probability:.2f} "
                f"shots={stats[DECODER_LABEL].shots}",
                flush=True,
            )
    return stats


def result_rows(decoder, distance, probability, stats, args):
    """Convert paired statistics to canonical CSV rows."""
    rows = []
    for label in (DECODER_LABEL, BPOSD_LABEL):
        item = stats[label]
        rows.append(
            {
                "family": "paper_bb",
                "f": str(F),
                "g": str(G),
                "distance": distance,
                "n": decoder.correction_size,
                "p": probability,
                "decoder": label,
                "shots": item.shots,
                "seed": args.seed,
                "bp_method": args.bp_method,
                "bp_schedule": args.bp_schedule,
                "ms_scaling_factor": args.ms_scaling_factor,
                "osd_order": args.osd_order,
                "max_iter": args.max_iter,
                "logical_failures": item.logical_failures,
                "decode_failures": item.decode_failures,
                "logical_error_rate": item.logical_error_rate,
                "standard_error": item.standard_error,
                "elapsed_seconds": item.elapsed_seconds,
            }
        )
    return rows


def write_outputs(rows, args):
    """Write the canonical CSV, Markdown table, and comparison plot."""
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = runtime_metadata(
        adapter=DECODER_LABEL,
        suite="paper_bb_family_vs_bp_osd",
        source_root=SOURCE_ROOT,
    )
    metadata.pop("source_root", None)
    for row in rows:
        row.update(metadata)
    columns = [
        "family",
        "f",
        "g",
        "distance",
        "n",
        "p",
        "decoder",
        "shots",
        "seed",
        "bp_method",
        "bp_schedule",
        "ms_scaling_factor",
        "osd_order",
        "max_iter",
        "logical_failures",
        "decode_failures",
        "logical_error_rate",
        "standard_error",
        "elapsed_seconds",
        "adapter",
        "suite",
        "source_revision",
        "source_tree_dirty",
        "sage_version",
        "python_version",
        "platform",
    ]
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    plot_path = output.with_suffix(".png")
    plot_logical_error_rates(
        rows,
        plot_path,
        "Paper BB Family: Unitary-Decouple-Based vs BP-OSD",
        decoder_markers={DECODER_LABEL: "o", BPOSD_LABEL: "s"},
    )
    markdown = output.with_suffix(".md")
    with markdown.open("w", encoding="utf-8") as stream:
        stream.write("# Paper BB Family: Unitary-Decouple-Based vs BP-OSD\n\n")
        stream.write(f"- H_X: `({F}, {G})`\n")
        stream.write(f"- Seed: `{args.seed}`\n")
        stream.write(f"- Shots per point: `{args.shots}`\n")
        stream.write(f"- BP method: `{args.bp_method}`\n")
        stream.write(f"- BP schedule: `{args.bp_schedule}`\n")
        stream.write(f"- Min-sum scaling: `{args.ms_scaling_factor}`\n")
        stream.write(f"- OSD order: `{args.osd_order}`\n")
        stream.write(f"- Maximum BP iterations: `{args.max_iter}`\n")
        stream.write(f"- Sage: `{metadata['sage_version']}`\n")
        stream.write(f"- Python: `{metadata['python_version']}`\n")
        stream.write(f"- Platform: `{metadata['platform']}`\n")
        stream.write(f"- Source revision: `{metadata['source_revision']}`\n")
        stream.write(f"- Source tree dirty: `{metadata['source_tree_dirty']}`\n\n")
        stream.write(f"![Logical error comparison]({plot_path.name})\n\n")
        stream.write(
            "| distance | n | p | decoder | shots | logical failures | "
            "decode failures | logical rate | stderr | seconds |\n"
        )
        stream.write("| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in rows:
            stream.write(
                f"| {row['distance']} | {row['n']} | {row['p']:.2f} | "
                f"{row['decoder']} | {row['shots']} | {row['logical_failures']} | "
                f"{row['decode_failures']} | {row['logical_error_rate']:.8f} | "
                f"{row['standard_error']:.8f} | {row['elapsed_seconds']:.6f} |\n"
            )
    return output, markdown, plot_path


def run(args):
    """Run all requested distances and probabilities."""
    distances = parse_positive_integers(args.distances)
    probabilities = parse_probabilities(args.probabilities)
    rows = []
    rng = numpy_rng(args.seed)
    for distance in distances:
        decoder = build_decoder(distance)
        classifier = decoder.logical_failure_classifier
        h_x_numpy = to_numpy_uint8(decoder.h_x_input_finite)
        for probability in probabilities:
            bp_osd = bp_osd_decoder(
                h_x_numpy,
                probability,
                args.max_iter,
                args.osd_order,
                bp_method=args.bp_method,
                ms_scaling_factor=args.ms_scaling_factor,
                schedule=args.bp_schedule,
            )
            stats = run_probability(
                decoder,
                classifier,
                bp_osd,
                distance,
                probability,
                rng,
                args,
            )
            rows.extend(
                result_rows(decoder, distance, probability, stats, args)
            )
    return rows


def main():
    """Run either the publication benchmark or a deterministic smoke check."""
    args = parse_args()
    if args.smoke:
        args.distances = "4"
        args.probabilities = "0.05"
        args.shots = 5
        args.progress_interval = 5
        args.output = Path("/tmp/paper_bb_family_smoke.csv")
    rows = run(args)
    outputs = write_outputs(rows, args)
    if args.smoke:
        if len(rows) != 2 or any(row["shots"] != 5 for row in rows):
            raise AssertionError("The smoke benchmark did not finish all paired trials.")
        print("paper BB benchmark smoke: ok")
    else:
        for output in outputs:
            print(f"wrote {output}")


if __name__ == "__main__":
    main()
