"""Shared plotting helpers for decoder benchmark result tables."""

from pathlib import Path

from .stats import binomial_rate_error_bar, binomial_wilson_interval


DEFAULT_DECODER_MARKERS = {
    "unitary_decouple": "o",
    "bposd": "s",
    "matching": "^",
    "bp_reweighted": "D",
}

DEFAULT_DECODER_LINESTYLES = {
    "unitary_decouple": "-",
    "bposd": "--",
    "matching": "-.",
    "bp_reweighted": ":",
}


def plot_rate_points(decoder_rows, *, interval_z=1.0):
    """Return log-safe plot coordinates and asymmetric error bars.

    Args:
        decoder_rows: Benchmark rows associated with one decoder.
        interval_z: Normal-approximation z value for Wilson intervals.

    Returns:
        object: Statistical estimate or interval data for benchmark reporting.
    """
    # Keep syndrome normalization separate from the matching solve.
    probabilities = []
    rates = []
    zero_rate = []
    lower_errors = []
    upper_errors = []
    for row in decoder_rows:
        shots = int(row["shots"])
        rate = float(row["logical_error_rate"])
        failures = int(row.get("logical_failures", round(rate * shots)))
        error_bar = binomial_rate_error_bar(
            failures,
            shots,
            z=interval_z,
        )
        probabilities.append(float(row["p"]))
        rates.append(error_bar.display_rate)
        zero_rate.append(error_bar.zero_rate)
        lower_errors.append(error_bar.lower_error)
        upper_errors.append(error_bar.upper_error)
    return probabilities, rates, zero_rate, (lower_errors, upper_errors)


def plot_logical_error_rates(
    rows,
    plot_path,
    title,
    *,
    distance_key="distance",
    decoder_key="decoder",
    decoder_markers=None,
    decoder_linestyles=None,
    interval_z=1.0,
):
    """Write a log-log logical-error-rate plot with binomial error bars.

    Args:
        rows: Rows or benchmark records to serialize.
        plot_path: Output path for the generated benchmark plot.
        title: Plot or report title.
        distance_key: Dictionary key that identifies code distance in benchmark rows.
        decoder_key: Dictionary key that identifies the decoder column in benchmark rows.
        decoder_markers: Marker mapping used when plotting decoder curves.
        decoder_linestyles: Line-style mapping used when plotting decoder curves.
        interval_z: Normal-approximation z value for Wilson intervals.

    Returns:
        object: Statistical estimate or interval data for benchmark reporting.
    """
    # Keep syndrome normalization separate from the matching solve.
    plot_path = Path(plot_path)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    markers = dict(DEFAULT_DECODER_MARKERS)
    if decoder_markers:
        markers.update(decoder_markers)
    linestyles = dict(DEFAULT_DECODER_LINESTYLES)
    if decoder_linestyles:
        linestyles.update(decoder_linestyles)
    distances = sorted({int(row[distance_key]) for row in rows})
    decoders = sorted({str(row[decoder_key]) for row in rows})
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    distance_colors = {
        distance: color_cycle[index % len(color_cycle)]
        for index, distance in enumerate(distances)
    }

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for distance in distances:
        color = distance_colors[distance]
        for decoder_name in decoders:
            decoder_rows = sorted(
                [
                    row for row in rows
                    if int(row[distance_key]) == distance
                    and str(row[decoder_key]) == decoder_name
                ],
                key=lambda row: row["p"],
            )
            if not decoder_rows:
                continue
            probabilities, rates, zero_rate, yerr = plot_rate_points(
                decoder_rows,
                interval_z=interval_z,
            )
            marker = markers.get(decoder_name, "o")
            linestyle = linestyles.get(decoder_name, "-")
            ax.errorbar(
                probabilities,
                rates,
                yerr=yerr,
                color=color,
                marker=marker,
                linestyle=linestyle,
                markersize=5.5,
                linewidth=1.6,
                elinewidth=1.2,
                capsize=4,
                capthick=1.2,
                label=f"d={distance} {decoder_name}",
            )
            zero_probabilities = [
                p for p, is_zero in zip(probabilities, zero_rate) if is_zero
            ]
            zero_rates = [r for r, is_zero in zip(rates, zero_rate) if is_zero]
            if zero_probabilities:
                ax.scatter(
                    zero_probabilities,
                    zero_rates,
                    marker=marker,
                    facecolors="white",
                    edgecolors=color,
                    zorder=3,
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Physical error rate")
    ax.set_ylabel("Logical error rate (0 shown at 0.5/shots)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    return plot_path
