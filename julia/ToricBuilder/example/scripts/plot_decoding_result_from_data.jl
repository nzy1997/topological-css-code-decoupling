using JSON
using CairoMakie
using LaTeXStrings

const TORICBUILDER_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPOSITORY_ROOT = normpath(joinpath(TORICBUILDER_ROOT, "..", ".."))
const EXPORTED_DATA_PATH = joinpath(TORICBUILDER_ROOT, "results", "decoding_benchmark.json")
const OUTPUT_PDF_PATH = joinpath(REPOSITORY_ROOT, "build", "reproduction", "decoding_benchmark.pdf")
const LOGICAL_ERROR_PANEL_LIMITS = (x = (1.0e-4, 0.25), y = (1.0e-8, 1.0))

distance_color(cols, d::Integer) = cols[d - 3]
distance_legend_label(d::Integer) = latexstring("d_{\\text{TC}} = ", d)

function logical_plot_series(series)
    interval = series["interval"]
    estimate = Float64.(interval["estimate"])
    low = Float64.(interval["low"])
    high = Float64.(interval["high"])
    likelihood_ratio = Float64(interval["likelihood_ratio"])
    center = map(estimate, high) do point, upper
        iszero(point) ? upper / likelihood_ratio : point
    end
    lower = map(estimate, center, low) do point, plotted, lower_bound
        iszero(point) ? 0.0 : plotted - lower_bound
    end
    upper = high .- center
    return (; center, lower, upper)
end

function draw_decoding_result_from_data(data_path::AbstractString=EXPORTED_DATA_PATH)
    payload = JSON.parsefile(data_path)
    ds = Int.(payload["distances"])
    cols = cgrad(:tab10, 20; categorical=true)

    fig = Figure(size = (1200, 400))
    ax = Axis(
        fig[1, 1],
        xscale = log10,
        yscale = log10,
        xlabel = "Physical error rate",
        ylabel = "Logical error rate",
    )
    ax2 = Axis(
        fig[1, 2],
        xscale = log10,
        yscale = log10,
        xlabel = "Physical error rate",
        ylabel = "Decoding time per sample(s)",
    )

    for d in ds
        key = string(d)
        color = distance_color(cols, d)

        bposd = payload["logical_error"][key]["bposd"]
        bposd_plot = logical_plot_series(bposd)
        errorbars!(ax, bposd["pvec"], bposd_plot.center, bposd_plot.lower, bposd_plot.upper; whiskerwidth=10, color=color)
        scatterlines!(ax, bposd["pvec"], bposd_plot.center; color=color, linestyle=:dash)

        unitary = payload["logical_error"][key]["unitary"]
        unitary_plot = logical_plot_series(unitary)
        errorbars!(ax, unitary["pvec"], unitary_plot.center, unitary_plot.lower, unitary_plot.upper; whiskerwidth=10, color=color)
        scatterlines!(ax, unitary["pvec"], unitary_plot.center; color=color)
    end

    for d in ds
        key = string(d)
        color = distance_color(cols, d)

        bposd = payload["decoding_time"][key]["bposd"]
        scatterlines!(ax2, bposd["pvec"], bposd["time_res"]; color=color, linestyle=:dash)

        unitary = payload["decoding_time"][key]["unitary"]
        scatterlines!(ax2, unitary["pvec"], unitary["time_res"]; color=color)
    end

    xlims!(ax, LOGICAL_ERROR_PANEL_LIMITS.x...)
    ylims!(ax, LOGICAL_ERROR_PANEL_LIMITS.y...)

    legend_elements = Any[
        LineElement(color=:black, linestyle=:solid, linewidth=3),
        LineElement(color=:black, linestyle=:dash, linewidth=3),
    ]
    legend_labels = Any["Unitary-decouple-based decoder", "BP-OSD"]
    for d in ds
        push!(legend_elements, LineElement(color=distance_color(cols, d), linewidth=3))
        push!(legend_labels, distance_legend_label(d))
    end

    fig[1, 3] = Legend(fig, legend_elements, legend_labels; framevisible=false)
    return fig
end

function save_decoding_result_from_data(
    data_path::AbstractString=EXPORTED_DATA_PATH,
    output_path::AbstractString=OUTPUT_PDF_PATH,
)
    mkpath(dirname(output_path))
    fig = draw_decoding_result_from_data(data_path)
    save(output_path, fig)
    return output_path
end


if abspath(PROGRAM_FILE) == @__FILE__
    save_decoding_result_from_data()
end
