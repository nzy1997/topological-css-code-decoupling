using JSON
using CairoMakie
using LaTeXStrings

const TORICBUILDER_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPOSITORY_ROOT = normpath(joinpath(TORICBUILDER_ROOT, "..", ".."))
const EXPORTED_DATA_PATH = joinpath(TORICBUILDER_ROOT, "results", "decoding_benchmark.json")
const OUTPUT_PDF_PATH = joinpath(REPOSITORY_ROOT, "build", "reproduction", "decoding_benchmark.pdf")

distance_color(cols, d::Integer) = cols[d - 3]
distance_legend_label(d::Integer) = latexstring("d_{\\text{TC}} = ", d)

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
        errorbars!(ax, bposd["pvec"], bposd["fit"]["avs_plot"], bposd["fit"]["ylow"], bposd["fit"]["yhigh"]; whiskerwidth=10, color=color)
        scatterlines!(ax, bposd["pvec"], bposd["fit"]["avs_plot"]; color=color, linestyle=:dash)

        unitary_decouple = payload["logical_error"][key]["unitary_decouple"]
        errorbars!(ax, unitary_decouple["pvec"], unitary_decouple["fit"]["avs_plot"], unitary_decouple["fit"]["ylow"], unitary_decouple["fit"]["yhigh"]; whiskerwidth=10, color=color)
        scatterlines!(ax, unitary_decouple["pvec"], unitary_decouple["fit"]["avs_plot"]; color=color)
    end

    for d in ds
        key = string(d)
        color = distance_color(cols, d)

        bposd = payload["decoding_time"][key]["bposd"]
        scatterlines!(ax2, bposd["pvec"], bposd["time_res"]; color=color, linestyle=:dash)

        unitary_decouple = payload["decoding_time"][key]["unitary_decouple"]
        scatterlines!(ax2, unitary_decouple["pvec"], unitary_decouple["time_res"]; color=color)
    end

    limits = payload["logical_error_panel_limits"]
    xlims!(ax, limits["x"][1], limits["x"][2])
    ylims!(ax, limits["y"][1], limits["y"][2])

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
