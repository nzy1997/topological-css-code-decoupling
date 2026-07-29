using CairoMakie
using GLM
using DataFrames

function default_input_paths()
    results_dir = joinpath(dirname(dirname(@__DIR__)), "results")
    return [
        joinpath(results_dir, "additional_bb_code_decoupling_results.md"),
        joinpath(results_dir, "bb_code_decoupling_results.md"),
    ]
end

function parse_markdown_table(filename)
    lines = readlines(filename)
    data_lines = lines[3:end]

    areas = Float64[]
    times = Float64[]
    degrees = Float64[]

    for line in data_lines
        parts = split(line, "|")

        if length(parts) < 13
            continue
        end

        area_str = strip(parts[7])
        time_str = strip(parts[11])
        max_term_q = strip(parts[12])
        degree_str = strip(parts[13])

        if max_term_q == "-"
            continue
        end

        if area_str == "-" || time_str == "-" || degree_str == "-"
            continue
        end

        try
            area = parse(Float64, area_str)
            time = parse(Float64, time_str)
            degree = parse(Float64, degree_str)
            push!(areas, area)
            push!(times, time)
            push!(degrees, degree)
        catch
            continue
        end
    end

    return (; areas, times, degrees)
end

function parse_markdown_tables(filenames)
    areas = Float64[]
    times = Float64[]
    degrees = Float64[]

    for filename in filenames
        data = parse_markdown_table(filename)
        append!(areas, data.areas)
        append!(times, data.times)
        append!(degrees, data.degrees)
    end

    return (; qs=2 .* areas, areas, times, degrees)
end

function positive_pairs(xs, ys)
    keep = [isfinite(x) && isfinite(y) && x > 0 && y > 0 for (x, y) in zip(xs, ys)]
    return (; xs=xs[keep], ys=ys[keep])
end

function fit_power_law(x, y)
    df = DataFrame(log_x=log10.(x), log_y=log10.(y))
    model = lm(@formula(log_y ~ log_x), df)
    return coef(model)[2], coef(model)[1]
end

function fitted_curve(x, slope, intercept)
    x_range = range(minimum(x), maximum(x), length=100)
    y_fit = 10 .^ (intercept .+ slope .* log10.(x_range))
    return x_range, y_fit
end

function main(;
    input_paths=default_input_paths(),
    input_path=nothing,
    output_path=joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "build", "reproduction", "area_comparison_plot.pdf"),
    display_figure::Bool=false,
)
    if input_path !== nothing
        input_paths = [input_path]
    end

    data = parse_markdown_tables(input_paths)
    qs = data.qs
    times = data.times
    degrees = data.degrees
    println("Extracted $(length(qs)) data points from $(length(input_paths)) result table(s)")
    time_data = positive_pairs(qs, times)
    degree_data = positive_pairs(qs, degrees)
    println("Using $(length(time_data.xs)) positive time points and $(length(degree_data.xs)) positive degree points")
    mkpath(dirname(output_path))

    fig = Figure(size=(700, 300))
    q_ticks = ([1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0], ["1", "10", "100", "1000", "10000", "100000"])

    ax_time = Axis(
        fig[1, 1],
        xlabel=L"$q$",
        ylabel=L"$T$ (s)",
        xscale=log10,
        yscale=log10,
        title=L"(a) Running time $T$ vs $q$",
        xticks=q_ticks,
        yticks=([0.1, 10.0, 1000.0], ["0.1", "10", "1000"]),
    )

    scatter!(ax_time, time_data.xs, time_data.ys, color=:blue, markersize=10, label=L"Data$$")

    time_slope, time_intercept = fit_power_law(time_data.xs, time_data.ys)
    time_x, time_y = fitted_curve(time_data.xs, time_slope, time_intercept)
    lines!(
        ax_time,
        time_x,
        time_y,
        color=:red,
        linewidth=2,
        label=L"Fit:$T$ \propto $q^{%$(round(time_slope, digits=2))}$",
    )
    axislegend(ax_time, position=:lt)
    text!(
        ax_time,
        0.02,
        0.98;
        text="(a)",
        space=:relative,
        align=(:left, :top),
        fontsize=24,
        color=:black,
    )

    println("Fit result: log10(time) = $(time_slope) * log10(q) + $(time_intercept)")
    println("That is: time = 10^$(time_intercept) * q^$(time_slope)")

    ax_degree = Axis(
        fig[1, 2],
        xlabel=L"$q$",
        ylabel=L"$\deg(\psi_1^{-1})$",
        xscale=log10,
        yscale=log10,
        title=L"(b) Locality deg$(\psi_1^{-1})$ vs $q$",
        xticks=q_ticks,
        yticks=([10.0, 100.0], ["10", "100"]),
    )

    scatter!(ax_degree, degree_data.xs, degree_data.ys, color=:blue, markersize=10, label=L"Data$$")

    degree_slope, degree_intercept = fit_power_law(degree_data.xs, degree_data.ys)
    degree_x, degree_y = fitted_curve(degree_data.xs, degree_slope, degree_intercept)
    lines!(
        ax_degree,
        degree_x,
        degree_y,
        color=:red,
        linewidth=2,
        label=L"Fit:deg$(\psi_1^{-1})$ \propto $q^{%$(round(degree_slope, digits=2))}$",
    )
    axislegend(ax_degree, position=:lt)
    text!(
        ax_degree,
        0.02,
        0.98;
        text="(b)",
        space=:relative,
        align=(:left, :top),
        fontsize=24,
        color=:black,
    )

    println("Fit result: log10(deg(psi_1^-1)) = $(degree_slope) * log10(q) + $(degree_intercept)")
    println("That is: deg(psi_1^-1) = 10^$(degree_intercept) * q^$(degree_slope)")

    save(output_path, fig)
    println("Saved figure to $(output_path)")

    if display_figure
        display(fig)
    end

    return (; fig, qs, areas=data.areas, times, degrees, input_paths, output_path, time_data, degree_data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
