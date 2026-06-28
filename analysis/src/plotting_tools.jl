function save_histogram(data, field, filename; bins=50, title="", xlabel="",
                         color=nothing)
    c = color === nothing ? (:steelblue, 0.7) : color
    values = [getindex(d, field) for d in data]
    finite_vals = [v for v in values if isfinite(v)]
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel="Count")
    hist!(ax, finite_vals, bins=bins, color=c)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved plot: $filename")
    return f
end

function save_scatter(data, xfield, yfield, filename; max_points=10_000, title="",
                      xlabel="", ylabel="", color=nothing)
    c = color === nothing ? (:steelblue, 0.5) : color
    xvals = [getindex(d, xfield) for d in data]
    yvals = [getindex(d, yfield) for d in data]
    finite = isfinite.(xvals) .& isfinite.(yvals)
    xvals, yvals = xvals[finite], yvals[finite]
    n = length(xvals)
    if n > max_points
        idx = rand(1:n, max_points)
        xvals, yvals = xvals[idx], yvals[idx]
    end
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    scatter!(ax, xvals, yvals, markersize=3, color=c)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved plot: $filename ($(n) original points)")
    return f
end

function save_scatter_color(data, xfield, yfield, cfield, filename;
                            max_points=10_000, title="",
                            xlabel="", ylabel="", clabel="",
                            markersize=3, colormap=:viridis)
    xvals = [getindex(d, xfield) for d in data]
    yvals = [getindex(d, yfield) for d in data]
    cvals = [getindex(d, cfield) for d in data]
    finite = isfinite.(xvals) .& isfinite.(yvals) .& isfinite.(cvals)
    xvals, yvals, cvals = xvals[finite], yvals[finite], cvals[finite]
    n = length(xvals)
    if n > max_points
        idx = rand(1:n, max_points)
        xvals, yvals, cvals = xvals[idx], yvals[idx], cvals[idx]
    end
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    sc = scatter!(ax, xvals, yvals, color=cvals, markersize=markersize,
                  colormap=colormap, transparency=true)
    Colorbar(f[1, 2], sc, label=clabel)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved colored scatter: $filename")
    return f
end

function save_population_comparison(populations, field, filename; bins=50, title="", xlabel="")
    labels = String[]
    all_vals = Float64[]
    pop_labels = String[]
    for (label, data) in pairs(populations)
        vals = [getindex(d, field) for d in data]
        finite_vals = [v for v in vals if isfinite(v)]
        append!(all_vals, finite_vals)
        append!(pop_labels, fill(string(label), length(finite_vals)))
    end
    f = Figure(size=(1000, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel="Density")
    colors = [:steelblue, :orange, :crimson]
    for (i, (label, data)) in enumerate(pairs(populations))
        vals = [getindex(d, field) for d in data]
        finite_vals = [v for v in vals if isfinite(v)]
        if !isempty(finite_vals)
            hist!(ax, finite_vals, bins=bins, color=(colors[(i-1)%3+1], 0.5),
                  label=string(label))
        end
    end
    axislegend(ax)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved population comparison: $filename")
    return f
end

function save_scaling_law_plot(x, y, fit, filename; title="", xlabel="", ylabel="",
                                max_points=10_000, color=nothing)
    c = color === nothing ? (:steelblue, 0.5) : color
    finite = isfinite.(x) .& isfinite.(y)
    xf, yf = x[finite], y[finite]
    n = length(xf)
    if n > max_points
        idx = rand(1:n, max_points)
        xf, yf = xf[idx], yf[idx]
    end
    x_range = collect(range(minimum(xf), maximum(xf); length=100))
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    scatter!(ax, xf, yf, markersize=3, color=c)
    lines!(ax, x_range, fit.intercept .+ fit.slope .* x_range, color=:red, linewidth=2,
           label="R²=$(round(fit.r2, digits=3))")
    axislegend(ax)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved scaling law plot: $filename")
    return f
end

function save_heatmap(data, xfield, yfield, filename; nbins_x=25, nbins_y=25,
                      title="", xlabel="", ylabel="", colormap=:hot)
    xvals = [getindex(d, xfield) for d in data]
    yvals = [getindex(d, yfield) for d in data]
    bins_x = range(minimum(xvals), maximum(xvals); length=nbins_x+1)
    bins_y = range(minimum(yvals), maximum(yvals); length=nbins_y+1)
    counts = zeros(Int, nbins_y, nbins_x)
    for i in eachindex(xvals)
        xi = clamp(searchsortedfirst(bins_x, xvals[i]) - 1, 1, nbins_x)
        yi = clamp(searchsortedfirst(bins_y, yvals[i]) - 1, 1, nbins_y)
        counts[yi, xi] += 1
    end
    f = Figure(size=(800, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    hm = heatmap!(ax, bins_x, bins_y, counts, colormap=colormap)
    Colorbar(f[1, 2], hm)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved heatmap: $filename")
    return f
end

function save_boxplot(data, groupfield, valuefield, filename; title="",
                      xlabel="", ylabel="", max_groups=20)
    group_vals = [getindex(d, groupfield) for d in data]
    value_vals = [getindex(d, valuefield) for d in data]
    groups = unique(group_vals)
    if length(groups) > max_groups
        sorted_groups = sort(groups)
        keep = Set(sorted_groups[1:max_groups])
        mask = [g in keep for g in group_vals]
        group_vals = group_vals[mask]
        value_vals = value_vals[mask]
        groups = sorted_groups[1:max_groups]
    end
    sorted_groups = sort(groups)
    plot_data = [Float64[value_vals[i] for i in eachindex(group_vals) if group_vals[i] == g] for g in sorted_groups]
    f = Figure(size=(900, 600))
    ax = Axis(f[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    boxplot!(ax, [float(g) for g in sorted_groups], plot_data)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved boxplot: $filename")
    return f
end

function save_pca_plot(pca_result, data, filename; color_field=:log10_abundance,
                       title="PCA Projection", max_points=10_000)
    proj = pca_result.projections
    n = size(proj, 1)
    colors = [getindex(d, color_field) for d in data]
    if n > max_points
        idx = rand(1:n, max_points)
        proj = proj[idx, :]
        colors = colors[idx]
    end
    var_pc1 = round(100 * pca_result.explained_var[1], digits=1)
    var_pc2 = round(100 * pca_result.explained_var[2], digits=1)
    f = Figure(size=(900, 700))
    ax = Axis(f[1, 1], title=title,
              xlabel="PC1 ($var_pc1% variance)",
              ylabel="PC2 ($var_pc2% variance)")
    sc = scatter!(ax, proj[:, 1], proj[:, 2], color=colors,
                  markersize=4, colormap=:viridis, transparency=true)
    Colorbar(f[1, 2], sc, label=string(color_field))
    if length(pca_result.explained_var) >= 3
        var_pc3 = round(100 * pca_result.explained_var[3], digits=1)
        ax2 = Axis(f[2, 1], title="PC3 vs PC1",
                   xlabel="PC1 ($var_pc1% variance)",
                   ylabel="PC3 ($var_pc3% variance)")
        scatter!(ax2, proj[:, 1], proj[:, 3], color=colors,
                 markersize=4, colormap=:viridis, transparency=true)
    end
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved PCA plot: $filename")
    return f
end

function save_explained_variance(pca_result, filename; title="Explained Variance Ratio")
    f = Figure(size=(800, 500))
    ax = Axis(f[1, 1], title=title,
              xlabel="Principal Component", ylabel="Explained Variance Ratio")
    n = length(pca_result.explained_var)
    barplot!(ax, 1:n, pca_result.explained_var, color=:steelblue)
    lines!(ax, 1:n, pca_result.explained_var_cumul, color=:red, linewidth=2)
    ylims!(ax, low=0, high=1.05)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved explained variance plot: $filename")
    return f
end

function save_geometry_abundance_plot(geometry_data, axion_data, filename; max_points=10_000)
    geom_abundances = Float64[]
    n_light_list = Int[]
    h11_list = Int[]
    for g in geometry_data
        axions = [d for d in axion_data
                  if d.h11 == g.h11 && d.polytope == g.polytope && d.frst == g.frst]
        total_abund = isempty(axions) ? 0.0 : log10(sum(10 .^ [a.log10_abundance for a in axions]))
        push!(geom_abundances, total_abund)
        push!(n_light_list, g.n_light_axions)
        push!(h11_list, g.h11)
    end
    n = length(geom_abundances)
    if n > max_points
        idx = rand(1:n, max_points)
        geom_abundances = geom_abundances[idx]
        n_light_list = n_light_list[idx]
        h11_list = h11_list[idx]
    end
    f = Figure(size=(1200, 800))
    ax1 = Axis(f[1, 1], title="Total Geometry Abundance Distribution",
               xlabel="log10(Total Ω_a h²)", ylabel="Count")
    hist!(ax1, geom_abundances, bins=50, color=:steelblue, alpha=0.7)
    ax2 = Axis(f[1, 2], title="Light Axions vs Geometry Abundance",
               xlabel="n_light_axions", ylabel="log10(Total Ω_a h²)")
    scatter!(ax2, n_light_list, geom_abundances, markersize=3, color=:steelblue, alpha=0.5)
    ax3 = Axis(f[2, 1], title="Geometry Abundance vs h11",
               xlabel="h11", ylabel="log10(Total Ω_a h²)")
    scatter!(ax3, h11_list, geom_abundances, markersize=3, color=:steelblue, alpha=0.5)
    ax4 = Axis(f[2, 2], title="Light Axions per Geometry",
               xlabel="n_light_axions", ylabel="Count")
    hist!(ax4, n_light_list, bins=50, color=:steelblue, alpha=0.7)
    save(joinpath(PLOT_DIR, filename), f)
    log_message("Saved geometry abundance plot: $filename")
    return f
end
