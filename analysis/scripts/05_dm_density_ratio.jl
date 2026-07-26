include(joinpath(@__DIR__, "00_setup.jl"))

log_message("="^60)
log_message("STAGE 5: Dark Matter Density Ratio Analysis")
log_message("="^60)

const DM_DENSITY = 0.12
const MAX_LOG10_MASS = 2.0
const H0_LOG10_MASS = -33.0
const CANDIDATE_MIN = 0.5
const CANDIDATE_MAX = 1.0

theta_vals = collect(range(0.01, π, length=10))
cache_label = ""

axion_data = load_dataset("axion_dataset$(cache_label).jls")
geometry_data = load_dataset("geometry_dataset$(cache_label).jls")
axions = axion_data.axions
geometries = geometry_data.geometries
log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

function gauss_sum(z::Float64)
    log2 = log(2.0)
    if abs(z) > 600.0
        return 0.5 * z + abs(0.5 * z)
    else
        return log2 + 0.5 * z + log(cosh(0.5 * z))
    end
end

function logsumexp_gaussian(log_values::Vector{Float64})
    finite_vals = filter(isfinite, log_values)
    isempty(finite_vals) && return -Inf
    sort!(finite_vals, rev=true)
    total_log = finite_vals[1]
    for v in finite_vals[2:end]
        total_log = total_log + gauss_sum(v - total_log)
    end
    return total_log
end

function log10sumexp_gaussian(log10_values::Vector{Float64})
    ln_values = log(10.0) .* log10_values
    ln_total = logsumexp_gaussian(ln_values)
    return ln_total / log(10.0)
end

function log10_density_ratio_one_axion(log10m::Float64, log10f::Float64, theta::Float64)
    theta == 0.0 && return -Inf
    return log10(11.3) +
           2.0 * log10(abs(theta / (π / 2.0))) +
           0.5 * (log10m + 17.0) +
           2.0 * (log10f - 16.0)
end

CairoMakie.activate!(type="png")
colors = (dodgerblue=(:dodgerblue, 0.4), tomato=(:tomato, 0.4), seagreen=(:seagreen, 0.4))

function log10_analysis_loop(axions, geometries, theta_vals, version_label, mask_fn)
    geom_ids_ratio = Int[]
    log10_total_ratios = Float64[]
    theta_used = Float64[]
    geom_ids_mass = Int[]
    mass_values = Float64[]
    geom_ids_f = Int[]
    f_values = Float64[]
    candidate_rows = []

    for (gid, g) in enumerate(geometries)
        matching = [a for a in axions
                    if a.h11 == g.h11 &&
                       a.polytope == g.polytope &&
                       a.frst == g.frst]

        for (i, a) in enumerate(matching)
            if isfinite(a.log10m) && isfinite(a.log10f) && mask_fn(a.log10m)
                push!(geom_ids_mass, gid)
                push!(mass_values, a.log10m)
                push!(geom_ids_f, gid)
                push!(f_values, a.log10f)
            end
        end

        for theta in theta_vals
            log10_ratios = Float64[]
            for a in matching
                if isfinite(a.log10m) && isfinite(a.log10f) && mask_fn(a.log10m)
                    push!(log10_ratios, log10_density_ratio_one_axion(a.log10m, a.log10f, theta))
                end
            end
            log10_total_ratio = log10sumexp_gaussian(log10_ratios)
            if isfinite(log10_total_ratio)
                ratio = 10.0^log10_total_ratio
                push!(geom_ids_ratio, gid)
                push!(log10_total_ratios, log10_total_ratio)
                push!(theta_used, theta)
                if CANDIDATE_MIN <= ratio <= CANDIDATE_MAX
                    push!(candidate_rows, (
                        geom_id=gid, h11=g.h11, theta=theta, ratio=ratio
                    ))
                end
            end
        end
    end

    return (geom_ids_ratio=geom_ids_ratio, log10_total_ratios=log10_total_ratios,
            theta_used=theta_used,
            geom_ids_mass=geom_ids_mass, mass_values=mass_values,
            geom_ids_f=geom_ids_f, f_values=f_values,
            candidates=candidate_rows)
end

function make_plots_and_reports(result, version_label, suffix)
    geom_ids_ratio = result.geom_ids_ratio
    log10_total_ratios = result.log10_total_ratios
    theta_used = result.theta_used
    geom_ids_mass = result.geom_ids_mass
    mass_values = result.mass_values
    geom_ids_f = result.geom_ids_f
    f_values = result.f_values
    candidate_rows = result.candidates

    fig = Figure(size=(1600, 1000), fontsize=18)

    ax1 = Axis(fig[1, 1:2],
        title="($version_label a) Sum of Density Ratios per Geometry",
        xlabel="Geometry index", ylabel="log10(Ω_geom h² / 0.12)",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)

    boxplot!(ax1, geom_ids_ratio, log10_total_ratios,
        width=0.6, color=colors.dodgerblue, strokecolor=:black, strokewidth=1.5,
        whiskerwidth=0.5, show_outliers=false)
    scatter!(ax1, geom_ids_ratio, log10_total_ratios,
        markersize=4, color=(:black, 0.3))
    hlines!(ax1, [0.0], color=:red, linestyle=:dash, linewidth=2.5)

    ax2 = Axis(fig[2, 1],
        title="($version_label b) Mass Distribution",
        xlabel="Geometry index", ylabel="log10(m_a / eV)",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)
    boxplot!(ax2, geom_ids_mass, mass_values,
        width=0.6, color=colors.tomato, strokecolor=:black, strokewidth=1.5,
        whiskerwidth=0.5, show_outliers=false)
    scatter!(ax2, geom_ids_mass, mass_values,
        markersize=3, color=(:black, 0.3))
    hlines!(ax2, [MAX_LOG10_MASS], color=:red, linestyle=:dash, linewidth=2,
        label="max mass cut (10² eV)")
    axislegend(ax2, position=:rb, titlesize=16, labelsize=14)

    ax3 = Axis(fig[2, 2],
        title="($version_label c) Decay Constant Distribution",
        xlabel="Geometry index", ylabel="log10(f / GeV)",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)
    boxplot!(ax3, geom_ids_f, f_values,
        width=0.6, color=colors.seagreen, strokecolor=:black, strokewidth=1.5,
        whiskerwidth=0.5, show_outliers=false)
    scatter!(ax3, geom_ids_f, f_values,
        markersize=3, color=(:black, 0.3))

    max_gid = length(geom_ids_ratio) > 0 ? maximum(geom_ids_ratio) : 0
    tick_step = max(1, ceil(Int, max(max_gid, 1) / 25))
    xticks = 1:tick_step:max(max_gid, 1)
    for ax in [ax1, ax2, ax3]
        ax.xticks = xticks
    end

    save(joinpath(PLOT_DIR, "dm_density_ratio_boxplots_$(suffix).png"), fig)
    log_message("Saved plot: dm_density_ratio_boxplots_$(suffix).png")

    candidate_path = joinpath(OUTPUT_DIR, "dm_density_ratio_candidates_$(suffix).txt")
    open(candidate_path, "w") do io
        println(io, "="^60)
        println(io, "  Dark Matter Candidates ($version_label)")
        println(io, "="^60)
        if isempty(candidate_rows)
            println(io, "  No candidates found.")
        else
            println(io, "  geom_id  h11  theta       ratio")
            println(io, "  " ^ 35)
            for c in candidate_rows
                println(io, "  $(lpad(c.geom_id, 3))      $(c.h11)   $(rpad(round(c.theta, digits=5), 10)) $(round(c.ratio, digits=6))")
            end
        end
        println(io, "="^60)
    end
    log_message("Saved candidate report: $(suffix)")

    best_geom_ids = Int[]
    best_theta_over_pi = Float64[]
    best_theta_values = Float64[]
    best_log10_ratios = Float64[]
    best_ratios = Float64[]

    unique_geoms = sort(unique(geom_ids_ratio))
    for gid in unique_geoms
        inds = findall(geom_ids_ratio .== gid)
        isempty(inds) && continue
        vals = log10_total_ratios[inds]
        thetas = theta_used[inds]
        local_best_index = argmin(abs.(vals .- 0.0))
        best_log = vals[local_best_index]
        best_theta = thetas[local_best_index]
        push!(best_geom_ids, gid)
        push!(best_theta_values, best_theta)
        push!(best_theta_over_pi, best_theta / π)
        push!(best_log10_ratios, best_log)
        push!(best_ratios, 10.0^best_log)
    end

    fig_best_theta = Figure(size=(1400, 700))
    ax_best_theta = Axis(fig_best_theta[1, 1],
        title="$version_label: Best theta per geometry",
        xlabel="Geometry index", ylabel="theta_best / pi",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)
    scatter!(ax_best_theta, best_geom_ids, best_theta_over_pi,
        markersize=8, color=best_log10_ratios, colormap=:viridis)
    Colorbar(fig_best_theta[1, 2],
        label="log10(Omega_geom / Omega_DM)",
        colormap=:viridis,
        limits=(minimum(best_log10_ratios), maximum(best_log10_ratios)),
        labelsize=18, ticklabelsize=14)
    max_gid_best = length(best_geom_ids) > 0 ? maximum(best_geom_ids) : 0
    tick_step_best = max(1, ceil(Int, max(max_gid_best, 1) / 20))
    ax_best_theta.xticks = 1:tick_step_best:max(max_gid_best, 1)
    save(joinpath(PLOT_DIR, "best_theta_per_geometry_$(suffix).png"), fig_best_theta)
    log_message("Saved plot: best_theta_per_geometry_$(suffix).png")

    fig_best_ratio = Figure(size=(1400, 700))
    ax_best_ratio = Axis(fig_best_ratio[1, 1],
        title="$version_label: Closest density ratio per geometry",
        xlabel="Geometry index", ylabel="log10(Omega_geom / Omega_DM)",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)
    scatter!(ax_best_ratio, best_geom_ids, best_log10_ratios,
        markersize=8, color=:black)
    hlines!(ax_best_ratio, [0.0], color=:red, linestyle=:dash, linewidth=3)
    ax_best_ratio.xticks = 1:tick_step_best:max(max_gid_best, 1)
    save(joinpath(PLOT_DIR, "best_ratio_per_geometry_$(suffix).png"), fig_best_ratio)
    log_message("Saved plot: best_ratio_per_geometry_$(suffix).png")

    save_dataset(
        (geom_ids=best_geom_ids, theta_best=best_theta_values,
         theta_best_over_pi=best_theta_over_pi,
         log10_best_ratio=best_log10_ratios, best_ratio=best_ratios),
        "best_theta_per_geometry_$(suffix).jls"
    )

    best_theta_report = joinpath(OUTPUT_DIR, "best_theta_per_geometry_$(suffix).txt")
    open(best_theta_report, "w") do io
        println(io, "$version_label: Best theta per geometry")
        println(io, "="^70)
        println(io, "geom_id, theta_best, theta_best/pi, log10_best_ratio, best_ratio")
        for i in eachindex(best_geom_ids)
            println(io, "$(best_geom_ids[i]), $(round(best_theta_values[i], digits=6)), $(round(best_theta_over_pi[i], digits=6)), $(round(best_log10_ratios[i], digits=6)), $(round(best_ratios[i], digits=6))")
        end
    end
    log_message("Saved report: best_theta_per_geometry_$(suffix).txt")

    fig_theta_hist = Figure(size=(900, 600))
    ax_theta_hist = Axis(fig_theta_hist[1, 1],
        title="$version_label: Distribution of best theta values",
        xlabel="theta_best / pi", ylabel="Number of geometries",
        titlesize=22, xlabelsize=18, ylabelsize=18,
        xticklabelsize=14, yticklabelsize=14)
    hist!(ax_theta_hist, best_theta_over_pi, bins=20,
        color=(:steelblue, 0.65), strokecolor=:black, strokewidth=1)
    vlines!(ax_theta_hist, [0.0, 1.0], color=:red, linestyle=:dash, linewidth=2)
    save(joinpath(PLOT_DIR, "best_theta_histogram_$(suffix).png"), fig_theta_hist)
    log_message("Saved plot: best_theta_histogram_$(suffix).png")

    candidate_theta_over_pi = Float64[]
    for i in eachindex(best_geom_ids)
        if 0.5 <= best_ratios[i] <= 1.0
            push!(candidate_theta_over_pi, best_theta_over_pi[i])
        end
    end

    if !isempty(candidate_theta_over_pi)
        fig_candidate_theta_hist = Figure(size=(900, 600))
        ax_candidate_theta_hist = Axis(fig_candidate_theta_hist[1, 1],
            title="$version_label: Best theta for DM candidates",
            xlabel="theta_best / pi", ylabel="Number of candidate geometries",
            titlesize=22, xlabelsize=18, ylabelsize=18,
            xticklabelsize=14, yticklabelsize=14)
        hist!(ax_candidate_theta_hist, candidate_theta_over_pi, bins=20,
            color=(:seagreen, 0.65), strokecolor=:black, strokewidth=1)
        vlines!(ax_candidate_theta_hist, [0.0, 1.0], color=:red, linestyle=:dash, linewidth=2)
        save(joinpath(PLOT_DIR, "best_theta_histogram_candidates_$(suffix).png"), fig_candidate_theta_hist)
        log_message("Saved plot: best_theta_histogram_candidates_$(suffix).png")
    end
end

# ═════════════════════════════════════════════════════════════════════
# Version A: upper cutoff only (m < 10² eV)
# ═════════════════════════════════════════════════════════════════════
log_message("Running Version A (upper cutoff only: m < 10² eV)...")
result_A = log10_analysis_loop(axions, geometries, theta_vals, "A",
    log10m -> isfinite(log10m) && log10m <= MAX_LOG10_MASS)

n_ax_A = sum(count(isfinite(a.log10m) && a.log10m <= MAX_LOG10_MASS for a in axions
    if a.h11 == g.h11 && a.polytope == g.polytope && a.frst == g.frst)
    for (gid, g) in enumerate(geometries))
log_message("Version A: $(n_ax_A) axions included, $(length(result_A.candidates)) candidates")

save_dataset(result_A, "dm_density_ratio_results_A.jls")
make_plots_and_reports(result_A, "Version A", "version_A")

# ═════════════════════════════════════════════════════════════════════
# Version B: physical mass window (H₀ < m < 10² eV)
# ═════════════════════════════════════════════════════════════════════
log_message("Running Version B (physical window: 10⁻³³ < m < 10² eV)...")
result_B = log10_analysis_loop(axions, geometries, theta_vals, "B",
    log10m -> isfinite(log10m) && H0_LOG10_MASS <= log10m <= MAX_LOG10_MASS)

n_frozen = count(a -> isfinite(a.log10m) && a.log10m < H0_LOG10_MASS, axions)
n_dm = count(a -> isfinite(a.log10m) && H0_LOG10_MASS <= a.log10m <= MAX_LOG10_MASS, axions)
n_heavy = count(a -> isfinite(a.log10m) && a.log10m > MAX_LOG10_MASS, axions)
log_message("Version B: $(n_dm) axions included, $(length(result_B.candidates)) candidates")
log_message("Classification: frozen=$n_frozen  DM=$n_dm  heavy=$n_heavy  total=$(n_frozen+n_dm+n_heavy)")

save_dataset(result_B, "dm_density_ratio_results_B.jls")
make_plots_and_reports(result_B, "Version B", "version_B")

log_message("="^60)
log_message("DM density ratio analysis complete!")
log_message("Files saved with _version_A and _version_B suffixes")
log_message("="^60)
