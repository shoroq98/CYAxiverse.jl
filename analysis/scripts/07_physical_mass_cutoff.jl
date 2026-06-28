include(joinpath(@__DIR__, "00_setup.jl"))

log_message("="^60)
log_message("STAGE 7: Physical Mass Cutoff Comparison")
log_message("Comparing Version A (upper only) vs Version B (physical window)")
log_message("="^60)

const DM_DENSITY = 0.12
const H0_LOG10 = -33.0
const MAX_LOG10_MASS = 2.0
const CANDIDATE_MIN = 0.5
const CANDIDATE_MAX = 1.0

theta_vals = collect(range(0.01, π, length=10))
cache_label = ""

axion_data = load_dataset("axion_dataset$(cache_label).jls")
geometry_data = load_dataset("geometry_dataset$(cache_label).jls")
axions = axion_data.axions
geometries = geometry_data.geometries
log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

function log10_density_ratio_one_axion(log10m::Float64, log10f::Float64, theta::Float64)
    theta == 0.0 && return -Inf
    return log10(0.4 / DM_DENSITY) +
           2.0 * log10(abs(theta / (π / 2.0))) +
           0.5 * (log10m + 17.0) +
           2.0 * (log10f - 16.0)
end

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

# ─────────────────────────────────────────────────────────────────────────────
# Run both versions
# ─────────────────────────────────────────────────────────────────────────────

function run_version(axions, geometries, theta_vals, label;
                     mass_mask_fn)
    geom_ids_ratio = Int[]
    log10_total_ratios = Float64[]
    theta_used = Float64[]
    candidate_rows = []
    n_axions_included = zeros(Int, length(geometries))

    for (gid, g) in enumerate(geometries)
        matching = [a for a in axions
                    if a.h11 == g.h11 &&
                       a.polytope == g.polytope &&
                       a.frst == g.frst]

        included = [a for a in matching
                    if isfinite(a.log10m) && isfinite(a.log10f) && mass_mask_fn(a.log10m)]
        n_axions_included[gid] = length(included)

        for theta in theta_vals
            log10_ratios = Float64[]
            for a in included
                push!(log10_ratios, log10_density_ratio_one_axion(a.log10m, a.log10f, theta))
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

        if gid % 50 == 0
            log_message("  $label: processed $gid/$(length(geometries)) geometries")
        end
    end

    return (geom_ids_ratio=geom_ids_ratio,
            log10_total_ratios=log10_total_ratios,
            theta_used=theta_used,
            candidates=candidate_rows,
            n_axions_included=n_axions_included)
end

log_message("Running Version A (upper cutoff only: m < 10² eV)...")
result_A = run_version(axions, geometries, theta_vals, "Version A";
    mass_mask_fn = log10m -> log10m <= MAX_LOG10_MASS)

log_message("Running Version B (physical window: 10⁻³³ < m < 10² eV)...")
result_B = run_version(axions, geometries, theta_vals, "Version B";
    mass_mask_fn = log10m -> (H0_LOG10 <= log10m <= MAX_LOG10_MASS))

# ─────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────────

n_geoms_A = length(unique(result_A.geom_ids_ratio))
n_geoms_B = length(unique(result_B.geom_ids_ratio))
n_cand_A = length(result_A.candidates)
n_cand_B = length(result_B.candidates)
zero_ax_A = count(iszero, result_A.n_axions_included)
zero_ax_B = count(iszero, result_B.n_axions_included)
n_total_ax_A = sum(result_A.n_axions_included)
n_total_ax_B = sum(result_B.n_axions_included)

log_message("="^60)
log_message("Comparison Summary:")
log_message("  $(lpad("", 35)) Version A (upper)   Version B (physical)")
log_message("  $(lpad("", 35)) $(lpad("─", 20)) $(lpad("─", 20))")
log_message("  Total axions included:  $(lpad(n_total_ax_A, 12))  $(lpad(n_total_ax_B, 12))")
log_message("  Geometries with DM ≥0:  $(lpad(n_geoms_A, 12))  $(lpad(n_geoms_B, 12))")
log_message("  DM candidates (R≈1):    $(lpad(n_cand_A, 12))  $(lpad(n_cand_B, 12))")
log_message("  Geoms w/ 0 axions:      $(lpad(zero_ax_A, 12))  $(lpad(zero_ax_B, 12))")
log_message("="^60)

result_comparison = (
    result_A = result_A, result_B = result_B,
    n_total_ax_A = n_total_ax_A, n_total_ax_B = n_total_ax_B,
    n_geoms_A = n_geoms_A, n_geoms_B = n_geoms_B,
    n_cand_A = n_cand_A, n_cand_B = n_cand_B,
    zero_ax_A = zero_ax_A, zero_ax_B = zero_ax_B
)

save_dataset(result_comparison, "physical_mass_cutoff_comparison$(cache_label).jls")

# ─────────────────────────────────────────────────────────────────────────────
# Load λ_self data if available
# ─────────────────────────────────────────────────────────────────────────────

lambda_available = isfile(joinpath(OUTPUT_DIR, "quartic_lambda_results.jls"))
lambda_geom_ids = Int[]
lambda_max_vals = Float64[]
if lambda_available
    lambda_data = deserialize(joinpath(OUTPUT_DIR, "quartic_lambda_results.jls"))
    lambda_geom_ids = lambda_data.geom_ids
    lambda_max_vals = lambda_data.lambda_max
    log_message("Loaded λ_self data for $(length(lambda_geom_ids)) geometries")
end

# ─────────────────────────────────────────────────────────────────────────────
# Best-theta per geometry for both versions
# ─────────────────────────────────────────────────────────────────────────────

function best_theta_per_geom(geom_ids_ratio, log10_total_ratios, theta_used)
    unique_geoms = sort(unique(geom_ids_ratio))
    best_ids = Int[]
    best_thetas = Float64[]
    best_logs = Float64[]
    best_ratios = Float64[]
    for gid in unique_geoms
        inds = findall(geom_ids_ratio .== gid)
        isempty(inds) && continue
        vals = log10_total_ratios[inds]
        thetas = theta_used[inds]
        best_idx = argmin(abs.(vals .- 0.0))
        push!(best_ids, gid)
        push!(best_thetas, thetas[best_idx])
        push!(best_logs, vals[best_idx])
        push!(best_ratios, 10.0^vals[best_idx])
    end
    return (geom_ids=best_ids, theta_best=best_thetas,
            log10_best=best_logs, best_ratio=best_ratios)
end

best_A = best_theta_per_geom(result_A.geom_ids_ratio, result_A.log10_total_ratios, result_A.theta_used)
best_B = best_theta_per_geom(result_B.geom_ids_ratio, result_B.log10_total_ratios, result_B.theta_used)

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────

CairoMakie.activate!(type="png")

colors = (A=(:dodgerblue, 0.5), B=(:tomato, 0.5))

# ── 1. Histogram: abundance ratio overlay ──────────────────────────────────

fig1 = Figure(size=(1200, 700), fontsize=16)
ax1 = Axis(fig1[1, 1],
    title="Abundance Distribution: Version A (upper) vs Version B (physical window)",
    xlabel="log10(Ω_geom / Ω_DM)", ylabel="Number of geometries",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)

hist!(ax1, best_A.log10_best, bins=30,
    color=colors.A, strokecolor=:dodgerblue, strokewidth=1.5,
    label="Version A (upper only, n=$(length(best_A.geom_ids)))")
hist!(ax1, best_B.log10_best, bins=30,
    color=colors.B, strokecolor=:tomato, strokewidth=1.5,
    label="Version B (physical window, n=$(length(best_B.geom_ids)))")
vlines!(ax1, [0.0], color=:black, linestyle=:dash, linewidth=2,
    label="Ω_geom = Ω_DM")
axislegend(ax1, position=:lt, labelsize=13)
save(joinpath(PLOT_DIR, "mass_cutoff_abundance_histogram.png"), fig1)
log_message("Saved: mass_cutoff_abundance_histogram.png")

# ── 2. Abundance vs h11 ────────────────────────────────────────────────────

h11_list = [g.h11 for g in geometries]

geom_h11_A = [geometries[gid].h11 for gid in best_A.geom_ids]
geom_h11_B = [geometries[gid].h11 for gid in best_B.geom_ids]

fig2 = Figure(size=(1200, 700), fontsize=16)
ax2 = Axis(fig2[1, 1],
    title="Best Density Ratio vs h11",
    xlabel="h11", ylabel="log10(Ω_geom / Ω_DM)",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)
scatter!(ax2, geom_h11_A, best_A.log10_best,
    markersize=8, color=(:dodgerblue, 0.5),
    label="Version A (upper only)")
scatter!(ax2, geom_h11_B .+ 0.1, best_B.log10_best,
    markersize=8, color=(:tomato, 0.5),
    label="Version B (physical window)")
hlines!(ax2, [0.0], color=:black, linestyle=:dash, linewidth=2)
axislegend(ax2, position=:rb, labelsize=13)
save(joinpath(PLOT_DIR, "mass_cutoff_abundance_vs_h11.png"), fig2)
log_message("Saved: mass_cutoff_abundance_vs_h11.png")

# ── 3. Abundance vs λ_self (if available) ──────────────────────────────────

if lambda_available
    lambda_lookup = Dict(lambda_geom_ids .=> lambda_max_vals)

    lambda_vals_A = Float64[get(lambda_lookup, gid, NaN) for gid in best_A.geom_ids]
    lambda_vals_B = Float64[get(lambda_lookup, gid, NaN) for gid in best_B.geom_ids]
    log10_best_lambda_A = best_A.log10_best
    log10_best_lambda_B = best_B.log10_best

    finite_mask_A = isfinite.(lambda_vals_A) .& isfinite.(log10_best_lambda_A)
    finite_mask_B = isfinite.(lambda_vals_B) .& isfinite.(log10_best_lambda_B)

    fig3 = Figure(size=(1200, 700), fontsize=16)
    ax3 = Axis(fig3[1, 1],
        title="Best Density Ratio vs log10(max λ_self)",
        xlabel="log10(Ω_geom / Ω_DM)", ylabel="log10(max λ_self)",
        titlesize=18, xlabelsize=16, ylabelsize=16,
        xticklabelsize=14, yticklabelsize=14)
    scatter!(ax3, log10_best_lambda_A[finite_mask_A], lambda_vals_A[finite_mask_A],
        markersize=8, color=(:dodgerblue, 0.5),
        label="Version A (upper only)")
    scatter!(ax3, log10_best_lambda_B[finite_mask_B], lambda_vals_B[finite_mask_B],
        markersize=8, color=(:tomato, 0.5),
        label="Version B (physical window)")
    vlines!(ax3, [0.0], color=:black, linestyle=:dash, linewidth=2)
    axislegend(ax3, position=:rb, labelsize=13)
    save(joinpath(PLOT_DIR, "mass_cutoff_abundance_vs_lambda.png"), fig3)
    log_message("Saved: mass_cutoff_abundance_vs_lambda.png")
end

# ── 4. Surviving axions per geometry: boxplot side-by-side ─────────────────

fig4 = Figure(size=(1200, 700), fontsize=16)
ax4 = Axis(fig4[1, 1],
    title="Number of Surviving Axions per Geometry",
    xlabel="Version", ylabel="Number of axions",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)

boxplot!(ax4, fill(1.0, length(result_A.n_axions_included)), Float64.(result_A.n_axions_included),
    width=0.5, color=colors.A, strokecolor=:dodgerblue, strokewidth=1.5,
    whiskerwidth=0.5, label="Version A (upper only)")
boxplot!(ax4, fill(2.0, length(result_B.n_axions_included)), Float64.(result_B.n_axions_included),
    width=0.5, color=colors.B, strokecolor=:tomato, strokewidth=1.5,
    whiskerwidth=0.5, label="Version B (physical window)")
ax4.xticks = ([1, 2], ["Version A\n(upper only)", "Version B\n(physical window)"])
axislegend(ax4, position=:rt, labelsize=13)
save(joinpath(PLOT_DIR, "mass_cutoff_surviving_axions.png"), fig4)
log_message("Saved: mass_cutoff_surviving_axions.png")

# ── 5. Bar chart: zero-axion geometries ────────────────────────────────────

fig5 = Figure(size=(800, 600), fontsize=16)
ax5 = Axis(fig5[1, 1],
    title="Geometries with Zero Surviving Axions",
    xlabel="Version", ylabel="Number of geometries",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)
barplot!(ax5, 1:2, [zero_ax_A, zero_ax_B],
    color=[:dodgerblue, :tomato])
ax5.xticks = ([1, 2], ["Version A\n(upper only)", "Version B\n(physical window)"])
save(joinpath(PLOT_DIR, "mass_cutoff_zero_axion_geoms.png"), fig5)
log_message("Saved: mass_cutoff_zero_axion_geoms.png")

# ── 6. Text report ─────────────────────────────────────────────────────────

report_path = joinpath(OUTPUT_DIR, "physical_mass_cutoff_report$(cache_label).txt")
open(report_path, "w") do io
    println(io, "="^70)
    println(io, "  Physical Mass Cutoff Comparison Report")
    println(io, "  H₀ ≈ 10⁻³³ eV  |  m_max = 10² eV")
    println(io, "="^70)
    println(io)
    println(io, "  Version A: upper cutoff only  (m < 10² eV)")
    println(io, "  Version B: physical window    (10⁻³³ < m < 10² eV)")
    println(io)
    println(io, "  $(lpad("Metric", 35))  $(lpad("Version A", 15))  $(lpad("Version B", 15))")
    println(io, "  $(lpad("─"^35, 35))  $(lpad("─"^15, 15))  $(lpad("─"^15, 15))")
    println(io, "  $(lpad("Total axions included", 35))  $(lpad(n_total_ax_A, 15))  $(lpad(n_total_ax_B, 15))")
    println(io, "  $(lpad("Geometries with DM ≥ 0", 35))  $(lpad(n_geoms_A, 15))  $(lpad(n_geoms_B, 15))")
    println(io, "  $(lpad("DM candidates (R ≈ 1)", 35))  $(lpad(n_cand_A, 15))  $(lpad(n_cand_B, 15))")
    println(io, "  $(lpad("Geometries w/ 0 axions", 35))  $(lpad(zero_ax_A, 15))  $(lpad(zero_ax_B, 15))")
    println(io)
    println(io, "  Legend:")
    println(io, "    Version A: only applies upper mass cut (m < 10² eV)")
    println(io, "    Version B: applies Hubble-scale lower cut + upper cut")
    println(io, "              (H₀ ≈ 10⁻³³ eV < m < 10² eV)")
    println(io)
    println(io, "  Dark Matter Candidates (Version A):")
    for c in result_A.candidates
        println(io, "    geom_id=$(c.geom_id) h11=$(c.h11) theta=$(round(c.theta, digits=5)) R=$(round(c.ratio, digits=6))")
    end
    println(io)
    println(io, "  Dark Matter Candidates (Version B):")
    for c in result_B.candidates
        println(io, "    geom_id=$(c.geom_id) h11=$(c.h11) theta=$(round(c.theta, digits=5)) R=$(round(c.ratio, digits=6))")
    end
end
log_message("Saved report: $report_path")

log_message("="^60)
log_message("Physical mass cutoff comparison complete!")
log_message("Saved 5+ comparison plots to $(PLOT_DIR)")
log_message("="^60)
