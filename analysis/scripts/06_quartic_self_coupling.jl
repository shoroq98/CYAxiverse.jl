using Pkg
Pkg.activate("/opt/CYAxiverse.jl")

using LinearAlgebra
using Statistics
using Serialization
using HDF5
using Printf
using Dates
using CairoMakie

ENV["MOSEKLM_LICENSE_FILE"] = "/home/cytools/mosek.lic"
ENV["newARGS"] = "docker"

PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
OUTPUT_DIR = joinpath(PROJECT_ROOT, "outputs")
PLOT_DIR = joinpath(PROJECT_ROOT, "plots")
LOG_DIR = joinpath(PROJECT_ROOT, "logs")
LOG_FILE = joinpath(LOG_DIR, "run_log.txt")
DB_DIR = "/scratch/database"

for d in [OUTPUT_DIR, PLOT_DIR, LOG_DIR]
    isdir(d) || mkpath(d)
end

function log_message(msg)
    open(LOG_FILE, "a") do io
        println(io, "[$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))] ", msg)
    end
    @info msg
end

const DM_DENSITY = 0.12
const THETA_VALS = collect(range(0.01, π, length=10))

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

function logsumexp_signed(log10_terms::Vector{Float64}, signs::Vector{Float64})
    finite = isfinite.(log10_terms)
    any(finite) || return (0.0, -Inf)
    vals = log10_terms[finite]
    sgns = signs[finite]
    max_val = maximum(vals)
    shifted = 10.0 .^ (vals .- max_val)
    total = sum(sgns .* shifted)
    if total == 0.0
        return (0, -Inf)
    elseif total > 0
        return (1, max_val + log10(total))
    else
        return (-1, max_val + log10(-total))
    end
end

function log10_axion_density_ratio(log10m::Float64, log10f::Float64, theta::Float64)
    theta == 0.0 && return -Inf
    return log10(0.4 / DM_DENSITY) +
           2.0 * log10(abs(theta / (π / 2.0))) +
           0.5 * (log10m + 17.0) +
           2.0 * (log10f - 16.0)
end

function compute_lambda_self_float64(L, Q, K)
    h11 = size(K, 1)
    signL = Float64.(L[:, 1])
    log10L = Float64.(L[:, 2])
    absL = signL .* 10.0 .^ log10L
    Qf = Matrix{Float64}(Q)
    @assert size(Qf, 2) == h11 "Q cols=$(size(Qf,2)) ≠ h11=$h11 (Q shape=$(size(Qf)))"
    @assert size(K, 1) == h11 "K=$(size(K)) ≠ h11=$h11"
    hess = Symmetric(Qf' * Diagonal(absL) * Qf)
    K_sym = Symmetric(Matrix{Float64}(K))
    eVals, eVecs = eigen(hess, K_sym)
    QMs = Qf * eVecs
    signQMs = sign.(QMs)
    log10QMs = log10.(abs.(QMs) .+ 1e-300)
    λ_log10 = zeros(Float64, h11)
    log10_2π = log10(2π)
    for k in 1:h11
        log_terms = log10L .+ 4.0 .* log10QMs[:, k]
        sgn, lse = logsumexp_signed(log_terms, signL)
        λ_log10[k] = sgn == 0 ? -Inf : lse + 4.0 * log10_2π
    end
    return λ_log10
end

function read_geometry_hdf5(h11, polytope, frst)
    h11_str = @sprintf("h11_%03d", h11)
    np_str = @sprintf("np_%07d", polytope)
    cy_str = @sprintf("cy_%07d", frst)
    fname = joinpath(DB_DIR, h11_str, np_str, cy_str, "cyax.h5")
    isfile(fname) || return nothing
    L = h5read(fname, "cytools/potential/L")
    Q = h5read(fname, "cytools/potential/Q")
    Kinv = h5read(fname, "cytools/geometric/Kinv")
    K = inv(Symmetric(Kinv))
    return (L=L, Q=Q, K=K)
end

log_message("="^60)
log_message("STAGE 6: Quartic Self-Coupling Analysis")
log_message("="^60)

log_message("Loading cached datasets...")
axion_data = deserialize(joinpath(OUTPUT_DIR, "axion_dataset.jls"))
axions = axion_data.axions
geom_data = deserialize(joinpath(OUTPUT_DIR, "geometry_dataset.jls"))
geometries = geom_data.geometries
log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

log_message("Computing λ_self and DM density ratios...")

geom_ids_all = Int[]
lambda_max_all = Float64[]
dm_ratio_best_all = Float64[]
candidate_flags = Bool[]

for (gid, g) in enumerate(geometries)
    result = read_geometry_hdf5(g.h11, g.polytope, g.frst)
    if result === nothing
        log_message("  WARNING: No data for geometry $gid (h11=$(g.h11))")
        continue
    end

    Q_ = result.Q
    K_ = result.K
    L_ = result.L

    if size(Q_, 2) != size(K_, 1)
        log_message("  SKIP geom $gid: size mismatch h11=$(g.h11) np=$(g.polytope) Q=$(size(Q_)) K=$(size(K_))")
        continue
    end

    λ_max = -Inf
    try
        λs = compute_lambda_self_float64(L_, Q_, K_)
        finite_λs = filter(isfinite, λs)
        λ_max = isempty(finite_λs) ? -Inf : maximum(finite_λs)
    catch e
        bt = catch_backtrace()
        msg = sprint(showerror, e, bt)
        log_message("  ERROR geom $gid h11=$(g.h11) np=$(g.polytope) frst=$(g.frst): $msg")
        continue
    end

    geom_axions = filter(a -> a.h11 == g.h11 && a.polytope == g.polytope && a.frst == g.frst, axions)
    log10_ratios = Float64[]
    for a in geom_axions
        for θ in THETA_VALS
            r = log10_axion_density_ratio(a.log10m, a.log10f, θ)
            isfinite(r) && push!(log10_ratios, r)
        end
    end
    best_log10_ratio = isempty(log10_ratios) ? -Inf : logsumexp_gaussian(log10_ratios)

    push!(geom_ids_all, gid)
    push!(lambda_max_all, λ_max)
    push!(dm_ratio_best_all, best_log10_ratio)

    ratio_linear = isfinite(best_log10_ratio) ? 10.0^best_log10_ratio : 0.0
    push!(candidate_flags, 0.5 <= ratio_linear <= 1.0)

    if gid % 50 == 0
        log_message("  Processed $gid/$(length(geometries))")
    end
end

log_message("Computed λ for $(length(geom_ids_all)) geometries")
n_candidates = count(candidate_flags)
log_message("DM candidates (0.5 <= R <= 1.0): $n_candidates")

finite_mask = isfinite.(lambda_max_all) .& isfinite.(dm_ratio_best_all)

CairoMakie.activate!(type="png")

fig = Figure(size=(1000, 800))
ax = Axis(fig[1, 1],
    title="Quartic Self-Coupling vs Best DM Density Ratio\nper Geometry (h11=4..20)",
    xlabel="best log10(Ω_geom / Ω_DM)", ylabel="log10(max λ_self)",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)
scatter!(ax, dm_ratio_best_all[finite_mask], lambda_max_all[finite_mask],
    markersize=6, color=(:steelblue, 0.5))
vlines!(ax, [0.0], color=:red, linestyle=:dash, linewidth=2,
    label="Ω_geom = Ω_DM")
axislegend(ax, position=:rb, labelsize=14)
save(joinpath(PLOT_DIR, "lambda_vs_dm_ratio.png"), fig)
log_message("Saved plot: lambda_vs_dm_ratio.png")

cand_λ = lambda_max_all[candidate_flags .& finite_mask]
rest_λ = lambda_max_all[.!candidate_flags .& finite_mask]

fig2 = Figure(size=(1000, 800))
ax2 = Axis(fig2[1, 1],
    title="λ_self Distribution: DM Candidates vs Rest\nh11=4..20",
    xlabel="log10(max λ_self)", ylabel="Number of geometries",
    titlesize=18, xlabelsize=16, ylabelsize=16,
    xticklabelsize=14, yticklabelsize=14)
if !isempty(cand_λ)
    hist!(ax2, cand_λ, bins=25, color=(:seagreen, 0.6),
        strokecolor=:black, strokewidth=1,
        label="DM candidates (R≈1) n=$(length(cand_λ))")
end
if !isempty(rest_λ)
    hist!(ax2, rest_λ, bins=25, color=(:tomato, 0.4),
        strokecolor=:black, strokewidth=1,
        label="Rest of landscape n=$(length(rest_λ))")
end
axislegend(ax2, position=:rt, labelsize=14)
save(joinpath(PLOT_DIR, "lambda_histogram_candidates.png"), fig2)
log_message("Saved plot: lambda_histogram_candidates.png")

serialize(joinpath(OUTPUT_DIR, "quartic_lambda_results.jls"),
    (geom_ids=geom_ids_all, lambda_max=lambda_max_all,
     dm_ratio_best=dm_ratio_best_all, candidate_flags=candidate_flags,
     finite_mask=finite_mask))

log_message("="^60)
log_message("Summary:")
log_message("  Geometries processed: $(length(geom_ids_all))")
log_message("  DM candidates: $n_candidates")
n_cand = length(cand_λ)
n_rest = length(rest_λ)
if n_cand > 0
    log_message("  Mean log10(λ) candidates: $(round(mean(cand_λ), digits=3)) ± $(round(std(cand_λ), digits=3))")
end
if n_rest > 0
    log_message("  Mean log10(λ) rest:      $(round(mean(rest_λ), digits=3)) ± $(round(std(rest_λ), digits=3))")
end
log_message("="^60)
log_message("Quartic self-coupling analysis complete!")
