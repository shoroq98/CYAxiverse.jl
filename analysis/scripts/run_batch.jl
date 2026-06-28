using Pkg
Pkg.activate("/opt/CYAxiverse.jl")

using CYAxiverse
using LinearAlgebra
using Statistics
using Serialization
using HDF5
using Printf
using Dates

include(joinpath(@__DIR__, "..", "src", "abundance_tools.jl"))

# ─────────────────────────────────────────────────────────────────────
# CLI args
# ─────────────────────────────────────────────────────────────────────
h11_min = 4
h11_max = 20
output_base = "/mnt/results"

for (i, arg) in enumerate(ARGS)
    if arg == "--h11-min" && i < length(ARGS)
        h11_min = parse(Int, ARGS[i+1])
    elseif arg == "--h11-max" && i < length(ARGS)
        h11_max = parse(Int, ARGS[i+1])
    elseif arg == "--output" && i < length(ARGS)
        output_base = ARGS[i+1]
    end
end

batch_label = @sprintf("batch_%03d_%03d", h11_min, h11_max)
output_dir = joinpath(output_base, "batches", batch_label)
mkpath(output_dir)

log_file = joinpath(output_dir, "batch.log")

function log(msg)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    line = "[$timestamp] [Batch h11=$h11_min-$h11_max] $msg"
    open(log_file, "a") do io; println(io, line); end
    @info line
end

# ─────────────────────────────────────────────────────────────────────
# STAGE A: Find geometries
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("STAGE A: Finding geometries in h11=$h11_min..$h11_max")
log("="^60)

full_list = CYAxiverse.filestructure.paths_cy()[2]
mask = h11_min .<= full_list[1, :] .<= h11_max
h11list = full_list[:, mask]
n_total = size(h11list, 2)
log("Found $n_total geometries")

# ─────────────────────────────────────────────────────────────────────
# STAGE B: Generate axion + geometry datasets
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("STAGE B: Dataset generation")
log("="^60)

axion_cache = joinpath(output_dir, "axion_dataset.jls")
geom_cache  = joinpath(output_dir, "geometry_dataset.jls")
failed_cache = joinpath(output_dir, "failed_geometries.jls")

if isfile(axion_cache) && isfile(geom_cache)
    log("Using cached datasets")
    axion_data = deserialize(axion_cache)
    geom_data  = deserialize(geom_cache)
else
    log("Building axion dataset...")
    axion_data = _build_axion_dataset(h11list)
    serialize(axion_cache, axion_data)
    log("Building geometry dataset...")
    geom_data  = _build_geometry_dataset(h11list)
    serialize(geom_cache, geom_data)
end

axions = axion_data.axions
geoms  = geom_data.geometries
log("Loaded $(length(axions)) axions from $(length(geoms)) geometries")

# ─────────────────────────────────────────────────────────────────────
# STAGE C: DM density ratio
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("STAGE C: Dark Matter density ratio analysis")
log("="^60)

const DM_DENSITY = 0.12
const MAX_LOG10_MASS = 2.0
const H0_LOG10_MASS = -33.0
theta_vals = collect(range(0.01, π, length=10))

function _gauss_sum(z::Float64)
    log2 = log(2.0)
    if abs(z) > 600.0
        return 0.5 * z + abs(0.5 * z)
    else
        return log2 + 0.5 * z + log(cosh(0.5 * z))
    end
end

function _logsumexp_gaussian(log_values::Vector{Float64})
    finite_vals = filter(isfinite, log_values)
    isempty(finite_vals) && return -Inf
    sort!(finite_vals, rev=true)
    total_log = finite_vals[1]
    for v in finite_vals[2:end]
        total_log = total_log + _gauss_sum(v - total_log)
    end
    return total_log
end

function _log10_density_ratio(log10m::Float64, log10f::Float64, theta::Float64)
    theta == 0.0 && return -Inf
    log10(0.4 / DM_DENSITY) + 2.0 * log10(abs(theta / (π / 2.0))) +
    0.5 * (log10m + 17.0) + 2.0 * (log10f - 16.0)
end

dm_cache = joinpath(output_dir, "dm_density_ratio_results.jls")

if isfile(dm_cache)
    log("Using cached DM results")
    dm_result = deserialize(dm_cache)
else
    geom_ids_ratio = Int[]
    log10_total_ratios = Float64[]
    theta_used = Float64[]
    candidate_rows = []

    for (gid, g) in enumerate(geoms)
        matching = [a for a in axions
                    if a.h11 == g.h11 && a.polytope == g.polytope && a.frst == g.frst]

        dm_mask = dm_relevant_axion_mask([a.log10m for a in matching],
                                         log10m_min=H0_LOG10_MASS,
                                         log10m_max=MAX_LOG10_MASS)

        for theta in theta_vals
            log10_ratios = Float64[]
            for (i, a) in enumerate(matching)
                if dm_mask[i] && isfinite(a.log10m) && isfinite(a.log10f)
                    push!(log10_ratios, _log10_density_ratio(a.log10m, a.log10f, theta))
                end
            end
            total = _logsumexp_gaussian(log10_ratios)
            if isfinite(total)
                push!(geom_ids_ratio, gid)
                push!(log10_total_ratios, total)
                push!(theta_used, theta)
                ratio_lin = 10.0^total
                if 0.5 <= ratio_lin <= 1.0
                    push!(candidate_rows, (
                        geom_id=gid, h11=g.h11, theta=theta, ratio=ratio_lin
                    ))
                end
            end
        end

        if gid % 50 == 0
            log("  DM: processed $gid/$(length(geoms)) geometries")
        end
    end

    dm_result = (geom_ids_ratio=geom_ids_ratio, log10_total_ratios=log10_total_ratios,
                 theta_used=theta_used, candidates=candidate_rows)
    serialize(dm_cache, dm_result)
    log("DM ratio done: $(length(candidate_rows)) candidates")
end

# ─────────────────────────────────────────────────────────────────────
# STAGE D: Quartic self-coupling via HDF5
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("STAGE D: Quartic self-coupling (λ_self) via HDF5")
log("="^60)

const DB_DIR = "/scratch/database"
const LOG10_2π = log10(2π)

function _logsumexp_signed(log10_terms::Vector{Float64}, signs::Vector{Float64})
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

function _compute_lambda_self(L, Q, K)
    h11 = size(K, 1)
    signL = Float64.(L[:, 1])
    log10L = Float64.(L[:, 2])
    absL = signL .* 10.0 .^ log10L
    Qf = Matrix{Float64}(Q)
    hess = Symmetric(Qf' * Diagonal(absL) * Qf)
    K_sym = Symmetric(Matrix{Float64}(K))
    eVals, eVecs = eigen(hess, K_sym)
    QMs = Qf * eVecs
    λ_log10 = zeros(Float64, h11)
    for k in 1:h11
        log_terms = log10L .+ 4.0 .* log10.(abs.(QMs[:, k]) .+ 1e-300)
        sgn, lse = _logsumexp_signed(log_terms, signL)
        λ_log10[k] = sgn == 0 ? -Inf : lse + 4.0 * LOG10_2π
    end
    return λ_log10
end

lambda_cache = joinpath(output_dir, "quartic_lambda_results.jls")

if isfile(lambda_cache)
    log("Using cached λ_self results")
    lambda_result = deserialize(lambda_cache)
else
    geom_ids_lambda = Int[]
    lambda_max_all = Float64[]
    skipped = 0

    for (gid, g) in enumerate(geoms)
        fname = joinpath(DB_DIR,
            @sprintf("h11_%03d", g.h11),
            @sprintf("np_%07d", g.polytope),
            @sprintf("cy_%07d", g.frst),
            "cyax.h5")
        if !isfile(fname)
            skipped += 1; continue
        end
        L = h5read(fname, "cytools/potential/L")
        Q = h5read(fname, "cytools/potential/Q")
        Kinv = h5read(fname, "cytools/geometric/Kinv")
        K = inv(Symmetric(Kinv))
        λs = _compute_lambda_self(L, Q, K)
        λ_max = maximum(filter(isfinite, λs); init=-Inf)
        push!(geom_ids_lambda, gid)
        push!(lambda_max_all, λ_max)
    end

    lambda_result = (geom_ids=geom_ids_lambda, lambda_max=lambda_max_all)
    serialize(lambda_cache, lambda_result)
    log("λ_self done: $(length(geom_ids_lambda)) geometries ($skipped skipped)")
end

# ─────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────
summary_path = joinpath(output_dir, "batch_info.json")
summary = Dict(
    "h11_min" => h11_min,
    "h11_max" => h11_max,
    "batch_label" => batch_label,
    "n_geometries" => length(geoms),
    "n_axions" => length(axions),
    "n_dm_candidates" => length(dm_result.candidates),
    "n_lambda_geometries" => length(lambda_result.geom_ids),
    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
)
open(summary_path, "w") do io
    println(io, "{\n")
    for (k, v) in summary
        println(io, "  \"$k\": $v,")
    end
    println(io, "}")
end

log("="^60)
log("BATCH COMPLETE: h11=$h11_min-$h11_max → $output_dir")
log("  Geometries: $(length(geoms))")
log("  Axions: $(length(axions))")
log("  DM candidates: $(length(dm_result.candidates))")
log("  λ_self computed: $(length(lambda_result.geom_ids))")
log("="^60)


# ═════════════════════════════════════════════════════════════════════
# Dataset building helpers (adapted from src/dataset_builders.jl)
# ═════════════════════════════════════════════════════════════════════

function _build_axion_dataset(h11list)
    n_total = size(h11list, 2)
    all_axions = NamedTuple{
        (:h11, :polytope, :frst, :axion_idx, :log10m, :log10f, :log10fK, :log10_abundance),
        Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64}
    }[]
    failed = Tuple{Int,Int,Int,String}[]
    for col_idx in 1:n_total
        col = h11list[:, col_idx]
        geom_idx = CYAxiverse.structs.GeometryIndex(col...)
        if col_idx % 20 == 0 || col_idx == 1
            log("  Dataset: $col_idx/$n_total (h11=$(geom_idx.h11), np=$(geom_idx.polytope), cy=$(geom_idx.frst))")
        end
        try
            spectrum = CYAxiverse.generate.pq_spectrum(geom_idx)
            ms = spectrum.m
            fs = spectrum.f
            fKs = spectrum.fK
            for i in 1:length(ms)
                log10_abundance = _abundance_ratio(ms[i], fs[i]) + log10(0.12)
                push!(all_axions, (
                    h11=geom_idx.h11, polytope=geom_idx.polytope, frst=geom_idx.frst,
                    axion_idx=i, log10m=ms[i], log10f=fs[i], log10fK=fKs[i],
                    log10_abundance=log10_abundance
                ))
            end
        catch e
            msg = string(e)
            push!(failed, (geom_idx.h11, geom_idx.polytope, geom_idx.frst, msg))
        end
    end
    log("Dataset built: $(length(all_axions)) axions, $(length(failed)) failed")
    return (axions=all_axions, failed=failed)
end

function _build_geometry_dataset(h11list)
    n_total = size(h11list, 2)
    all_geoms = NamedTuple{
        (:h11, :polytope, :frst, :n_axions, :cy_volume, :h21, :n_light_axions),
        Tuple{Int,Int,Int,Int,Float64,Int,Int}
    }[]
    failed = Tuple{Int,Int,Int,String}[]
    for col_idx in 1:n_total
        col = h11list[:, col_idx]
        geom_idx = CYAxiverse.structs.GeometryIndex(col...)
        if col_idx % 20 == 0 || col_idx == 1
            log("  Geom data: $col_idx/$n_total (h11=$(geom_idx.h11))")
        end
        try
            geo_data = CYAxiverse.read.geometry(geom_idx)
            pq = CYAxiverse.generate.pq_spectrum(geom_idx)
            n_light = sum(pq.m .< -27.0)
            push!(all_geoms, (
                h11=geom_idx.h11, polytope=geom_idx.polytope, frst=geom_idx.frst,
                n_axions=geom_idx.h11, cy_volume=geo_data.cy_volume,
                h21=geo_data.h21, n_light_axions=n_light
            ))
        catch e
            msg = string(e)
            push!(failed, (geom_idx.h11, geom_idx.polytope, geom_idx.frst, msg))
        end
    end
    log("Geometry dataset: $(length(all_geoms)) geometries, $(length(failed)) failed")
    return (geometries=all_geoms, failed=failed)
end

function _abundance_ratio(log10m, log10f; θ=nothing)
    log10_θ = θ === nothing ? 0.0 : 2.0 * log10(θ / π)
    H_eq = -27.0
    if log10m > H_eq
        return log10_θ + 1.5*(log10f - 12.0) - 0.5*(log10m + 22.0)
    else
        return log10_θ + 2.0*(log10f - 12.0) + 0.5*(log10m + 22.0)
    end
end
