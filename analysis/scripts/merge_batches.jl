using Serialization
using Printf
using Dates

# ─────────────────────────────────────────────────────────────────────
# CLI args
# ─────────────────────────────────────────────────────────────────────
input_base = "/mnt/results"
output_base = "/mnt/results"

for (i, arg) in enumerate(ARGS)
    if arg == "--input" && i < length(ARGS)
        input_base = ARGS[i+1]
    elseif arg == "--output" && i < length(ARGS)
        output_base = ARGS[i+1]
    end
end

batch_dir = joinpath(input_base, "batches")
@assert isdir(batch_dir) "Batch directory not found: $batch_dir"

function log(msg)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    println("[$ts] $msg")
end

log("="^60)
log("Merging batch results")
log("Input: $batch_dir")
log("Output: $output_base")
log("="^60)

batch_names = sort(readdir(batch_dir))
log("Found $(length(batch_names)) batches")

# ─────────────────────────────────────────────────────────────────────
# Merge axion datasets
# ─────────────────────────────────────────────────────────────────────
all_axions = []
all_ax_failed = []
axion_total = 0

for bname in batch_names
    f = joinpath(batch_dir, bname, "axion_dataset.jls")
    if !isfile(f); log("  WARNING: no axion data in $bname"); continue; end
    data = deserialize(f)
    append!(all_axions, data.axions)
    append!(all_ax_failed, data.failed)
    axion_total += length(data.axions)
end

serialize(joinpath(output_base, "axion_dataset.jls"),
    (axions=all_axions, failed=all_ax_failed))
log("Merged axions: $axion_total (across $(length(all_ax_failed)) failures)")

# ─────────────────────────────────────────────────────────────────────
# Merge geometry datasets
# ─────────────────────────────────────────────────────────────────────
all_geoms = []
all_geom_failed = []
geom_total = 0

for bname in batch_names
    f = joinpath(batch_dir, bname, "geometry_dataset.jls")
    if !isfile(f); log("  WARNING: no geometry data in $bname"); continue; end
    data = deserialize(f)
    append!(all_geoms, data.geometries)
    append!(all_geom_failed, data.failed)
    geom_total += length(data.geometries)
end

serialize(joinpath(output_base, "geometry_dataset.jls"),
    (geometries=all_geoms, failed=all_geom_failed))
log("Merged geometries: $geom_total (across $(length(all_geom_failed)) failures)")

# ─────────────────────────────────────────────────────────────────────
# Merge DM ratio results
# ─────────────────────────────────────────────────────────────────────
all_geom_ids = Int[]
all_log10_ratios = Float64[]
all_thetas = Float64[]
all_candidates = []

for bname in batch_names
    f = joinpath(batch_dir, bname, "dm_density_ratio_results.jls")
    if !isfile(f); continue; end
    data = deserialize(f)
    append!(all_geom_ids, data.geom_ids_ratio)
    append!(all_log10_ratios, data.log10_total_ratios)
    append!(all_thetas, data.theta_used)
    append!(all_candidates, data.candidates)
end

serialize(joinpath(output_base, "dm_density_ratio_results.jls"),
    (geom_ids_ratio=all_geom_ids, log10_total_ratios=all_log10_ratios,
     theta_used=all_thetas, candidates=all_candidates))
log("Merged DM results: $(length(all_geom_ids)) entries, $(length(all_candidates)) candidates")

# ─────────────────────────────────────────────────────────────────────
# Merge λ_self results
# ─────────────────────────────────────────────────────────────────────
all_lambda_geom_ids = Int[]
all_lambda_max = Float64[]

for bname in batch_names
    f = joinpath(batch_dir, bname, "quartic_lambda_results.jls")
    if !isfile(f); continue; end
    data = deserialize(f)
    append!(all_lambda_geom_ids, data.geom_ids)
    append!(all_lambda_max, data.lambda_max)
end

serialize(joinpath(output_base, "quartic_lambda_results.jls"),
    (geom_ids=all_lambda_geom_ids, lambda_max=all_lambda_max))
log("Merged λ_self results: $(length(all_lambda_geom_ids)) geometries")

# ─────────────────────────────────────────────────────────────────────
# Merge batch logs into a single pipeline log
# ─────────────────────────────────────────────────────────────────────
log_path = joinpath(output_base, "pipeline_merge.log")
open(log_path, "w") do io
    println(io, "="^60)
    println(io, "  Pipeline merge completed")
    println(io, "  $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "="^60)
    println(io)
    println(io, "Batches processed: $(length(batch_names))")
    println(io, "Total axions: $axion_total")
    println(io, "Total geometries: $geom_total")
    println(io, "DM candidates: $(length(all_candidates))")
    println(io, "λ_self geometries: $(length(all_lambda_geom_ids))")
end

log("="^60)
log("Merge complete!")
log("  Batches: $(length(batch_names))")
log("  Axions: $axion_total")
log("  Geometries: $geom_total")
log("  DM candidates: $(length(all_candidates))")
log("  λ_self: $(length(all_lambda_geom_ids)) geometries")
log("="^60)
