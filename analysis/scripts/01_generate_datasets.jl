include(joinpath(@__DIR__, "00_setup.jl"))

log_message("Stage 1: Dataset Generation")

full_run = "--full" in ARGS
h11list_all = CYAxiverse.filestructure.paths_cy()[2]
n_total = size(h11list_all, 2)
log_message("Total geometries in database: $n_total")
log_message("h11 range in DB: $(minimum(h11list_all[1,:])) to $(maximum(h11list_all[1,:]))")

if full_run
    h11list = h11list_all
    log_message("Mode: FULL — processing ALL geometries")
    cache_label = ""
else
    h11_min = 4
    h11_max = 10
    mask = h11_min .<= h11list_all[1, :] .<= h11_max
    h11list = h11list_all[:, mask]
    log_message("Mode: TEST — processing h11=$h11_min to h11=$h11_max ($(size(h11list,2)) geometries)")
    cache_label = "_test"
end

axion_cache = "axion_dataset$(cache_label).jls"
geom_cache = "geometry_dataset$(cache_label).jls"

if cache_exists(axion_cache) && cache_exists(geom_cache)
    log_message("Cached datasets found. Skipping generation.")
    log_message("  $axion_cache: $(format_bytes_filesize(joinpath(OUTPUT_DIR, axion_cache)))")
    log_message("  $geom_cache: $(format_bytes_filesize(joinpath(OUTPUT_DIR, geom_cache)))")
else
    axion_data = build_axion_dataset(h11list; verbose=true)
    save_dataset(axion_data, axion_cache)
    if !isempty(axion_data.failed)
        save_dataset(axion_data.failed, "failed_geometries_axion$(cache_label).jls")
        log_message("Failed geometries (axion): $(length(axion_data.failed))")
    end
    geometry_data = build_geometry_dataset(h11list; verbose=true)
    save_dataset(geometry_data, geom_cache)
    if !isempty(geometry_data.failed)
        save_dataset(geometry_data.failed, "failed_geometries_geometry$(cache_label).jls")
        log_message("Failed geometries (geometry): $(length(geometry_data.failed))")
    end
end

log_message("Dataset generation complete.")
log_message("Output files:")
for f in readdir(OUTPUT_DIR)
    fp = joinpath(OUTPUT_DIR, f)
    isfile(fp) && log_message("  $f ($(format_bytes_filesize(fp)))")
end
log_message("="^60)
