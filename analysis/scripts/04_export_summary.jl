include(joinpath(@__DIR__, "00_setup.jl"))

log_message("Stage 4: Summary Export")

cache_label = "--full" in ARGS ? "" : "_test"
axion_data = load_dataset("axion_dataset$(cache_label).jls")
geometry_data = load_dataset("geometry_dataset$(cache_label).jls")
stats = load_dataset("summary_statistics.jls")

axions = axion_data.axions
geometries = geometry_data.geometries

summary_path = joinpath(OUTPUT_DIR, "pipeline_summary.txt")
open(summary_path, "w") do io
    println(io, "="^60)
    println(io, "  CYAxiverse Analysis Pipeline Summary")
    println(io, "  Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "="^60)
    println(io)

    println(io, "Dataset Overview:")
    println(io, "  Total axions: $(length(axions))")
    println(io, "  Total geometries: $(length(geometries))")
    println(io, "  Failed geometries (axion): $(length(axion_data.failed))")
    println(io, "  Failed geometries (geometry): $(length(geometry_data.failed))")

    h11s = unique([d.h11 for d in axions])
    println(io, "  h11 range: $(minimum(h11s)) - $(maximum(h11s))")
    println(io)

    println(io, "Summary Statistics:")
    for (field, s) in stats.summary
        println(io, "  $field:")
        println(io, "    n = $(s.n), mean = $(round(s.mean, digits=4)), std = $(round(s.std, digits=4))")
        println(io, "    min = $(round(s.min, digits=4)), max = $(round(s.max, digits=4))")
    end
    println(io)

    println(io, "PCA:")
    for i in eachindex(stats.pca.explained_var)
        println(io, "  PC$i: $(round(100*stats.pca.explained_var[i], digits=2))% variance explained")
    end
    println(io)

    println(io, "Output Files:")
    println(io, "  Datasets:")
    for f in sort(readdir(OUTPUT_DIR))
        fp = joinpath(OUTPUT_DIR, f)
        isfile(fp) && println(io, "    $f ($(format_bytes_filesize(fp)))")
    end
    println(io, "  Plots:")
    for f in sort(readdir(PLOT_DIR))
        fp = joinpath(PLOT_DIR, f)
        isfile(fp) && println(io, "    $f ($(format_bytes_filesize(fp)))")
    end
end

log_message("Summary exported to: $summary_path")
println("\n" * read(summary_path, String))
log_message("="^60)
log_message("Pipeline complete.")
