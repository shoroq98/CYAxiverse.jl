if !@isdefined(OUTPUT_DIR)
    PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
    const OUTPUT_DIR = joinpath(PROJECT_ROOT, "outputs")
    const PLOT_DIR = joinpath(PROJECT_ROOT, "plots")
    const LOG_DIR = joinpath(PROJECT_ROOT, "logs")
    const LOG_FILE = joinpath(LOG_DIR, "run_log.txt")
    for d in [OUTPUT_DIR, PLOT_DIR, LOG_DIR]
        isdir(d) || mkpath(d)
    end
end

function log_message(msg)
    open(LOG_FILE, "a") do io
        println(io, "[$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))] ", msg)
    end
    @info msg
end

function extract_axion_data(geom_idx, spectrum; θ=nothing)
    h11, polytope, frst = geom_idx.h11, geom_idx.polytope, geom_idx.frst
    ms = spectrum.m
    fs = spectrum.f
    fKs = spectrum.fK
    n_axions = length(ms)
    records = Vector{NamedTuple{(:h11, :polytope, :frst, :axion_idx, :log10m, :log10f, :log10fK, :log10_abundance), Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64}}}(undef, n_axions)
    for i in 1:n_axions
        log10m = ms[i]
        log10f = fs[i]
        log10fK = fKs[i]
        log10_abundance = log10_axion_density_ratio(log10m, log10f; θ=θ) + log10(0.12)
        records[i] = (h11=h11, polytope=polytope, frst=frst, axion_idx=i,
                      log10m=log10m, log10f=log10f, log10fK=log10fK,
                      log10_abundance=log10_abundance)
    end
    return records
end

function extract_geometry_data(geom_idx)
    h11, polytope, frst = geom_idx.h11, geom_idx.polytope, geom_idx.frst
    geo_data = CYAxiverse.read.geometry(geom_idx)
    n_light = sum(CYAxiverse.generate.pq_spectrum(geom_idx).m .< -27.0)
    return (h11=h11, polytope=polytope, frst=frst, n_axions=h11,
            cy_volume=geo_data.cy_volume, h21=geo_data.h21,
            n_light_axions=n_light)
end

function build_axion_dataset(h11list; max_geoms=typemax(Int), verbose=true)
    n_total = min(size(h11list, 2), max_geoms)
    all_axions = NamedTuple{(:h11, :polytope, :frst, :axion_idx, :log10m, :log10f, :log10fK, :log10_abundance), Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64}}[]
    failed = Tuple{Int,Int,Int,String}[]
    for col_idx in 1:n_total
        col = h11list[:, col_idx]
        geom_idx = CYAxiverse.structs.GeometryIndex(col...)
        if verbose
            @printf("Processing axion data: geometry %d / %d (h11=%d, polytope=%d, frst=%d)\n",
                    col_idx, n_total, geom_idx.h11, geom_idx.polytope, geom_idx.frst)
            flush(stdout)
        end
        try
            spectrum = CYAxiverse.generate.pq_spectrum(geom_idx)
            records = extract_axion_data(geom_idx, spectrum)
            append!(all_axions, records)
        catch e
            msg = "Error: $(e)"
            @warn "Skipping geometry ($(geom_idx.h11), $(geom_idx.polytope), $(geom_idx.frst)): $msg"
            push!(failed, (geom_idx.h11, geom_idx.polytope, geom_idx.frst, msg))
        end
    end
    log_message("Axion dataset: $(length(all_axions)) axions from $(n_total) geometries, $(length(failed)) failed")
    return (axions=all_axions, failed=failed)
end

function build_geometry_dataset(h11list; max_geoms=typemax(Int), verbose=true)
    n_total = min(size(h11list, 2), max_geoms)
    all_geoms = NamedTuple{(:h11, :polytope, :frst, :n_axions, :cy_volume, :h21, :n_light_axions), Tuple{Int,Int,Int,Int,Float64,Int,Int}}[]
    failed = Tuple{Int,Int,Int,String}[]
    for col_idx in 1:n_total
        col = h11list[:, col_idx]
        geom_idx = CYAxiverse.structs.GeometryIndex(col...)
        if verbose
            @printf("Processing geometry data: %d / %d (h11=%d, polytope=%d, frst=%d)\n",
                    col_idx, n_total, geom_idx.h11, geom_idx.polytope, geom_idx.frst)
            flush(stdout)
        end
        try
            geom_record = extract_geometry_data(geom_idx)
            push!(all_geoms, geom_record)
        catch e
            msg = "Error: $(e)"
            @warn "Skipping geometry ($(geom_idx.h11), $(geom_idx.polytope), $(geom_idx.frst)): $msg"
            push!(failed, (geom_idx.h11, geom_idx.polytope, geom_idx.frst, msg))
        end
    end
    log_message("Geometry dataset: $(length(all_geoms)) geometries, $(length(failed)) failed")
    return (geometries=all_geoms, failed=failed)
end

function build_interaction_dataset(h11list; max_geoms=typemax(Int), verbose=true)
    n_total = min(size(h11list, 2), max_geoms)
    all_interactions = NamedTuple{(:h11, :polytope, :frst, :axion_idx, :log10m, :log10fK, :log10_fpert, :log10_lambda_self), Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64}}[]
    failed = Tuple{Int,Int,Int,String}[]
    for col_idx in 1:n_total
        col = h11list[:, col_idx]
        geom_idx = CYAxiverse.structs.GeometryIndex(col...)
        if verbose
            @printf("Processing HP spectrum: %d / %d (h11=%d, polytope=%d, frst=%d)\n",
                    col_idx, n_total, geom_idx.h11, geom_idx.polytope, geom_idx.frst)
            flush(stdout)
        end
        try
            hp_data = CYAxiverse.generate.hp_spectrum(geom_idx)
            msign = hp_data["msign"]
            m_vals = hp_data["m"]
            fK = hp_data["fK"]
            fpert = hp_data["fpert"]
            lself = hp_data["λself"]
            n = length(m_vals)
            for i in 1:n
                push!(all_interactions, (
                    h11=geom_idx.h11, polytope=geom_idx.polytope, frst=geom_idx.frst,
                    axion_idx=i,
                    log10m=m_vals[i],
                    log10fK=fK[i],
                    log10_fpert=fpert[i],
                    log10_lambda_self=lself[i]
                ))
            end
        catch e
            msg = "Error: $(e)"
            @warn "Skipping HP spectrum ($(geom_idx.h11), $(geom_idx.polytope), $(geom_idx.frst)): $msg"
            push!(failed, (geom_idx.h11, geom_idx.polytope, geom_idx.frst, msg))
        end
    end
    log_message("Interaction dataset: $(length(all_interactions)) axions from $(n_total) geometries, $(length(failed)) failed")
    return (interactions=all_interactions, failed=failed)
end

function cache_exists(name)
    path = joinpath(OUTPUT_DIR, name)
    return isfile(path)
end

function save_dataset(data, filename)
    path = joinpath(OUTPUT_DIR, filename)
    serialize(path, data)
    log_message("Saved dataset: $filename ($(format_bytes_filesize(path)))")
    return path
end

function load_dataset(filename)
    path = joinpath(OUTPUT_DIR, filename)
    if !isfile(path)
        error("Dataset not found: $path. Run 01_generate_datasets.jl first.")
    end
    return deserialize(path)
end

function format_bytes_filesize(path)
    bytes = filesize(path)
    if bytes < 1024
        return "$bytes B"
    elseif bytes < 1024^2
        return "$(round(bytes / 1024, digits=1)) KB"
    elseif bytes < 1024^3
        return "$(round(bytes / 1024^2, digits=1)) MB"
    else
        return "$(round(bytes / 1024^3, digits=1)) GB"
    end
end
