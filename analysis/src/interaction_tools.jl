function analyze_interactions(interaction_data)
    lself_vals = [d.log10_lambda_self for d in interaction_data]
    n = length(lself_vals)
    median_lself = median(lself_vals)
    max_lself = maximum(lself_vals)
    min_lself = minimum(lself_vals)
    mean_lself = mean(lself_vals)
    std_lself = std(lself_vals)
    strong = sort(interaction_data; by=d -> d.log10_lambda_self, rev=true)
    n_strong = min(10, n)
    strongest = strong[1:n_strong]
    return (n=n, median=median_lself, max=max_lself, min=min_lself,
            mean=mean_lself, std=std_lself, strongest=strongest)
end

function print_interaction_summary(name, stats, label="log10(λ_self)")
    println("\n"^2, "="^60)
    println("  $name")
    println("="^60)
    println("  n        = $(stats.n)")
    println("  mean     = $(round(stats.mean, digits=4))")
    println("  std      = $(round(stats.std, digits=4))")
    println("  median   = $(round(stats.median, digits=4))")
    println("  max      = $(round(stats.max, digits=4))")
    println("  min      = $(round(stats.min, digits=4))")
    println("="^60)
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
            m_vals = hp_data["m"]
            fK = hp_data["fK"]
            fpert = hp_data["fpert"]
            lself = hp_data["λself"]
            n = length(m_vals)
            for i in 1:n
                push!(all_interactions, (
                    h11=geom_idx.h11, polytope=geom_idx.polytope, frst=geom_idx.frst,
                    axion_idx=i, log10m=m_vals[i], log10fK=fK[i],
                    log10_fpert=fpert[i], log10_lambda_self=lself[i]))
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

function collect_abundance_vs_lambda(axion_data, interaction_data)
    axion_dict = Dict(
        (d.h11, d.polytope, d.frst, d.axion_idx) => d
        for d in axion_data
    )
    combined = NamedTuple{(:h11, :polytope, :frst, :axion_idx, :log10m, :log10f, :log10_abundance, :log10_lambda_self), Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64}}[]
    for ia in interaction_data
        key = (ia.h11, ia.polytope, ia.frst, ia.axion_idx)
        if haskey(axion_dict, key)
            ad = axion_dict[key]
            push!(combined, (
                h11=ad.h11, polytope=ad.polytope, frst=ad.frst,
                axion_idx=ad.axion_idx, log10m=ad.log10m, log10f=ad.log10f,
                log10_abundance=ad.log10_abundance,
                log10_lambda_self=ia.log10_lambda_self))
        end
    end
    log_message("Combined abundance-lambda dataset: $(length(combined)) axions")
    return combined
end

function geometry_interaction_stats(interaction_data)
    geom_keys = unique([(d.h11, d.polytope, d.frst) for d in interaction_data])
    results = []
    for gk in geom_keys
        vals = [d.log10_lambda_self for d in interaction_data
                if d.h11 == gk[1] && d.polytope == gk[2] && d.frst == gk[3]]
        push!(results, (
            h11=gk[1], polytope=gk[2], frst=gk[3],
            median_lambda=median(vals), max_lambda=maximum(vals),
            mean_lambda=mean(vals), n_axions=length(vals)))
    end
    return results
end
