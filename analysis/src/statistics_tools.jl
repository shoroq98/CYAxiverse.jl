function extract_vector(data, field)
    [getindex(d, field) for d in data]
end

function quantile_manual(vals, p)
    n = length(vals)
    sorted = sort(vals)
    idx = max(1, min(n, ceil(Int, p * n)))
    return sorted[idx]
end

function summary_statistics(values)
    n = length(values)
    m = mean(values)
    s = std(values)
    q = [quantile_manual(values, p) for p in [0.05, 0.25, 0.50, 0.75, 0.95]]
    return (n=n, mean=m, std=s, median=q[3], min=minimum(values),
            max=maximum(values), q05=q[1], q25=q[2], q75=q[4], q95=q[5])
end

function print_summary(name, stats)
    println("\n"^2, "="^60)
    println("  $name")
    println("="^60)
    println("  n        = $(stats.n)")
    println("  mean     = $(round(stats.mean, digits=4))")
    println("  std      = $(round(stats.std, digits=4))")
    println("  median   = $(round(stats.median, digits=4))")
    println("  min      = $(round(stats.min, digits=4))")
    println("  max      = $(round(stats.max, digits=4))")
    println("  5%%       = $(round(stats.q05, digits=4))")
    println("  25%%      = $(round(stats.q25, digits=4))")
    println("  75%%      = $(round(stats.q75, digits=4))")
    println("  95%%      = $(round(stats.q95, digits=4))")
    println("="^60)
end

function build_feature_matrix(data, fields)
    n = length(data)
    p = length(fields)
    X = zeros(n, p)
    for (j, f) in enumerate(fields)
        X[:, j] = extract_vector(data, f)
    end
    return X, fields
end

function correlation_table(data, fields)
    X, _ = build_feature_matrix(data, fields)
    C = cor(X)
    println("\n"^2, "="^60)
    println("  Correlation Matrix")
    println("="^60)
    print("  $(rpad("", 16))")
    for f in fields
        print("$(rpad(f, 16))")
    end
    println()
    for i in 1:length(fields)
        print("  $(rpad(fields[i], 16))")
        for j in 1:length(fields)
            print("$(rpad(round(C[i,j], digits=4), 16))")
        end
        println()
    end
    println("="^60)
    return C
end

function run_pca(data, fields)
    X, _ = build_feature_matrix(data, fields)
    n, p = size(X)
    Xnorm = (X .- mean(X, dims=1)) ./ std(X, dims=1)
    C = cov(Xnorm)
    eig = eigen(C)
    perm = sortperm(eig.values, rev=true)
    eigenvalues = eig.values[perm]
    eigenvectors = eig.vectors[:, perm]
    explained_var = eigenvalues / sum(eigenvalues)
    explained_var_cumul = cumsum(explained_var)
    projections = Xnorm * eigenvectors
    result = (eigenvalues=eigenvalues, eigenvectors=eigenvectors,
              explained_var=explained_var, explained_var_cumul=explained_var_cumul,
              projections=projections, fields=fields)
    println("\n"^2, "="^60)
    println("  PCA Results")
    println("="^60)
    for i in 1:p
        println("  PC$(i): λ=$(round(eigenvalues[i], digits=4)), "
                * "var=$(round(100*explained_var[i], digits=2))%, "
                * "cumul=$(round(100*explained_var_cumul[i], digits=2))%")
    end
    println("\n  Loadings:")
    for i in 1:min(p, 3)
        println("  PC$(i): " * join(["$(fields[j])=$(round(eigenvectors[j,i], digits=4))" for j in 1:p], ", "))
    end
    println("="^60)
    return result
end

function binned_statistic(data, group_field, value_field; nbins=10)
    group_vals = extract_vector(data, group_field)
    value_vals = extract_vector(data, value_field)
    bin_edges = range(minimum(group_vals), maximum(group_vals); length=nbins+1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in 1:nbins]
    binned = [Float64[] for _ in 1:nbins]
    for i in eachindex(group_vals)
        for b in 1:nbins
            if bin_edges[b] <= group_vals[i] < bin_edges[b+1] || (b == nbins && group_vals[i] == bin_edges[end])
                push!(binned[b], value_vals[i])
                break
            end
        end
    end
    bin_means = [isempty(b) ? NaN : mean(b) for b in binned]
    bin_stds = [isempty(b) ? NaN : std(b) for b in binned]
    bin_counts = [length(b) for b in binned]
    return (bin_centers=bin_centers, bin_means=bin_means, bin_stds=bin_stds,
            bin_counts=bin_counts, bin_edges=bin_edges)
end

function distribution_analysis(values)
    n = length(values)
    m = mean(values)
    s = std(values)
    med = median(values)
    q = [quantile_manual(values, p) for p in [0.05, 0.25, 0.50, 0.75, 0.95]]
    skew = sum((values .- m).^3) / (n * s^3)
    kurt = sum((values .- m).^4) / (n * s^4) - 3.0
    return (n=n, mean=m, std=s, median=med, min=minimum(values), max=maximum(values),
            q05=q[1], q25=q[2], q75=q[4], q95=q[5], skewness=skew, kurtosis=kurt)
end

function print_distribution(name, d)
    println("\n"^2, "="^60)
    println("  $name")
    println("="^60)
    println("  n        = $(d.n)")
    println("  mean     = $(round(d.mean, digits=4))")
    println("  std      = $(round(d.std, digits=4))")
    println("  median   = $(round(d.median, digits=4))")
    println("  skewness = $(round(d.skewness, digits=4))")
    println("  kurtosis = $(round(d.kurtosis, digits=4))")
    println("  [5% , 95%] = [$(round(d.q05, digits=2)), $(round(d.q95, digits=2))]")
    println("  range    = [$(round(d.min, digits=2)), $(round(d.max, digits=2))]")
    println("="^60)
end

function fit_scaling_law(x, y)
    X = hcat(ones(length(x)), x)
    β = X \ y
    a, b = β[1], β[2]
    y_pred = X * β
    residuals = y - y_pred
    r2 = 1.0 - sum(residuals.^2) / sum((y .- mean(y)).^2)
    return (intercept=a, slope=b, r2=r2, y_pred=y_pred)
end

function print_scaling_law(name, fit)
    println("\n"^2, "="^60)
    println("  Scaling Law Fit: $name")
    println("="^60)
    println("  log10(Ω) = $(round(fit.intercept, digits=4)) + $(round(fit.slope, digits=4)) × log10(X)")
    println("  R² = $(round(fit.r2, digits=4))")
    println("="^60)
end

function outlier_analysis(values; n_sigma=3)
    m = mean(values)
    s = std(values)
    outliers = findall(abs.(values .- m) .> n_sigma * s)
    return (n_outliers=length(outliers), outlier_indices=outliers,
            threshold_low=m - n_sigma*s, threshold_high=m + n_sigma*s,
            n_sigma=n_sigma)
end

function population_separation(axion_data; mass_cuts=[-27.0, 10.0])
    ms = extract_vector(axion_data, :log10m)
    ultra_light = [axion_data[i] for i in eachindex(ms) if ms[i] < mass_cuts[1]]
    light = [axion_data[i] for i in eachindex(ms) if mass_cuts[1] <= ms[i] < mass_cuts[2]]
    heavy = [axion_data[i] for i in eachindex(ms) if ms[i] >= mass_cuts[2]]
    log_message("Population separation at log10m=$mass_cuts: ultra-light=$(length(ultra_light)), light=$(length(light)), heavy=$(length(heavy))")
    return (ultra_light=ultra_light, light=light, heavy=heavy)
end

function print_outlier_summary(name, outliers, total_n)
    println("\n"^2, "="^60)
    println("  Outlier Analysis: $name")
    println("="^60)
    println("  Total points: $total_n")
    println("  Outliers: $(outliers.n_outliers) ($(round(100*outliers.n_outliers/total_n, digits=2))%)")
    println("  Threshold: [$(round(outliers.threshold_low, digits=2)), $(round(outliers.threshold_high, digits=2))]")
    println("  Sigma: $(outliers.n_sigma)")
    println("="^60)
end

function conditional_statistics_abundance_by_h11(data)
    h11_vals = extract_vector(data, :h11)
    abundance_vals = extract_vector(data, :log10_abundance)
    h11s = sort(unique(h11_vals))
    println("\n"^2, "="^60)
    println("  Abundance by h11")
    println("="^60)
    println("  $(rpad("h11", 6)) $(rpad("count", 8)) $(rpad("mean", 12)) $(rpad("std", 12)) $(rpad("median", 12))")
    results = []
    for h in h11s
        idx = h11_vals .== h
        vals = abundance_vals[idx]
        if length(vals) > 0
            m = mean(vals)
            s = std(vals)
            med = median(vals)
            println("  $(rpad(h, 6)) $(rpad(length(vals), 8)) $(rpad(round(m, digits=4), 12)) $(rpad(round(s, digits=4), 12)) $(rpad(round(med, digits=4), 12))")
            push!(results, (h11=h, count=length(vals), mean=m, std=s, median=med))
        end
    end
    println("="^60)
    return results
end

function count_light_axions_per_geometry(axion_data)
    geoms = unique([(d.h11, d.polytope, d.frst) for d in axion_data])
    println("\n"^2, "="^60)
    println("  Light Axions per Geometry (m < 10^-27 eV)")
    println("="^60)
    println("  $(rpad("h11", 6)) $(rpad("total_geoms", 14)) $(rpad("mean_n_light", 14)) $(rpad("median_n_light", 14))")
    results = []
    for h in sort(unique([g[1] for g in geoms]))
        geoms_h = [g for g in geoms if g[1] == h]
        n_light_list = Int[]
        for g in geoms_h
            axions_for_geom = [d for d in axion_data if d.h11 == g[1] && d.polytope == g[2] && d.frst == g[3]]
            n_light = count(d -> d.log10m < -27.0, axions_for_geom)
            push!(n_light_list, n_light)
        end
        m = mean(n_light_list)
        med = median(n_light_list)
        println("  $(rpad(h, 6)) $(rpad(length(geoms_h), 14)) $(rpad(round(m, digits=2), 14)) $(rpad(round(med, digits=2), 14))")
        push!(results, (h11=h, n_geoms=length(geoms_h), mean_n_light=m, median_n_light=med))
    end
    println("="^60)
    return results
end

function save_statistics(stats, filename)
    path = joinpath(OUTPUT_DIR, filename)
    serialize(path, stats)
    log_message("Saved statistics: $filename")
    return path
end

function load_statistics(filename)
    path = joinpath(OUTPUT_DIR, filename)
    if !isfile(path)
        error("Statistics not found: $path. Run 02_statistical_analysis.jl first.")
    end
    return deserialize(path)
end
