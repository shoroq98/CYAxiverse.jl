include(joinpath(@__DIR__, "00_setup.jl"))

log_message("Stage 2: Statistical Analysis")

cache_label = "--full" in ARGS ? "" : "_test"
axion_data = load_dataset("axion_dataset$(cache_label).jls")
geometry_data = load_dataset("geometry_dataset$(cache_label).jls")
axions = axion_data.axions
geometries = geometry_data.geometries

log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

println("\n" * "="^60)
println("  1. Summary Statistics")
println("="^60)

axion_fields = [
    (:log10m, "log10(mass / eV)"),
    (:log10f, "log10(decay constant / GeV)"),
    (:log10_abundance, "log10(Ω_a h²)"),
]

summary_results = Dict{Symbol, Any}()
for (field, label) in axion_fields
    vals = extract_vector(axions, field)
    stats = summary_statistics(vals)
    summary_results[field] = stats
    print_summary(label, stats)
end

geometry_fields = [
    (:h11, "h11"),
    (:cy_volume, "log10(CY Volume)"),
]

for (field, label) in geometry_fields
    vals = extract_vector(geometries, field)
    if field == :cy_volume
        vals = log10.(vals)
    end
    stats = summary_statistics(vals)
    summary_results[field] = stats
    print_summary(label, stats)
end

println("\n" * "="^60)
println("  2. Correlation Matrix")
println("="^60)

corr_fields = [:log10m, :log10f, :log10_abundance]
C_axion = correlation_table(axions, corr_fields)

println("\n" * "="^60)
println("  3. Conditional Statistics: Abundance by Mass Bins")
println("="^60)
mass_bins = binned_statistic(axions, :log10m, :log10_abundance; nbins=15)
println("  Mass bins (log10 eV) | Mean abundance | Std | Count")
for i in eachindex(mass_bins.bin_centers)
    isfinite(mass_bins.bin_means[i]) || continue
    println("  $(round(mass_bins.bin_centers[i], digits=2)) ± $(round((mass_bins.bin_edges[i+1]-mass_bins.bin_edges[i])/2, digits=2)) | " *
            "$(round(mass_bins.bin_means[i], digits=4)) | $(round(mass_bins.bin_stds[i], digits=4)) | $(mass_bins.bin_counts[i])")
end

println("\n" * "="^60)
println("  4. Conditional Statistics: Abundance by h11")
println("="^60)
abund_by_h11 = conditional_statistics_abundance_by_h11(axions)

println("\n" * "="^60)
println("  5. Light Axions per Geometry")
println("="^60)
light_axion_stats = count_light_axions_per_geometry(axions)

println("\n" * "="^60)
println("  6. PCA")
println("="^60)
pca_result = run_pca(axions, [:log10m, :log10f, :log10_abundance])

println("\n" * "="^60)
println("  Saving Results")
println("="^60)

all_statistics = (summary=summary_results, correlation=C_axion, pca=pca_result,
                  mass_bins=mass_bins, abund_by_h11=abund_by_h11,
                  light_axion_stats=light_axion_stats,
                  n_axions=length(axions), n_geometries=length(geometries))

save_dataset(all_statistics, "summary_statistics.jls")
save_dataset(C_axion, "correlation_matrix.jls")
save_dataset(pca_result, "pca_results.jls")

log_message("="^60)
log_message("Statistical analysis complete.")
