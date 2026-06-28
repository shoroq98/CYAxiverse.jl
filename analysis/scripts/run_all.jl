include(joinpath(@__DIR__, "00_setup.jl"))

# Stage 1: Generate Datasets
log_message("="^60)
log_message("STAGE 1: Dataset Generation")
log_message("="^60)

h11list_all = CYAxiverse.filestructure.paths_cy()[2]
n_total = size(h11list_all, 2)
log_message("Total geometries in database: $n_total")

full_run = "--full" in ARGS
if full_run
    h11list = h11list_all
    log_message("Mode: FULL — processing ALL geometries")
    cache_label = ""
else
    h11_min = 4; h11_max = 20
    mask = h11_min .<= h11list_all[1, :] .<= h11_max
    h11list = h11list_all[:, mask]
    log_message("Mode: DEFAULT — processing h11=$h11_min to $h11_max ($(size(h11list,2)) geometries)")
    cache_label = ""
end

axion_cache = "axion_dataset$(cache_label).jls"
geom_cache = "geometry_dataset$(cache_label).jls"

if cache_exists(axion_cache) && cache_exists(geom_cache)
    log_message("Using cached datasets: $axion_cache, $geom_cache")
else
    axion_data = build_axion_dataset(h11list; verbose=true)
    save_dataset(axion_data, axion_cache)
    geometry_data = build_geometry_dataset(h11list; verbose=true)
    save_dataset(geometry_data, geom_cache)
end

# Stage 2: Statistical Analysis
log_message("="^60)
log_message("STAGE 2: Statistical Analysis")
log_message("="^60)

axion_data = load_dataset(axion_cache)
geometry_data = load_dataset(geom_cache)
axions = axion_data.axions
geometries = geometry_data.geometries
log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

axion_fields = [
    (:log10m, "log10(mass / eV)"),
    (:log10f, "log10(decay constant / GeV)"),
    (:log10_abundance, "log10(Ω_a h²)"),
]
summary_results = Dict{Symbol, Any}()
for (field, label) in axion_fields
    vals = extract_vector(axions, field)
    summary_results[field] = summary_statistics(vals)
    print_summary(label, summary_results[field])
end

for (field, label) in [(:h11, "h11"), (:cy_volume, "log10(CY Volume)")]
    vals = extract_vector(geometries, field)
    if field == :cy_volume; vals = log10.(vals); end
    summary_results[field] = summary_statistics(vals)
    print_summary(label, summary_results[field])
end

C_axion = correlation_table(axions, [:log10m, :log10f, :log10_abundance])
abund_by_h11 = conditional_statistics_abundance_by_h11(axions)
light_axion_stats = count_light_axions_per_geometry(axions)
pca_result = run_pca(axions, [:log10m, :log10f, :log10_abundance])

log_message("Distribution Analysis (skewness, kurtosis)...")
for (field, label) in axion_fields
    vals = extract_vector(axions, field)
    d = distribution_analysis(vals)
    print_distribution(label, d)
    summary_results[Symbol("$(field)_dist")] = d
end

log_message("Scaling Law Fits...")
abund_vals = extract_vector(axions, :log10_abundance)
mass_vals = extract_vector(axions, :log10m)
f_vals = extract_vector(axions, :log10f)
finite_all = isfinite.(abund_vals) .& isfinite.(mass_vals) .& isfinite.(f_vals)

fit_mass = fit_scaling_law(mass_vals[finite_all], abund_vals[finite_all])
print_scaling_law("log10(Ω) vs log10(m)", fit_mass)

fit_f = fit_scaling_law(f_vals[finite_all], abund_vals[finite_all])
print_scaling_law("log10(Ω) vs log10(f)", fit_f)

log_message("Outlier Analysis...")
outlier_abund = outlier_analysis(abund_vals[isfinite.(abund_vals)])
print_outlier_summary("Abundance", outlier_abund, length(abund_vals))

log_message("Population Separation...")
populations = population_separation(axions; mass_cuts=[-27.0, 10.0])

all_statistics = (summary=summary_results, correlation=C_axion, pca=pca_result,
                  abund_by_h11=abund_by_h11, light_axion_stats=light_axion_stats,
                  fit_mass=fit_mass, fit_f=fit_f, outlier_abund=outlier_abund,
                  populations=populations,
                  n_axions=length(axions), n_geometries=length(geometries))

save_dataset(all_statistics, "summary_statistics$(cache_label).jls")
save_dataset(C_axion, "correlation_matrix$(cache_label).jls")
save_dataset(pca_result, "pca_results$(cache_label).jls")

# Stage 3: Plots
log_message("="^60)
log_message("STAGE 3: Plot Generation")
log_message("="^60)

CairoMakie.activate!(type="png")

save_histogram(axions, :log10m, "hist_mass.png"; bins=80,
               title="Axion Mass Distribution", xlabel="log10(mass / eV)")
save_histogram(axions, :log10f, "hist_decay.png"; bins=80,
               title="Axion Decay Constant Distribution", xlabel="log10(decay constant / GeV)")
save_histogram(axions, :log10_abundance, "hist_abundance.png"; bins=80,
               title="Axion Abundance Distribution", xlabel="log10(Ω_a h²)")
save_scatter(axions, :log10m, :log10_abundance, "scatter_mass_abundance.png";
             title="Mass vs Abundance", xlabel="log10(mass / eV)", ylabel="log10(Ω_a h²)")
save_scatter(axions, :log10f, :log10_abundance, "scatter_decay_abundance.png";
             title="Decay Constant vs Abundance", xlabel="log10(decay constant / GeV)", ylabel="log10(Ω_a h²)")
save_scatter_color(axions, :log10m, :log10f, :log10_abundance, "scatter_mass_decay_color.png";
                   title="Mass vs Decay Constant (colored by abundance)",
                   xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)", clabel="log10(Ω_a h²)")
save_heatmap(axions, :log10m, :log10f, "heatmap_mass_decay.png"; nbins_x=30, nbins_y=30,
             title="Mass vs Decay Constant Density", xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)")
save_scatter(axions, :log10m, :log10f, "scatter_mass_decay.png";
             title="Mass vs Decay Constant", xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)")

save_histogram(geometries, :n_light_axions, "hist_light_axions.png"; bins=50,
               title="Light Axions per Geometry", xlabel="Number of light axions (m < 10^-27 eV)")

DM_OBS = 0.12
geom_abundances = Float64[]
geom_h11 = Int[]
geom_n_light = Int[]
for g in geometries
    matching = [d for d in axions if d.h11 == g.h11 && d.polytope == g.polytope && d.frst == g.frst]
    total_abund = isempty(matching) ? -Inf : log10(sum(10 .^ [a.log10_abundance for a in matching]) / DM_OBS)
    push!(geom_abundances, total_abund); push!(geom_h11, g.h11); push!(geom_n_light, g.n_light_axions)
end

finite_abund = geom_abundances[isfinite.(geom_abundances)]
f_geom_abund = Figure(size=(800, 600))
ax_geom_abund = Axis(f_geom_abund[1, 1], title="Total Geometry Abundance Distribution", xlabel="log10(Total Ω_a h²)", ylabel="Count")
hist!(ax_geom_abund, finite_abund, bins=50, color=(:steelblue, 0.7))
save(joinpath(PLOT_DIR, "hist_geometry_abundance.png"), f_geom_abund)

finite_mask = isfinite.(geom_abundances)
f_scatter_geom = Figure(size=(800, 600))
ax_scatter_geom = Axis(f_scatter_geom[1, 1], title="Geometry Abundance vs h11", xlabel="h11", ylabel="log10(Total Ω_a h²)")
scatter!(ax_scatter_geom, geom_h11[finite_mask], geom_abundances[finite_mask], markersize=3, color=(:steelblue, 0.5))
save(joinpath(PLOT_DIR, "scatter_geom_abundance_h11.png"), f_scatter_geom)

f_scatter_DM = Figure(size=(800, 600))
ax_scatter_DM = Axis(f_scatter_DM[1, 1], title="Geometry DM Density Ratio vs h11", xlabel="h11", ylabel="log10(Ωₐ h² / 0.12)")
scatter!(ax_scatter_DM, geom_h11[finite_mask], geom_abundances[finite_mask], markersize=3, color=(:steelblue, 0.5))
hlines!(ax_scatter_DM, [0], color=:black, linestyle=:dash, linewidth=2)
save(joinpath(PLOT_DIR, "scatter_geom_DMratio_h11.png"), f_scatter_DM)

f_abund_h11 = Figure(size=(1200, 800))
ax_abund1 = Axis(f_abund_h11[1, 1], title="Geometry Abundance vs h11", xlabel="h11", ylabel="log10(Total Ω_a h²)")
scatter!(ax_abund1, geom_h11[finite_mask], geom_abundances[finite_mask], markersize=3, color=(:steelblue, 0.5))
ax_abund2 = Axis(f_abund_h11[1, 2], title="Light Axions vs Geometry Abundance", xlabel="n_light_axions", ylabel="log10(Total Ω_a h²)")
scatter!(ax_abund2, geom_n_light[finite_mask], geom_abundances[finite_mask], markersize=3, color=(:steelblue, 0.5))
ax_abund3 = Axis(f_abund_h11[2, 1], title="Light Axions Distribution", xlabel="n_light_axions", ylabel="Count")
hist!(ax_abund3, geom_n_light, bins=50, color=(:steelblue, 0.7))
save(joinpath(PLOT_DIR, "geometry_abundance_h11.png"), f_abund_h11)

pca_file = "pca_results$(cache_label).jls"
if isfile(joinpath(OUTPUT_DIR, pca_file))
    pca_result = deserialize(joinpath(OUTPUT_DIR, pca_file))
    save_pca_plot(pca_result, axions, "pca_projection.png"; title="PCA: Axion Parameter Space", color_field=:log10_abundance)
    save_explained_variance(pca_result, "pca_explained_variance.png"; title="PCA Explained Variance Ratio")
end

log_message("Generating advanced analysis plots...")
save_population_comparison(populations, :log10_abundance, "population_abundance.png";
                           bins=60, title="Abundance by Population",
                           xlabel="log10(Ω_a h²)")
save_population_comparison(populations, :log10m, "population_mass.png";
                           bins=60, title="Mass by Population",
                           xlabel="log10(mass / eV)")

save_scaling_law_plot(mass_vals, abund_vals, fit_mass, "scaling_mass_abundance.png";
                      title="Scaling Law: Abundance vs Mass (R²=$(round(fit_mass.r2, digits=3)))",
                      xlabel="log10(mass / eV)", ylabel="log10(Ω_a h²)")
save_scaling_law_plot(f_vals, abund_vals, fit_f, "scaling_decay_abundance.png";
                      title="Scaling Law: Abundance vs Decay Constant (R²=$(round(fit_f.r2, digits=3)))",
                      xlabel="log10(decay constant / GeV)", ylabel="log10(Ω_a h²)")

log_message("Theta-scan abundance for first 10 geometries...")
theta_scan_results = []
for g in geometries[1:min(10, length(geometries))]
    matching = [d for d in axions if d.h11 == g.h11 && d.polytope == g.polytope && d.frst == g.frst]
    ms = [d.log10m for d in matching]
    fs = [d.log10f for d in matching]
    scan = theta_scan_geometry_density(ms, fs)
    push!(theta_scan_results, (geom=g, scan=scan))
end
save_dataset(theta_scan_results, "theta_scan_results$(cache_label).jls")
log_message("Theta scan saved for $(length(theta_scan_results)) geometries")

log_message("All plots saved to: $PLOT_DIR")
for f in sort(readdir(PLOT_DIR))
    fp = joinpath(PLOT_DIR, f); isfile(fp) && log_message("  $f ($(format_bytes_filesize(fp)))")
end

# Stage 4: Summary
log_message("="^60)
log_message("STAGE 4: Summary Export")
log_message("="^60)

summary_path = joinpath(OUTPUT_DIR, "pipeline_summary.txt")
open(summary_path, "w") do io
    println(io, "="^60)
    println(io, "  CYAxiverse Analysis Pipeline Summary")
    println(io, "  Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "="^60); println(io)
    println(io, "Dataset Overview:")
    println(io, "  Total axions: $(length(axions))")
    println(io, "  Total geometries: $(length(geometries))")
    println(io, "  Failed geometries: $(length(axion_data.failed))")
    h11s = unique([d.h11 for d in axions])
    println(io, "  h11 range: $(minimum(h11s)) - $(maximum(h11s))")
    println(io)
    println(io, "Output Files:")
    println(io, "  Datasets:")
    for f in sort(readdir(OUTPUT_DIR))
        fp = joinpath(OUTPUT_DIR, f); isfile(fp) && println(io, "    $f ($(format_bytes_filesize(fp)))")
    end
    println(io, "  Plots:")
    for f in sort(readdir(PLOT_DIR))
        fp = joinpath(PLOT_DIR, f); isfile(fp) && println(io, "    $f ($(format_bytes_filesize(fp)))")
    end
end
log_message("Summary exported to: $summary_path")
log_message("="^60)
log_message("Pipeline complete!")
