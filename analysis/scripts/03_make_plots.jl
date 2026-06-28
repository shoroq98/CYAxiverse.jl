include(joinpath(@__DIR__, "00_setup.jl"))

log_message("Stage 3: Plot Generation")

cache_label = "--full" in ARGS ? "" : "_test"
axion_data = load_dataset("axion_dataset$(cache_label).jls")
geometry_data = load_dataset("geometry_dataset$(cache_label).jls")
axions = axion_data.axions
geometries = geometry_data.geometries

log_message("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

CairoMakie.activate!(type="png")

log_message("Generating axion-level plots...")

save_histogram(axions, :log10m, "hist_mass.png"; bins=80,
               title="Axion Mass Distribution",
               xlabel="log10(mass / eV)")

save_histogram(axions, :log10f, "hist_decay.png"; bins=80,
               title="Axion Decay Constant Distribution",
               xlabel="log10(decay constant / GeV)")

save_histogram(axions, :log10_abundance, "hist_abundance.png"; bins=80,
               title="Axion Abundance Distribution",
               xlabel="log10(Ω_a h²)")

save_scatter(axions, :log10m, :log10_abundance, "scatter_mass_abundance.png";
             title="Mass vs Abundance",
             xlabel="log10(mass / eV)", ylabel="log10(Ω_a h²)")

save_scatter(axions, :log10f, :log10_abundance, "scatter_decay_abundance.png";
             title="Decay Constant vs Abundance",
             xlabel="log10(decay constant / GeV)", ylabel="log10(Ω_a h²)")

save_scatter_color(axions, :log10m, :log10f, :log10_abundance,
                   "scatter_mass_decay_color.png";
                   title="Mass vs Decay Constant (colored by abundance)",
                   xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)",
                   clabel="log10(Ω_a h²)")

save_heatmap(axions, :log10m, :log10f, "heatmap_mass_decay.png";
             nbins_x=30, nbins_y=30,
             title="Mass vs Decay Constant Density",
             xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)")

save_scatter(axions, :log10m, :log10f, "scatter_mass_decay.png";
             title="Mass vs Decay Constant",
             xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)")

log_message("Generating geometry-level plots...")

save_histogram(geometries, :n_light_axions, "hist_light_axions.png"; bins=50,
               title="Light Axions per Geometry",
               xlabel="Number of light axions (m < 10^-27 eV)")

n_geom = length(geometries)
geom_abundances = Float64[]
geom_h11 = Int[]
geom_n_light = Int[]
for g in geometries
    matching = [d for d in axions
                if d.h11 == g.h11 && d.polytope == g.polytope && d.frst == g.frst]
    total_abund = isempty(matching) ? -Inf : log10(sum(10 .^ [a.log10_abundance for a in matching]))
    push!(geom_abundances, total_abund)
    push!(geom_h11, g.h11)
    push!(geom_n_light, g.n_light_axions)
end

geom_abundances_plot = [Float64[]]
for a in geom_abundances
    if isfinite(a)
        push!(geom_abundances_plot[1], a)
    end
end
f_geom_abund = Figure(size=(800, 600))
ax_geom_abund = Axis(f_geom_abund[1, 1], title="Total Geometry Abundance Distribution",
                     xlabel="log10(Total Ω_a h²)", ylabel="Count")
hist!(ax_geom_abund, geom_abundances_plot[1], bins=50, color=(:steelblue, 0.7))
save(joinpath(PLOT_DIR, "hist_geometry_abundance.png"), f_geom_abund)
log_message("Saved plot: hist_geometry_abundance.png")

finite_mask = isfinite.(geom_abundances)
sampled_geom_h11 = geom_h11[finite_mask]
sampled_geom_abund = geom_abundances[finite_mask]

f_scatter_geom = Figure(size=(800, 600))
ax_scatter_geom = Axis(f_scatter_geom[1, 1], title="Geometry Abundance vs h11",
                       xlabel="h11", ylabel="log10(Total Ω_a h²)")
scatter!(ax_scatter_geom, sampled_geom_h11, sampled_geom_abund, markersize=3, color=(:steelblue, 0.5))
save(joinpath(PLOT_DIR, "scatter_geom_abundance_h11.png"), f_scatter_geom)
log_message("Saved plot: scatter_geom_abundance_h11.png")

f_abund_h11 = Figure(size=(1200, 800))
ax_abund1 = Axis(f_abund_h11[1, 1], title="Geometry Abundance vs h11",
                 xlabel="h11", ylabel="log10(Total Ω_a h²)")
scatter!(ax_abund1, sampled_geom_h11, sampled_geom_abund, markersize=3, color=(:steelblue, 0.5))
ax_abund2 = Axis(f_abund_h11[1, 2], title="Light Axions vs Geometry Abundance",
                 xlabel="n_light_axions", ylabel="log10(Total Ω_a h²)")
n_light_finite = geom_n_light[finite_mask]
scatter!(ax_abund2, n_light_finite, sampled_geom_abund, markersize=3, color=(:steelblue, 0.5))
ax_abund3 = Axis(f_abund_h11[2, 1], title="Light Axions Distribution",
                 xlabel="n_light_axions", ylabel="Count")
hist!(ax_abund3, geom_n_light, bins=50, color=(:steelblue, 0.7))
save(joinpath(PLOT_DIR, "geometry_abundance_h11.png"), f_abund_h11)
log_message("Saved plot: geometry_abundance_h11.png")

log_message("Loading PCA results and generating PCA plots...")
pca_result = load_dataset("pca_results.jls")
save_pca_plot(pca_result, axions, "pca_projection.png";
              title="PCA: Axion Parameter Space",
              color_field=:log10_abundance)
save_explained_variance(pca_result, "pca_explained_variance.png";
                        title="PCA Explained Variance Ratio")

log_message("All plots saved to: $PLOT_DIR")
for f in sort(readdir(PLOT_DIR))
    fp = joinpath(PLOT_DIR, f)
    isfile(fp) && log_message("  $f ($(format_bytes_filesize(fp)))")
end
log_message("="^60)
log_message("Plot generation complete.")
