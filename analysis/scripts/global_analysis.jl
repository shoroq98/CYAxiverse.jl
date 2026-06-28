using Pkg
Pkg.activate("/opt/CYAxiverse.jl")

using CYAxiverse
using LinearAlgebra
using Statistics
using Serialization
using Printf
using Dates
using CairoMakie

ENV["MOSEKLM_LICENSE_FILE"] = "/home/cytools/mosek.lic"
ENV["newARGS"] = "docker"

# ─────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────
data_dir = length(ARGS) >= 1 ? ARGS[1] : "/mnt/results"
output_dir = data_dir
plot_dir = joinpath(output_dir, "plots")
log_dir = joinpath(output_dir, "logs")
mkpath.(plot_dir); mkpath(log_dir)

log_file = joinpath(log_dir, "global_analysis.log")

function log(msg)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    line = "[$ts] [Global] $msg"
    open(log_file, "a") do io; println(io, line); end
    @info line
end

# ─────────────────────────────────────────────────────────────────────
# Load merged datasets
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("Global Analysis — Loading merged data")
log("="^60)

axion_data = deserialize(joinpath(data_dir, "axion_dataset.jls"))
geom_data  = deserialize(joinpath(data_dir, "geometry_dataset.jls"))
axions = axion_data.axions
geometries = geom_data.geometries
log("Loaded $(length(axions)) axions from $(length(geometries)) geometries")

# ─────────────────────────────────────────────────────────────────────
# Stage 1: Summary statistics
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("Stage 1: Summary Statistics")
log("="^60)

function extract_vector(data, field; filter_fn=isfinite)
    [getfield(d, field) for d in data if filter_fn(getfield(d, field))]
end

function summary_stats(vals)
    n = length(vals)
    n == 0 && return nothing
    (n=n, mean=mean(vals), std=std(vals),
     min=minimum(vals), max=maximum(vals),
     q25=quantile(vals, 0.25), q50=quantile(vals, 0.5), q75=quantile(vals, 0.75))
end

function print_stats(label, stats)
    stats === nothing && return
    log("  $label: n=$(stats.n)  mean=$(round(stats.mean, digits=3))  std=$(round(stats.std, digits=3))")
    log("    min=$(round(stats.min, digits=3))  q25=$(round(stats.q25, digits=3))  median=$(round(stats.q50, digits=3))  q75=$(round(stats.q75, digits=3))  max=$(round(stats.max, digits=3))")
end

axion_fields = [
    (:log10m, "log10(mass / eV)"),
    (:log10f, "log10(decay constant / GeV)"),
    (:log10_abundance, "log10(Ω_a h²)"),
]

results = Dict()
for (field, label) in axion_fields
    vals = extract_vector(axions, field)
    s = summary_stats(vals)
    results[field] = s
    print_stats(label, s)
end

for (field, label) in [(:h11, "h11"), (:cy_volume, "log10(CY Volume)")]
    vals = extract_vector(geometries, field)
    if field == :cy_volume; vals = log10.(vals); end
    s = summary_stats(vals)
    results[field] = s
    print_stats(label, s)
end

# ─────────────────────────────────────────────────────────────────────
# Stage 2: Correlation + PCA
# ─────────────────────────────────────────────────────────────────────
log("Correlation table (mass, decay, abundance)...")
m_vals = extract_vector(axions, :log10m)
f_vals = extract_vector(axions, :log10f)
a_vals = extract_vector(axions, :log10_abundance)
finite_all = isfinite.(m_vals) .& isfinite.(f_vals) .& isfinite.(a_vals)
M = [m_vals[finite_all] f_vals[finite_all] a_vals[finite_all]]
C = cor(M)
log("  ρ(m,f)=$(round(C[1,2], digits=3))  ρ(m,Ω)=$(round(C[1,3], digits=3))  ρ(f,Ω)=$(round(C[2,3], digits=3))")
serialize(joinpath(data_dir, "correlation_matrix.jls"), C)

log("PCA...")
using LinearAlgebra
M_std = (M .- mean(M, dims=1)) ./ std(M, dims=1)
SVD = svd(M_std)
explained = (SVD.S .^ 2) / sum(SVD.S .^ 2)
log("  PC1: $(round(explained[1]*100, digits=1))%  PC2: $(round(explained[2]*100, digits=1))%  PC3: $(round(explained[3]*100, digits=1))%")
pca_result = (projection=M_std * SVD.Vt, loadings=SVD.Vt, explained=explained)
serialize(joinpath(data_dir, "pca_results.jls"), pca_result)

# ─────────────────────────────────────────────────────────────────────
# Stage 3: Plots
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("Stage 3: Generating Plots")
log("="^60)

CairoMakie.activate!(type="png")

function save_hist(vals, filename; title="", xlabel="", bins=80)
    finite_vals = filter(isfinite, vals)
    isempty(finite_vals) && return
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel="Count")
    hist!(ax, finite_vals, bins=bins, color=(:steelblue, 0.7))
    save(joinpath(plot_dir, filename), fig)
    log("  Saved: $filename")
end

function save_scatter(x, y, filename; title="", xlabel="", ylabel="")
    finite = isfinite.(x) .& isfinite.(y)
    count(finite) < 2 && return
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel=ylabel)
    scatter!(ax, x[finite], y[finite], markersize=3, color=(:steelblue, 0.5))
    save(joinpath(plot_dir, filename), fig)
    log("  Saved: $filename")
end

save_hist(m_vals, "hist_mass.png", title="Axion Mass Distribution", xlabel="log10(mass / eV)")
save_hist(f_vals, "hist_decay.png", title="Axion Decay Constant Distribution", xlabel="log10(decay constant / GeV)")
save_hist(a_vals, "hist_abundance.png", title="Axion Abundance Distribution", xlabel="log10(Ω_a h²)")
save_scatter(m_vals, a_vals, "scatter_mass_abundance.png",
    title="Mass vs Abundance", xlabel="log10(mass / eV)", ylabel="log10(Ω_a h²)")
save_scatter(f_vals, a_vals, "scatter_decay_abundance.png",
    title="Decay Constant vs Abundance", xlabel="log10(decay constant / GeV)", ylabel="log10(Ω_a h²)")
save_scatter(m_vals, f_vals, "scatter_mass_decay.png",
    title="Mass vs Decay Constant", xlabel="log10(mass / eV)", ylabel="log10(decay constant / GeV)")

# Light axions per geometry
n_light = [g.n_light_axions for g in geometries]
save_hist(n_light, "hist_light_axions.png",
    title="Light Axions per Geometry", xlabel="Number of light axions (m < 10^-27 eV)")

# Geometry abundance
DM_OBS = 0.12
geom_abundances = Float64[]
for g in geometries
    matching = [a for a in axions
        if a.h11 == g.h11 && a.polytope == g.polytope && a.frst == g.frst]
    total = isempty(matching) ? -Inf : log10(sum(10.0 .^ [a.log10_abundance for a in matching]) / DM_OBS)
    push!(geom_abundances, total)
end
finite_ga = isfinite.(geom_abundances)
h11_vals = [g.h11 for g in geometries]
save_hist(geom_abundances[finite_ga], "hist_geometry_abundance.png",
    title="Total Geometry Abundance", xlabel="log10(Total Ω_a h²)")

fig_scatter = Figure(size=(800, 600))
ax_scatter = Axis(fig_scatter[1,1],
    title="Geometry DM Density Ratio vs h11",
    xlabel="h11", ylabel="log10(Ωₐ h² / 0.12)")
scatter!(ax_scatter, h11_vals[finite_ga], geom_abundances[finite_ga],
    markersize=3, color=(:steelblue, 0.5))
hlines!(ax_scatter, [0], color=:black, linestyle=:dash, linewidth=2)
save(joinpath(plot_dir, "scatter_geom_DMratio_h11.png"), fig_scatter)
log("  Saved: scatter_geom_DMratio_h11.png")

# ─────────────────────────────────────────────────────────────────────
# Stage 4: Summary export
# ─────────────────────────────────────────────────────────────────────
log("="^60)
log("Stage 4: Summary Export")
log("="^60)

summary_path = joinpath(output_dir, "pipeline_summary.txt")
open(summary_path, "w") do io
    println(io, "="^60)
    println(io, "  CYAxiverse Pipeline Summary")
    println(io, "  Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "="^60); println(io)
    println(io, "Dataset Overview:")
    println(io, "  Total axions: $(length(axions))")
    println(io, "  Total geometries: $(length(geometries))")
    println(io, "  Failed geometries: $(length(axion_data.failed))")
    h11s = sort(unique([d.h11 for d in axions]))
    println(io, "  h11 range: $(minimum(h11s)) - $(maximum(h11s))")
    println(io)
    println(io, "Mass (log10 eV):")
    if results[:log10m] !== nothing
        s = results[:log10m]
        println(io, "  Mean: $(round(s.mean, digits=2))  Std: $(round(s.std, digits=2))")
        println(io, "  Range: [$(round(s.min, digits=2)), $(round(s.max, digits=2))]")
    end
    println(io)
    println(io, "Decay Constant (log10 GeV):")
    if results[:log10f] !== nothing
        s = results[:log10f]
        println(io, "  Mean: $(round(s.mean, digits=2))  Std: $(round(s.std, digits=2))")
        println(io, "  Range: [$(round(s.min, digits=2)), $(round(s.max, digits=2))]")
    end
    println(io)
    println(io, "PCA:")
    println(io, "  PC1: $(round(explained[1]*100, digits=1))%")
    println(io, "  PC2: $(round(explained[2]*100, digits=1))%")
    println(io, "  PC3: $(round(explained[3]*100, digits=1))%")
end

log("Summary exported to: $summary_path")
log("="^60)
log("Global analysis complete!")
log("="^60)
