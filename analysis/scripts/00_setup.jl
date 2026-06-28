using Pkg
Pkg.activate("/opt/CYAxiverse.jl")

using CYAxiverse
using LinearAlgebra
using Statistics

using CairoMakie
using Serialization
using Printf
using Dates

ENV["MOSEKLM_LICENSE_FILE"] = "/home/cytools/mosek.lic"
ENV["newARGS"] = "docker"

PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
OUTPUT_DIR = joinpath(PROJECT_ROOT, "outputs")
PLOT_DIR = joinpath(PROJECT_ROOT, "plots")
LOG_DIR = joinpath(PROJECT_ROOT, "logs")
LOG_FILE = joinpath(LOG_DIR, "run_log.txt")
for d in [OUTPUT_DIR, PLOT_DIR, LOG_DIR]
    isdir(d) || mkpath(d)
end

include(joinpath(@__DIR__, "..", "src", "dataset_builders.jl"))
include(joinpath(@__DIR__, "..", "src", "abundance_tools.jl"))
include(joinpath(@__DIR__, "..", "src", "interaction_tools.jl"))
include(joinpath(@__DIR__, "..", "src", "statistics_tools.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting_tools.jl"))

log_message("="^60)
log_message("CYAxiverse Analysis Pipeline Started")
log_message("Host: $(gethostname())")
log_message("Julia: $(VERSION)")
log_message("CYAxiverse loaded successfully")
log_message("="^60)
