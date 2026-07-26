function light_axion_mask(ms; log10m_max=-27.0)
    return [m < log10m_max for m in ms]
end

"""
    dm_relevant_axion_mask(log10m; log10m_min=-33.0, log10m_max=2.0)

Physical mass window for CDM misalignment calculation.

Lower bound:  H₀ ≈ 10⁻³³ eV — axions below this haven't started oscillating.
Upper bound:  10² eV — heavy axions excluded by standard criterion.

Returns a boolean mask: true for axions in [log10m_min, log10m_max].
"""
function dm_relevant_axion_mask(log10m;
                                log10m_min::Real=-33.0,
                                log10m_max::Real=2.0)
    return (log10m .>= log10m_min) .& (log10m .<= log10m_max)
end

"""
    classify_axion(log10m; frozen_max=-33.0, heavy_min=2.0)

Classify an axion mass into three groups:
  :frozen   — m < H₀  (not yet oscillating, like dark energy)
  :dm       — H₀ < m < 10² eV  (CDM candidate)
  :heavy    — m > 10² eV  (excluded)
"""
function classify_axion(log10m; frozen_max::Real=-33.0, heavy_min::Real=2.0)
    if log10m < frozen_max
        return :frozen
    elseif log10m > heavy_min
        return :heavy
    else
        return :dm
    end
end

function log10_axion_density_ratio(log10m, log10f; θ=nothing)
    log10_θ = θ === nothing ? 0.0 : 2.0 * log10(θ / (π / 2))
    return log10(11.3) + log10_θ + 0.5*log10m + 2.0*log10f - 23.5
end

function log10sumexp10(vals)
    isempty(vals) && return -Inf
    m = maximum(vals)
    return m + log10(sum(10.0 .^ (vals .- m)))
end

function log10_axion_density_ratios(ms, fs; θ=nothing)
    [log10_axion_density_ratio(m, f; θ=θ) for (m, f) in zip(ms, fs)]
end

function log10_geometry_density_ratio(ms, fs; θ=nothing)
    log10sumexp10(log10_axion_density_ratios(ms, fs; θ=θ))
end

function theta_scan_geometry_density(ms, fs; thetas=10.0 .^ range(-2, 0, length=20))
    [(θ=θ, log10_ratio=log10_geometry_density_ratio(ms, fs; θ=θ)) for θ in thetas]
end
