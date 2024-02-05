
##---------------------------------------------------------------
## WATER_ONLY_SIMULATOR -- utility_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## December 2023
##
## This script contains general utility functions, called by many different
## functions for various purposes.
##---------------------------------------------------------------

using DataFrames: DataFrameColumns
## define parameters


"""
    generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, n_ht::Int64 = 1, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.0,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.001,
                           spread_multiplier::Float64 = 25.0, min_ht::Float64 = 0.3, max_ht::Float64 = 0.7,
                           mechanistic = false, Xmax::Float64 = 8.0,
                           aₘ::Float64 = 0.03, γ::Float64 = 0.5,
                           rᵣ::Float64 = 0.01, cₓ::Float64 = 1.5)

Generates a community of species for use in simulations or equilibrium calculations. The size of the
community is specified through the `Nspp` parameter. The physical characteristics
of the species are randomly assigned. Specifically, the non-water-limited growth rate of each speices is
drawn from a uniform distribution. The break-even times of the species are then calculated and then their
critical water contents are calculated to adhere to a coexistence-maintaining tradeoff. The user has the
option to add scatter around this tradeoff by setting `tradeoff_sd` to a value greater than 0. The height
strategy of each species is also drawn from a uniform distribution. The number of height strategies is specified
by the user in `n_ht.`

"""
function generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, n_ht::Int64 = 1, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.0,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.001,
                           spread_multiplier::Float64 = 25.0, min_ht::Float64 = 0.3, max_ht::Float64 = 0.7,
                           mechanistic = false, Xmax::Float64 = 8.0,
                           aₘ::Float64 = 0.03, γ::Float64 = 0.5,
                           rᵣ::Float64 = 0.01, cₓ::Float64 = 1.5)

    if mechanistic

        ## generate a list of species with X (wood density) drawn from a uniform distribution
        X = rand(Uniform(0.1, Xmax), Nspp)
        spp_data = DataFrame(X = X,
                             aₘ = repeat([aₘ], Nspp),
                             C₁ = (aₘ - (γ * rᵣ)) ./ (cₓ .* X),
                             C₂ = repeat([C₂max], Nspp))

    else

        ## generate a list of species with C₁ drawn from a uniform distribution
        spp_data = DataFrame(C₁ = rand(Uniform(0.0000005, C₁max), Nspp) .* 365.0,
                             C₂ = repeat([C₂max], Nspp) .* 365.0);

    end

    ## calculate the break-even time of each species
    spp_data.τ = broadcast(calc_τ, spp_data.C₁, spp_data.C₂, F, μ, T, b);
    ## calculate the critical water content for each species from an exponential tradeoff
    spp_data.Wᵢ = Wmax .* exp.(-tradeoff_exp .* (5 .- spread_multiplier * spp_data.C₁)) .+
        rand(Normal(0, tradeoff_sd), Nspp) ## add random noise, if specified by user

    ## sort species in order of drought tolerance from least to most
    spp_data = sort(spp_data, :Wᵢ, rev = true)
    ## ensure characteristics meet threshold requirements
    spp_data.Wᵢ = replace(x -> isless(x, 0) ? 1e-4 : x, spp_data.Wᵢ)
    spp_data.Wᵢ = replace(x -> isless(Wmax, x) ? Wmax : x, spp_data.Wᵢ)
    ## create column for spp names
    spp_data.spp = Int[1:1:Nspp;];

    ## assign height strategies at random
    hts = rand(Uniform(min_ht, max_ht), n_ht)
    spp_data.ht = rand(hts, nrow(spp_data))

    return spp_data

end;

"""
    calc_vpd(t::Float64, rh::Float64)

Calculate vapor pressure deficit (`dₛ`) from temperature in C and relative humidity
as a percent (i.e. 20 = 20 percent).

"""
function calc_vpd(t::Float64, rh::Float64)
    es = 0.6108 * exp((17.27 * t) / (t + 237.3))
    ea = (rh / 100) * es
    (es - ea) * 1000
end;


"""
    calc_aₘ(cₐ::Float64, dₛ::Float64, V::Float64 = 1.04545,
                 d₀::Float64 = 1000, Γ::Float64 = 50, m::Float64 = 5.6,
                 rₗ::Float64 = 0.1, ω::Float64 = 450)

Calculate maximum carbon assimilation rate (`aₘ`) from atmospheric carbon concentration `cₐ`,
vapor pressure deficit `dₛ`, and physiological parameters.

"""
function calc_aₘ(cₐ::Float64, dₛ::Float64, V::Float64 = 1.04545,
                 d₀::Float64 = 1000.0, Γ::Float64 = 50.0, m::Float64 = 5.6,
                 ω::Float64 = 450.0, rₗ::Float64 = 0.001)

    ((V * (cₐ - ((cₐ - Γ) * (1 + (dₛ / d₀)) / m))) /
        (ω + (cₐ - ((cₐ - Γ) * (1 + (dₛ / d₀)) / m)))) - rₗ

end;


"""
    calc_V(aₘ::Float64, cₐ::Float64, dₛ::Float64, d₀::Float64 = 1000,
           Γ:Float64 = 50, m::Float64 = 5.6, rₗ::Float64 = 0.1, ω::Float64 = 450)

Calculate maximum photosynthetic rate (`Vmax`) from

"""
function calc_V(aₘ::Float64, cₐ::Float64, dₛ::Float64, d₀::Float64 = 1000.0,
                Γ::Float64 = 50.0, m::Float64 = 5.6, ω::Float64 = 450.0, rₗ::Float64 = 0.001)

    (aₘ + rₗ) * (ω + (cₐ - ((cₐ - Γ) * (1 + (dₛ / d₀)) / m))) / (cₐ - ((cₐ - Γ) * (1 + (dₛ / d₀)) / m))

end;


"""
    adjust_spp_data(spp_data::DataFrame, t::Float64, rh::Float64, cₐ::Float64,
                         default_t::Float64 = 25, default_rh::Float64 = 0.3,
                         default_cₐ::Float64 = 419, b::Float64 = 3)

Adjust C₁ values in an spp_data object for changes in vpd (temperature and relative humidity) and
atmospheric carbon concentration.

"""
function adjust_spp_data!(spp_data::DataFrame, t::Float64, rh::Float64, cₐ::Float64,
                          default_t::Float64 = 24.0, default_rh::Float64 = 30.0,
                          default_cₐ::Float64 = 280.0, b::Float64 = 3.0,
                          γ::Float64 = 0.5, rᵣ::Float64 = 0.01, cₓ::Float64 = 1.5,
                          vpd = missing)

    V = calc_V.(spp_data.aₘ, default_cₐ, calc_vpd(default_t, default_rh))
    if ismissing(vpd)
        aₘ = calc_aₘ.(cₐ, calc_vpd(t, rh), V)
    else
        aₘ = calc_aₘ.(cₐ, vpd, V)
    end

    spp_data.C₁ = (aₘ .- (γ * rᵣ)) ./ (cₓ .* spp_data.X)
    return(spp_data)

end;



"""
    sigma(μ::Float64, ex::Float64)

Calculates polylogarithm term that re-occurs in many expressions.
"""
function sigma(μ::Float64, ex::Float64)
    Float64(polylog(-ex, exp(-μ)))
end;


"""
    calc_g(t::Float64, C₁::Float64, C₂::Float64, T::Float64)

Calculates the average growth rate of a species over an inter-storm interval.
"""
function calc_g(t::Float64, C₁::Float64, C₂::Float64, T::Float64)
    (t * C₁ - (T-t) * C₂) / T
end;

"""
    calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)

Calculates the equilibrium leaf area given a vector of equilibrium abundances `eqN.`
"""
function calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)
    eqN ./ (μ .* F)
end


"""
    calc_τ(C₁::Float64, C₂::Float64, F::Float64, μ::Float64, T::Float64, b::Float64 = 2.5)

Calculates break-even time for species when strictly water limited (open canopy system).
"""
function calc_τ(C₁::Float64, C₂::Float64, F::Float64, μ::Float64, T::Float64, b::Float64 = 2.5)
    s = sigma(μ*T, b-1)
    (T/(C₁ + C₂)) * ((1/(F * T^b * s))^(1/(b-1)) + C₂)
end

##---------------------------------------------------------------
## Plot generating functions
##---------------------------------------------------------------

"""
    plot_simulation_dynamics(results, save::Bool = false, filename = "")

Generates a plot of population density through time from simulation results as
generated by `sim_water_ppa`, `sim_water_only`, or `sim_ppa.`

"""
function plot_simulation_dynamics(results, save::Bool = false, filename = "")
    pdata = stack(results[2])
    pdata.variable = parse.(Int64, pdata.variable)

    t = zeros(length(results[13]))
    t[1] = results[13][1]
    for i in 2:length(results[13])
        t[i] = t[i-1] + results[13][i]
    end
    t = repeat(t, outer = length(unique(pdata.variable)))

    p = plot(t, pdata.value, group = pdata.variable, line_z = pdata.variable,
             ylim = [0, round(maximum(pdata.value))+1.0], xlim = [0, maximum(t)],
             seriescolor = my_cgrad, seriestype = :line,
             legend = :bottomright, frame = :box, grid = false, linewidth = 2.5,
             xlab = "Time (years)", ylab = "Population density")

    pdata.value

    if save
        savefig(p, filename)
    end

    return p
end

"""
    plot_simulation_dynamics_stochastic(results, rainfall_regime, save::Bool = false, filename = "")

Generates plot of population density over time for stochastic simulation results
"""
function plot_simulation_dynamics_stochastic(results, rainfall_regime, save::Bool = false, filename = "")

    t = Vector{Float64}(undef, length(rainfall_regime[1]))
    t[1] = rainfall_regime[2][1]
    for i in 2:length(t)
        t[i] = t[i-1] + rainfall_regime[2][i]
    end

    pdata = stack(results[2])
    pdata.variable = parse.(Int64, pdata.variable)

    p = plot(repeat(t, outer = length(unique(pdata.variable))), pdata.value, group = pdata.variable,line_z = pdata.variable,
             ylim = [0, round(maximum(pdata.value))+1.0], xlim = [0, maximum(t)],
             seriescolor = my_cgrad, seriestype = :line,
             legend = :none, frame = :box, grid = false, linewidth = 2.5)
    p2 = plot_rainfall_regime(rainfall_regime)
    p2 = plot!(t, out[5], color = :blue, ylim = [0, 0.4], xlim = [0, maximum(t)],)

    p3 = plot(p, p2, layout = (2,1))

    if save
        savefig(p, filename)
    end

    return p3
end

"""
    plot_canopy_cover(result::Any)


"""
function plot_canopy_cover(result)
    P = plot(1:nrow(result[8]), result[8][:,2], line_z = 1, linewidth = 2.5, frame = :box,
             seriescolor = my_cgrad, colorbar = :none,
             ylim = [0.0, maximum([1.0, maximum(Matrix(result[8][:,2:ncol(result[8])]))])+0.1],
             label = "species 1", legend = :topright)
    if ncol(result[8]) > 2

        for i in 3:ncol(result[8])
            p = plot!(1:nrow(result[8]), result[8][:,i], line_z = i-1, colorbar = :none,
                      seriescolor = my_cgrad, linewidth = 2.5, label = "species " * string(i-1))
        end
        rowsums = zeros(nrow(result[8]))

        for i in 1:nrow(result[8])
            for j in 2:ncol(result[8])
                rowsums[i] = rowsums[i] + result[8][i,j]
            end
        end
        p = plot!(1:nrow(result[8]), rowsums, linewidth = 2.5,
              color = "black", label = "total")
    end
    return p
end


"""
    plot_rainfall_regime(rainfall_regime)

"""
function plot_rainfall_regime(rainfall_regime, max_t = missing, ymax = missing)

    t = Vector{Float64}(undef, length(rainfall_regime[1]))
    t[1] = rainfall_regime[2][1]
    for i in 2:length(t)
        t[i] = t[i-1] + rainfall_regime[2][i]
    end

    w = copy(rainfall_regime[1])

    if !ismissing(max_t)
        w = w[t .< max_t]
        t = t[t .< max_t]
    end

    if !ismissing(ymax)
        ylims = [0.0, ymax]
    else
        ylims = [0.0, maximum(w)+(0.1*maximum(w))]
    end

    p = plot(t, w, seriestype = :scatter, color = :black, markersize = 0,
             frame = :box, legend = :none, grid = false, xlim = [0.0, maximum(t)],
             ylim = ylims,
             xlab = "Time (years)", ylab = "Storm size")

    for i in 1:length(t)
       p = plot!(vcat(t[i], t[i]), vcat(0.0, w[i]), color = :black, linewidth = 2)
    end

    return p

end
