
##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- alternate_forcings.jl
##
## Diversity effects of variation in non-precipitation climate
## forcings -- temperature (mortality and VPD) and atmospheric
## carbon concentration.
##
## by: Jacob Levine - jacoblevine@princeton.edu
## November 2022
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. define includes and parameters
##---------------------------------------------------------------

using Plots, QuadGK, DataFrames, Distributions,
    SpecialFunctions, NLsolve, AutoPreallocation,
    Profile, PProf, PlotThemes, CSV, GeoArrays, Random, Debugger,
    ForwardDiff, Polylogarithms, ProgressBars, Roots

## simulation parameters

## ecological parameters
μ::Float64 = 0.31 ## mortality rate
E::Float64 = 0.5 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 3.0  ## biomass allometric constant
F::Float64 = 200.0  ## fecundity per unit biomass
W₀::Float64 = 0.4 ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity

## include function headers
include("utility_functions.jl")
include("simulation_functions.jl")
include("eq_functions.jl")
include("meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)


## View relationship between temp and VPD
x = collect(range(100.0, 1500.0, 100))
plot(x, calc_aₘ.(x, calc_vpd(20.0, 30.0)))
plot!(x, calc_aₘ.(x, calc_vpd(20.0, 30.0), 1.04545, 1000.0, 50.0, 10.0))

##---------------------------------------------------------------
## carbon and MAP
##---------------------------------------------------------------

meq = multi_eq_carbon_map(50, 30, 0.2, 100.0, 750.0, 100, 0.2, 0.4, 5, 10.0, F, μ, 2, 0.1)
smeq = summarize_multi_eq_carbon_map(meq)
smeq.cₐ .= round.(smeq.cₐ, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.cₐ, sub.mean, group = sub.map ./ 10.0, line_z = sub.map ./ 10.0,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.map ./ 10.0, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Atmospheric carbon concentration (ppm)", ylab = "Species richness")


smeq_b = summarize_multi_eq_biomass_carbon_map(meq)
smeq_b.cₐ .= round.(smeq_b.cₐ, digits = 2)
bp = plot(smeq_b.cₐ, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0,
          fill_z = smeq_b.map ./ 10.0,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Atmospheric carbon concentration (ppm)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_carbon.pdf")

summarize_multi_eq_transpiration(meq)

plot(sub.mean, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0, fill_z = smeq_b.map ./ 10.0,
     linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Species Richness", ylab = "Total biomass")
savefig("figures/div_es_total.pdf")



##---------------------------------------------------------------
## carbon and storm freq
##---------------------------------------------------------------

calc_aₘ(500, 1000)

x = collect(range, )
plot(x, calc_aₘ.(500, calc_vpd.(x, 30.0)))

x = collect(range(100.0, 1500.0, 100))
plot(x, calc_aₘ.(x, calc_vpd(20.0, 30.0)))
plot!(x, calc_aₘ.(x, calc_vpd(20.0, 30.0), 1.04545, 1000.0, 50.0, 10.0))

meq = multi_eq_carbon_P(50, 30, 0.2, 100.0, 750.0, 100, 7.0, 15.0, 5, 0.3, F, μ, 2, 0.1)
meq[1]
smeq = summarize_multi_eq_carbon_P(meq)
smeq.cₐ .= round.(smeq.cₐ, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.cₐ, sub.mean, group = sub.P, line_z = sub.P,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.P, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Atmospheric carbon concentration (ppm)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_carbon_P(meq)
smeq_b.cₐ .= round.(smeq_b.cₐ, digits = 2)
bp = plot(smeq_b.cₐ, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P,
          fill_z = smeq_b.P,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Atmospheric carbon concentration (ppm)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/freq_carbon.pdf")

summarize_multi_eq_transpiration(meq)

plot(sub.mean, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P, fill_z = smeq_b.P,
     linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Species Richness", ylab = "Total biomass")
savefig("figures/div_es_total.pdf")


##---------------------------------------------------------------
## VPD
##---------------------------------------------------------------

meq = multi_eq_vpd_map(50, 5, 0.2, 500.0, 10000.0, 100, 0.2, 0.4, 5, 10.0, F, μ, 2, 0.1)
smeq = summarize_multi_eq_vpd_map(meq)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.vpd, sub.mean, group = sub.map ./ 10.0, line_z = sub.map ./ 10.0,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.map ./ 10.0, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_map(meq)
bp = plot(smeq_b.vpd, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0,
          fill_z = smeq_b.map ./ 10.0,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_temp.pdf")

summarize_multi_eq_transpiration(meq)

plot(sub.mean, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0, fill_z = smeq_b.map ./ 10.0,
     linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Species Richness", ylab = "Total biomass")
savefig("figures/div_es_total.pdf")


##---------------------------------------------------------------
## VPD and storm freq
##---------------------------------------------------------------

meq = multi_eq_vpd_P(50, 30, 0.2, 500.0, 10000.0, 100, 7.0, 15.0, 5, 0.3, F, μ, 2, 0.1)
smeq = summarize_multi_eq_temp_P(meq)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.vpd, sub.mean, group = sub.P, line_z = sub.P,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.P, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_P(meq)
bp = plot(smeq_b.vpd, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P,
          fill_z = smeq_b.P,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_freq.pdf")

##---------------------------------------------------------------
## Mortality (season length (temperature))
##---------------------------------------------------------------

meq = multi_eq_μ_map(50, 5, 0.2, 0.05, 0.2, 100, 0.2, 0.4, 5, 10.0, F, μ, 2, 0.1)
smeq = summarize_multi_eq_vpd_map(meq)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.vpd, sub.mean, group = sub.map ./ 10.0, line_z = sub.map ./ 10.0,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.map ./ 10.0, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_map(meq)
bp = plot(smeq_b.vpd, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0,
          fill_z = smeq_b.map ./ 10.0,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_temp.pdf")

##---------------------------------------------------------------
## VPD and storm freq
##---------------------------------------------------------------

meq = multi_eq_μ_P(50, 30, 0.2, 0.05, 0.2, 100, 7.0, 15.0, 5, 0.3, F, μ, 2, 0.1)
smeq = summarize_multi_eq_temp_P(meq)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.vpd, sub.mean, group = sub.P, line_z = sub.P,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.P, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_P(meq)
bp = plot(smeq_b.vpd, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P,
          fill_z = smeq_b.P,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_freq.pdf")
