
##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- main.jl
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
Nspp::Int64 = 4; Nyr::Int64 = 400; Ninit::Float64 = 1.0;

## ecological parameters
E::Float64 = 0.5 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 3.0  ## biomass allometric constant
F::Float64 = 200.0  ## fecundity per unit biomass
W₀::Float64 = 0.4 ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity
μ::Float64 = 0.31 ## mortality rate
uf::Float64 = 0.1 ## understory growth reduction factor

## include function headers
include("utility_functions.jl")
include("simulation_functions.jl")
include("eq_functions.jl")
include("meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

## generate plant community for simulations
Random.seed!(4)
spp_data = generate_spp_data(4, 0.7, 1, 1.0 / P, F, μ, b, 0.4, 0.0, 0.0001, 0.00005, 11.0)

P::Float64 = 5.0 ## storm frequency
mp::Float64 = 0.2 ## mean precipitation

## check that stochastic simulator and standard simulator return same results
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.0, true, false)
plot_rainfall_regime(r)
length(r[1])

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

plot_simulation_dynamics(out)

out2 = sim_water_ppa(spp_data, Int(length(r[1]) / P), nrow(spp_data), 1.0, μ, F, P, mp, θ_fc, zeros(1,1),
                     true, 0.4, 3.0, 0.1)
plot_simulation_dynamics(out2)


##---------------------------------------------------------------
## Variation in mean annual precip
##---------------------------------------------------------------

##---------------------------------------------------------------
## First for moderate mean map (all coexisting)
##---------------------------------------------------------------

P = 5.0
mp = 0.2

sd1 = 0.0
sd2 = 0.1
sd3 = 0.2
sd4 = 0.4

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 1000.0, mp, sd1, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(400, P, 1000.0, mp, sd2, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(400, P, 1000.0, mp, sd3, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)


r_4 = generate_rainfall_regime(400, P, 1000.0, mp, sd4, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, missing, ymax)
plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_mid_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, missing, ymax)
plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_mid_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, missing, ymax)
plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_mid_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, missing, ymax)
plot(plot(rplot_4, xlab = "", title = "σ map = " * string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_mid_4.pdf")

##---------------------------------------------------------------
## For low mean map (only latest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 5.0
mp = 0.05

sd1 = 0.0
sd2 = 0.033
sd3 = 0.1
sd4 = 0.45

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 1000.0, mp, sd1, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(400, P, 1000.0, mp, sd2, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(400, P, 1000.0, mp, sd3, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)


r_4 = generate_rainfall_regime(400, P, 1000.0, mp, sd4, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, missing, ymax)
plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_low_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, missing, ymax)
plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_low_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, missing, ymax)
plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_low_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, missing, ymax)
plot(plot(rplot_4, xlab = "", title = "σ map = " * string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_low_4.pdf")


##---------------------------------------------------------------
## For high mean map (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 5.0
mp = 0.5

sd1 = 0.0
sd2 = 0.13
sd3 = 0.25
sd4 = 0.5

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 1000.0, mp, sd1, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(400, P, 1000.0, mp, sd2, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(400, P, 1000.0, mp, sd3, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)


r_4 = generate_rainfall_regime(400, P, 1000.0, mp, sd4, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, missing, ymax)
plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_high_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, missing, ymax)
plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_high_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, missing, ymax)
plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_high_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, missing, ymax)
plot(plot(rplot_4, xlab = "", title = "σ map = " * string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_map_high_4.pdf")



##---------------------------------------------------------------
## For moderate mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 10.0
mp = 0.3

sd1 = 0.0
sd2 = 5.0
sd3 = 10.0
sd4 = 15.0

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(200, P, 1000.0, mp, 0.0, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(200, P, sd2, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(200, P, sd3, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)


r_4 = generate_rainfall_regime(200, P, sd4, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
plot(plot(rplot_1, xlab = "", title = "σ P = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_mid_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
plot(plot(rplot_2, xlab = "", title = "σ P = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_mid_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
plot(plot(rplot_3, xlab = "", title = "σ P = " *  string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_mid_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, 10, ymax)
plot(plot(rplot_4, xlab = "", title = "σ P = " *  string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_mid_4.pdf")


##---------------------------------------------------------------
## For small mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 3.0
mp = 0.3

sd1 = 0.0
sd2 = 5.0
sd3 = 10.0
sd4 = 15.0

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(200, P, 1000.0, mp, 0.0, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(200, P, sd2, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(200, P, disp3, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)


r_4 = generate_rainfall_regime(200, P, disp4, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
plot(plot(rplot_1, xlab = "", title = "σ P = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_low_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
plot(plot(rplot_2, xlab = "", title = "σ P = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_low_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
plot(plot(rplot_3, xlab = "", title = "σ P = " * string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_low_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, 10, ymax)
plot(plot(rplot_4, xlab = "", title = "σ P = " * string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_low_4.pdf")


##---------------------------------------------------------------
## For high mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 20.0
mp = 0.3

sd1 = 0.0
sd2 = 5.0
sd3 = 10.0
sd4 = 20.0

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(200, P, 0.0, mp, 0.0, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics(out)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(200, P, sd2, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics(out)


r_3 = generate_rainfall_regime(200, P, sd3, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics(out)

r_4 = generate_rainfall_regime(200, P, sd4, mp, 0.0, false, true, false)

out = sim_water_ppa_stochastic(spp_data, length(r_4[1]), nrow(spp_data), 1.0, r_4, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_4 = plot_simulation_dynamics(out)

ymax = maximum(r_4[1]) + (0.1 * maximum(r_4[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
plot(plot(rplot_1, xlab = "", title = "σ P = " * string(sd1)), plot(dynplot_1, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_high_1.pdf")

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
plot(plot(rplot_2, xlab = "", title = "σ P = " * string(sd2)), plot(dynplot_2, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_high_2.pdf")

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
plot(plot(rplot_3, xlab = "", title = "σ P = " * string(sd3)), plot(dynplot_3, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_high_3.pdf")

rplot_4 = plot_rainfall_regime(r_4, 10, ymax)
plot(plot(rplot_4, xlab = "", title = "σ P = " * string(sd4)), plot(dynplot_4, colorbar = :none), layout = [1,1])
savefig("figures/stochastic_P_high_4.pdf")


##---------------------------------------------------------------
## multi_eq methods
##
##---------------------------------------------------------------


out = multi_eq_variable_map(15, 5, 300, 0.4, 0.05, 0.5, 20, 0.0, 0.5, 4,
                            5.0, 100.0, 2, F, μ, false, b)

smeq = summarize_multi_eq_variable_map(out)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.mapmean, sub.mean, group = sub.mapsd, line_z = sub.mapsd, ribbon = sub.sd, fill_z = sub.mapsd, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Mean annual precip", ylab = "Species richness")


out2 = multi_eq_variable_P(15, 5, 300, 0.4, 3.0, 20.0, 20, -3.0, 3.0, 3,
                            0.3, 2, F, μ, false, b)

smeq = summarize_multi_eq_variable_P(out2)
sub = smeq[smeq.var .== "n", :]
sub.Pdisp .= round.(sub.Pdisp, digits = 2)
np = plot(sub.Pmean, sub.mean, group = sub.Pdisp, line_z = log.(sub.Pdisp), ribbon = sub.sd,
          fill_z = log.(sub.Pdisp), linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Storm frequency (storms per month)", ylab = "Species richness")
savefig("figures/P_stochastic.pdf")
