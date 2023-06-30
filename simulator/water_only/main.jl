##---------------------------------------------------------------
## WATER_ONLY_SIMULATOR -- main.jl
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
    ForwardDiff, Polylogarithms, ProgressBars

## simulation parameters
Nspp::Int64 = 10; Nyr::Int64 = 1000; Ninit::Float64 = 1.0;

## ecological parameters
μ::Float64 = 0.1 ## mortality rate
E::Float64 = 0.02 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 2.5  ## biomass allometric constant
F::Float64 = 10.0  ## fecundity per unit biomass
W₀::Float64 = 0.4 ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity

## include function headers
include("simulation_functions.jl")
include("eq_functions.jl")
include("meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

##---------------------------------------------------------------
## 0. simulator checks
##---------------------------------------------------------------
μ = 0.5
F = 10.0
P = 10
Nyr = 400

spp_data = generate_spp_data(5, 0.6, 1.0 / P, F, μ, 2.5, 0.05, 0.0)
@time out = sim_water_only(spp_data, Nyr, 5, Ninit, μ, F, P, 16.0, θ_fc)

y = calc_eqN(spp_data, F, E, 0.4, θ_fc, μ, 1.0 / P, P, false)
Vector(out[2][Nyr*P,2:ncol(out[2])]) ./ y.eqN
plot(Vector(out[2][Nyr*P,2:ncol(out[2])]), y.eqN, seriestype = :scatter)
plot_simulation_dynamics(out)

##---------------------------------------------------------------
## 1. check agreement between simulator and equilibrium
##---------------------------------------------------------------

## to run in parallel must specify JULIA_NUM_THREADS environment variable
## can be slow depending on how many threads you have set up.
out = check_eq_agreement(10, 10, 3000, W₀, 10, 4.0, μ, F)
plot_eq_agreement(out, true, "figures/sim_vs_an.pdf")

##---------------------------------------------------------------
## 2. For constant MAP and storm frequencies
##---------------------------------------------------------------

μ = 0.2
multi_1 = multi_eq(30, 20, 0.4,
                   0.3*P, 0.8*P, 6, 4, 50, 30,
                   F, μ)
summary_1 = summarize_multi_eq(multi_1)
sub = summary_1[summary_1.var .== "n",:]

plot(sub.T, sub.mean, group = sub.total, line_z = sub.total,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.T), maximum(sub.T)],
     legend = :none, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")

sub = summarize_multi_eq_leafarea(multi_1)
plot(sub.T, sub.mean, group = sub.total, line_z = sub.total,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.T), maximum(sub.T)],
     legend = :none, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")

sub = summarize_multi_eq_transpiration(multi_1)
plot(sub.T, sub.mean, group = sub.total, line_z = sub.total,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.T), maximum(sub.T)],
     legend = :none, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")

##---------------------------------------------------------------
## 3. For variable MAP and storm frequencies
##---------------------------------------------------------------
P = 10
μ = 0.5
Nyr = 100

spp_data = generate_spp_data(10, 0.6, 1.0 / P, F, μ, 2.5, 0.05, 0.0)

out = sim_water_only(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 16.0, θ_fc)
plot_simulation_dynamics(out)

rr = generate_rainfall_regime(Nyr, P, 2.5, 0.5, 0.5)
out = sim_water_only_stochastic(spp_data, nrow(spp_data), 1.0, rr,
                                0.4, μ, F)
plot_simulation_dynamics_stochastic(out, rr)

rr = generate_rainfall_regime(Nyr, P, 2.5, 0.5, 0.5, true)
out = sim_water_only_stochastic(spp_data, nrow(spp_data), 1.0, rr,
                                0.4, μ, F)
plot_simulation_dynamics_stochastic(out, rr)


la = Matrix(out[1])
mean(sum(la[size(la)[1]-(0.1*size(la)[1]):size(la)[1],2:size(la)[2]], dims = 2))

transpiration = (out[6] - out[5]) ./ rr[2]
l = length(transpiration)
mean(transpiration[length(transpiration)-Int(round(0.1*length(transpiration))):length(transpiration)])

length(transpiration)
##---------------------------------------------------------------
## multi
##---------------------------------------------------------------

var_t = multi_eq_variable_total(20, 10, 400, 0.6, 0.015, 8.0, 12,
                              0.0, 4.5, 5, 10, 10.0, 10.0, 0.5, false)
CSV.write("../../data/var_t_water_only.csv", var_t)
summary_var_t = summarize_multi_eq_variable_total(var_t)

sub = summary_var_t[summary_var_t.var .== "n",:]
plot(sub.totalmean, sub.mean, group = sub.totalsd, line_z = sub.totalsd,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.totalmean), maximum(sub.totalmean)],
     ylim = [0.0, 20], frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")

var_p = multi_eq_variable_P(20, 10, 400, 0.6, 2, 40, 12,
                            0.5, 10.0, 5, 4.0, 1.5, 10.0, 0.5, false)
CSV.write("../../data/var_p_water_only.csv", var_p)
summary_var_p = summarize_multi_eq_variable_total(var_p)

sub = summary_var_p[summary_var_p.var .== "n",:]
plot(sub.totalmean, sub.mean, group = sub.totalsd, line_z = sub.totalsd,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.totalmean), maximum(sub.totalmean)],
     ylim = [0.0, 20], frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")
