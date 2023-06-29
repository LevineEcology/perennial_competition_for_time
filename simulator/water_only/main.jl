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
multi_1 = multi_eq(30, 20, 0.3,
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
Nyr = 200
spp_data = generate_spp_data(5, 0.6, 1.0 / P, F, μ, 2.5, 0.05, 0.0)

out = sim_water_only(spp_data, Nyr, 5, Ninit, μ, F, P, 16.0, θ_fc)
plot_simulation_dynamics(out)

rr = generate_rainfall_regime(Nyr, P, 2.5, 4.0, 0.5)
out = sim_water_only_stochastic(spp_data, nrow(spp_data), 1.0, rr,
                                0.4, μ, F)
plot_simulation_dynamics_stochastic(out, rr)


##---------------------------------------------------------------
## 2. calculate equilibrium densities under constant, initial
##    water contents and interrain interval lengths
##---------------------------------------------------------------

## First for dependence on T with 4 values of W₀
@time results = multi_eq_constant_water(20, 20, 0.4, 0.2, 0.4, 2, 0.1, 20.0, 10,
                                        0.5, 0.01);
summary = summarize_multi_eq(results)

## plot
plot_multi_eq(summary, "n", :T, results[4])
plot_multi_eq(summary, "n", :T, results[4], true, "figures/eq_nfeas_T.pdf")
plot_multi_eq(summary, "min", :T, results[4],true, "figures/eq_minfeas_T.pdf")
plot_multi_eq(summary, "avg", :T, results[4], true, "figures/eq_avgfeas_T.pdf")

## Then for dependence on W₀ with 4 values of T
@time results = multi_eq_constant_water(30, 50, 0.4, 0.1, 0.6, 18, 20.0, 140.0, 4);
summary = summarize_multi_eq(results)

plot_multi_eq(summary, "n", :W₀, results[4], true, "figures/eq_nfeas_W0.pdf")
plot_multi_eq(summary, "min", :W₀, results[4], true, "figures/eq_minfeas_W0.pdf")
plot_multi_eq(summary, "avg", :W₀, results[4], true, "figures/eq_avgfeas_W0.pdf")

##---------------------------------------------------------------
## 3. simulate population dynamics under variable W₀ and T
##---------------------------------------------------------------
spp_data = generate_spp_data(10, 0.7, 100.0, 100.0, 0.1, 2.5, 0.4,
                             0.0)

plot(spp_data.τ, spp_data.Wᵢ, seriestype = :scatter)

out1 = sim_water_only(spp_data, 3000, nrow(spp_data), 1.0, 0.1, 0.6, 100.0)
out2 = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, true, false,
                                0.6, 0.2, 100.0)

plot_simulation_dynamics(out1)
savefig("figures/simulation_novar.pdf")
plot_simulation_dynamics(out2)
savefig("figures/simulation_var.pdf")

out3 = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.0, 100.0, 30.0)

plot_simulation_dynamics(out3)

out = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.1, 40.0, 20.0)

plot_simulation_dynamics(out)
savefig("figures/simulation_varT.pdf")

## simulations done on cluster for speed. load results CSV from file
variable_results = CSV.read("simulator_runs.csv", DataFrame)

unique(variable_results.Tmean)
p_data = variable_results[variable_results.Tmean .== 110.0, :]
p_data = p_data[p_data.Tsd .== 0.0, :]
p_data = p_data[p_data.var .== "n", :]

plot(p_data.W₀mean, p_data.W₀sd, p_data.mean,
             st = :surface,
             #surfacecolor = subdata[:,groupvar],
             seriescolor = my_cgrad,
             xflip = true,
             legend = :topleft, frame = :box,  linewidth = 3, fillalpha = 0.7, colorbar = false,
              zlab = "# species persisting")

unique(variable_results.Tmean)
p_data = variable_results[variable_results.W₀mean .== 0.2, :]
p_data = p_data[p_data.W₀sd .== 0.0, :]
p_data = p_data[p_data.var .== "n", :]

plot(p_data.Tmean, p_data.Tsd, p_data.mean,
             st = :surface,
             #surfacecolor = subdata[:,groupvar],
             seriescolor = my_cgrad,
             xflip = true,
             legend = :topleft, frame = :box,  linewidth = 3, fillalpha = 0.7, colorbar = false,
              zlab = "# species persisting")

plot_multi_eq_variable(variable_results[variable_results.Tmean .== 50.0, :], "n", :W₀, 30, true, "figures/eq_nfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "n", :T, 30, true, "figures/eq_nfeas_T_var.pdf")

plot_multi_eq_variable(variable_results, "min", :W₀, 30, true, "figures/eq_minfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "min", :T, 30, true, "figures/eq_minfeas_T_var.pdf")

plot_multi_eq_variable(variable_results, "avg", :W₀, 30, true, "figures/eq_avgfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "avg", :T, 30, true, "figures/eq_avgfeas_T_var.pdf")

data = variable_results; yvar = "avg"; xvar = "T"
