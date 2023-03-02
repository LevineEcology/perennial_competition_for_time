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
    SpecialFunctions, LazySets, NLsolve, AutoPreallocation,
    Profile, PProf, PlotThemes, CSV, GeoArrays

## simulation parameters
Nspp::Int64 = 10; Nyr::Int64 = 1000; Ninit::Float64 = 1.0;

## ecological parameters
μ::Float64 = 0.1 ## mortality rate
E::Float64 = 0.02 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 2.5  ## biomass allometric constant
F::Float64 = 100  ## fecundity per unit biomass
T::Float64 = 40.0 ## length of interrain interval (default)
W₀::Float64 = 0.6 ## initial water content (default)

## include function headers
include("simulation_functions.jl")
include("eq_functions.jl")
include("meta_functions.jl")

## set plotting variables
theme(:dark)
my_cgrad = cgrad(:acton)

##---------------------------------------------------------------
## 1. check agreement between simulator and equilibrium calculations
##---------------------------------------------------------------

## to run in parallel must specify JULIA_NUM_THREADS environment variable
## should take about ~1min after compilation if running on four threads
out = check_eq_agreement(80, 10)
out = check_eq_agreement(40, 30)
plot_eq_agreement(out, true, "figures/sim_vs_an.pdf")

##---------------------------------------------------------------
## 2. calculate equilibrium densities under constant, initial
##    water contents and interrain interval lengths
##---------------------------------------------------------------

## First for dependence on T with 4 values of W₀
@time results = multi_eq_constant_water(30, 15, 0.4, 0.1, 0.4, 4, 1.0, 140.0, 14);
summary = summarize_multi_eq(results)

## plot
plot_multi_eq(summary, "n", :T, results[4], true, "figures/eq_nfeas_T.pdf")
plot_multi_eq(summary, "min", :T, results[4],true, "figures/eq_minfeas_T.pdf")
plot_multi_eq(summary, "avg", :T, results[4], true, "figures/eq_avgfeas_T.pdf")

## Then for dependence on W₀ with 4 values of T
@time results = multi_eq_constant_water(30, 15, 0.4, 0.1, 0.6, 14, 1.0, 140.0, 5);
summary = summarize_multi_eq(results)

plot_multi_eq(summary, "n", :W₀, results[4], true, "figures/eq_nfeas_W0.pdf")
plot_multi_eq(summary, "min", :W₀, results[4], true, "figures/eq_minfeas_W0.pdf")
plot_multi_eq(summary, "avg", :W₀, results[4], true, "figures/eq_avgfeas_W0.pdf")

##---------------------------------------------------------------
## 3. simulate population dynamics under variable W₀ and T
##---------------------------------------------------------------
spp_data = generate_spp_data(10, 0.7, 40.0, 100.0, 0.1, 2.5, 0.4,
                             0.0)

plot(spp_data.τ, spp_data.Wᵢ, seriestype = :scatter)

out1 = sim_water_only(spp_data, 3000, nrow(spp_data), 1.0)
out2 = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, true, false,
                                0.6, 0.1)

plot_simulation_dynamics(out1)
savefig("figures/simulation_novar.pdf")
plot_simulation_dynamics(out2)
savefig("figures/simulation_var.pdf")

out3 = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.0, 110.0, 30.0)

plot_simulation_dynamics(out3)

out = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.1, 40.0, 20.0)

plot_simulation_dynamics(out)
savefig("figures/simulation_varT.pdf")


@time out = multi_eq_variable_water(30, 20, 4000,
                                    0.1, 0.6, 5,
                                    0.0, 0.2, 5,
                                    1.0, 40.0, 5,
                                    0.0, 30.0, 5,
                                    0.4)

## simulations done on cluster for speed. load results CSV from file
variable_results = CSV.read("simulator_runs.csv", DataFrame)

plot_multi_eq_variable(variable_results, "n", :W₀, 30, true, "figures/eq_nfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "n", :T, 30, true, "figures/eq_nfeas_T_var.pdf")

plot_multi_eq_variable(variable_results, "min", :W₀, 30, true, "figures/eq_minfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "min", :T, 30, true, "figures/eq_minfeas_T_var.pdf")

plot_multi_eq_variable(variable_results, "avg", :W₀, 30, true, "figures/eq_avgfeas_W0_var.pdf")
plot_multi_eq_variable(variable_results, "avg", :T, 30, true, "figures/eq_avgfeas_T_var.pdf")

data = variable_results; yvar = "avg"; xvar = "T"


##---------------------------------------------------------------
## 4. Spatial patterns
##---------------------------------------------------------------

function predict_soil_moisture(elevation)
    0.9 ./ (exp.(0.0008 .* ((elevation .- minimum(elevation) ./ (maximum(elevation) - minimum(elevation))))))
end

function degrade_by_aspect(sm, aspect)
    sm / (exp(0.8 * aspect / 180))
end

import GeoArrays
elevation_grid = GeoArrays.read("elevation/dem.tif")
aspect_grid = GeoArrays.read("elevation/aspect.tif")
elevation_grid = elevation_grid[600:1800,2500:3100,:]
aspect_grid = aspect_grid[600:1800,2500:3100,:]
soil_moisture_grid = copy(elevation_grid)
soil_moisture_grid[:,:,1] = predict_soil_moisture(soil_moisture_grid[:,:,1])
soil_moisture_grid[:,:,1] = degrade_by_aspect.(soil_moisture_grid[:,:,1], aspect_grid[:,:,1])

elp = plot(elevation_grid, c = :thermal)
savefig(elp, "figures/elevation_grid.pdf")
plot(aspect_grid, c = :thermal)

plot(soil_moisture_grid)

transect = elevation_grid[:, 1500, :]

sm_sub = soil_moisture_grid[500:700,200:300,:]
plot(sm_sub)

diversity_grid = @time multi_eq_geography(soil_moisture_grid, 40)
dv_grid = copy(soil_moisture_grid)
dv_grid[:,:,1] = diversity_grid
plot(dv_grid, c = :greens)

GeoArrays.write("diversity_grid.tif", dv_grid)
dv_grid = GeoArrays.read("diversity_grid.tif")
dvp = plot(dv_grid, c = :greens)
savefig(dvp, "figures/diversity_grid.pdf")
## takes forever
@time diversity_grid2 = multi_eq_geography(sm_sub, 40)

bm_grid = copy(sm_sub)
bm_grid[:,:,1] = diversity_grid2[2]
plot(bm_grid, c = :greens)

bm_grid = copy(sm_sub)
bm_grid[:,:,1] = diversity_grid2[2]
plot(bm_grid, c = :greens)

ph_grid = copy(sm_sub)
ph_grid[:,:,1] = diversity_grid2[3]
plot(ph_grid, c = :blues)



##---------------------------------------------------------------
## 5. Simulations for fourth year talk
##---------------------------------------------------------------

spp_data = generate_spp_data(20)
