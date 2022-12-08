using Profile: profile_printing_listener
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
    Profile, PProf, PlotThemes, CSV

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
@time results = multi_eq_constant_water(30,8, 0.1, 0.6, 6, 1.0, 100.0, 6);
summary = summarize_multi_eq(results)

theme(:dark)
my_cgrad = cgrad(:acton)

plot_multi_eq(summary, "n", :W₀, results[4], true, "figures/eq_nfeas_W0.pdf")
plot_multi_eq(summary, "n", :T, results[4], true, "figures/eq_nfeas_T.pdf")

plot_multi_eq(summary, "min", :W₀, results[4], true, "figures/eq_minfeas_W0.pdf")
plot_multi_eq(summary, "min", :T, results[4],true, "figures/eq_minfeas_T.pdf")

plot_multi_eq(summary, "avg", :W₀, results[4], true, "figures/eq_avgfeas_W0.pdf")
plot_multi_eq(summary, "avg", :T, results[4], true, "figures/eq_avgfeas_T.pdf")

##---------------------------------------------------------------
## 3. simulate population dynamics under variable W₀ and T
##---------------------------------------------------------------
spp_data = generate_spp_data(10, 0.8, 40.0, 100.0, 0.1, 2.5, 0.4,
                             0.0)

out = sim_water_only(spp_data, 3000, nrow(spp_data), 1.0)
out = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, true, false,
                                0.6, 0.1)

plot_simulation_dynamics(out)

out = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.1)

plot_simulation_dynamics(out)

out = sim_water_only_stochastic(spp_data, 3000, nrow(spp_data), 1.0, false, true,
                                0.6, 0.1, 40.0, 20.0)

plot_simulation_dynamics(out)

CSV.read()


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

plot(elevation_grid, c = :thermal)
plot(soil_moisture_grid)

transect = elevation_grid[:, 1500, :]

sm_sub = soil_moisture_grid[500:700,200:300,:]
plot(sm_sub)

diversity_grid = @time multi_eq_geography(soil_moisture_grid, 40)
dv_grid = copy(soil_moisture_grid)
dv_grid[:,:,1] = diversity_grid
plot(dv_grid, c = :reds)

GeoArrays.write("diversity_grid.tif", dv_grid)

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
