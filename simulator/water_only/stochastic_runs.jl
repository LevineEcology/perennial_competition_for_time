using Pkg
Pkg.activate("../../..//home/jl104/.julia/perennial_simulators")
Pkg.instantiate()

using Plots, DataFrames, Distributions,
    LazySets, NLsolve, PlotThemes, CSV

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

out = multi_eq_variable_water(30, 20, 4000,
                              0.1, 0.6, 5,
                              0.0, 0.2, 5,
                              1.0, 40.0, 5,
                              0.0, 30.0, 5,
                              0.4)
CSV.write("simulator_runs.csv", out)

