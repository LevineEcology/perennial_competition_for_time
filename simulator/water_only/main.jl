using Base: NonReshapedReinterpretArray
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
    Profile, PProf, PlotThemes

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


theme(:dark)
my_cgrad = cgrad(:acton)


x = :W₀
typeof(x)

typeof("x")

function plot_multi_eq(data::DataFrame, yvar::String = "n", xvar::Symbol = :W₀)

    if xvar == :W₀
        groupvar = :T
    else
        groupvar = :W₀
    end

    subdata = data[data.var .== yvar, :]
    p = plot(subdata[:,xvar], subdata.mean, group = subdata[:,groupvar], line_z = subdata[:,groupvar],
             fill_z = subdata[:,groupvar], ribbon = subdata.sd, seriescolor = my_cgrad,
             seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata[:,xvar]), maximum(subdata[:,xvar])],
             legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)
    return p

end;

plot_multi_eq(summary)

plot(subdata.T, subdata.mean, group = subdata.W₀, line_z = subdata.W₀, fill_z = subdata.W₀,
     ribbon = subdata.sd, seriescolor = my_cgrad,
     seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata.T), maximum(subdata.T)],
     legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)


subdata = summary[summary.var .== "avg", :]
plot(subdata.W₀, subdata.mean, group = subdata.T, line_z = subdata.T, fill_z = subdata.T,
     ribbon = subdata.sd, seriescolor = my_cgrad,
     seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata.W₀), maximum(subdata.W₀)],
     legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)

plot(subdata.T, subdata.mean, group = subdata.W₀, line_z = subdata.W₀, fill_z = subdata.W₀,
     ribbon = subdata.sd, seriescolor = my_cgrad,
     seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata.T), maximum(subdata.T)],
     legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)


subdata = summary[summary.var .== "min", :]
plot(subdata.W₀, subdata.mean, group = subdata.T, line_z = subdata.T, fill_z = subdata.T,
     ribbon = subdata.sd, seriescolor = my_cgrad,
     seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata.W₀), maximum(subdata.W₀)],
     legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)

plot(subdata.T, subdata.mean, group = subdata.W₀, line_z = subdata.W₀, fill_z = subdata.W₀,
     ribbon = subdata.sd, seriescolor = my_cgrad,
     seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata.T), maximum(subdata.T)],
     legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3)
