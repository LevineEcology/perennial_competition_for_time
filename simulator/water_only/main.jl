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
    Profile, PProf

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

function dens_est_constant_water(Nspp::Int64 = 10, Niter::Int64 = 10,
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)

    W₀_list = collect(range(minW₀, stop = maxW₀, length = lengthW₀));
    T_list = collect(range(minT, stop = maxT, length = lengthT));
    params = collect(Base.product(W₀_list, T_list));

    spp_data = generate_spp_data(Nspp, 0.6) ## initialize spp_data
    full_results = Array{Float64}(undef, Nspp*Niter, length(params))
    Threads.@threads for i in 1:length(params)
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        for j in 1:Niter
            spp_data = generate_spp_data(Nspp, 0.8, params[i][2], 100.0, 0.1, 2.5, 0.4, 0.0)
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] = calc_eqN(spp_data, F, E, params[i][1])[:,:eqN]
        end
        full_results[:,i] = sub_results
    end

    return [full_results, Niter, params]

end

@time results = dens_est_constant_water(30,8);

nfeas = Vector{Float64}(undef, (size(results[1])[2]))
for i in 1:size(results[1])[2]
    nfeas[i] = count(results[1][:,i] .> 0.0) / results[2]

end

plot([results[3][i][j] for i = 1:100, j = 1:2][:,1], nfeas,
     group = [results[3][i][j] for i = 1:100, j = 1:2][:,2], seriestype = :scatter, ylim = [0,30], legend = :bottomright)

plot([results[3][i][j] for i = 1:100, j = 1:2][:,2], nfeas,
     group = [results[3][i][j] for i = 1:100, j = 1:2][:,1], seriestype = :scatter, ylim = [0,30], legend = :bottomright)

plot(results[1][:,12], seriestype = :histogram)












spp_data = generate_spp_data(Nspp, 0.8, 56.0, 100.0, 0.1, 2.5, 0.4, 0.0)
calc_eqN(spp_data, F, E, 0.1)



spp_data = generate_spp_data(8, 0.6)
Profile.Allocs.@profile out = sim_water_only(spp_data, 1000, 8, Ninit);
PProf.Allocs.pprof(from_c = false)

biomass_dynamics = out[1]
n_dynamics = out[2]
n_data = out[3]
check_eq(biomass_dynamics, 1e-7)
calc_eqN(spp_data, F)
Vector(out[3][nrow(out[3]),:])

plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

plot(stack(biomass_dynamics).rowkey, stack(biomass_dynamics).value.+1, group = stack(biomass_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, xlims = (0,Nyr), yscale = :log10)

plot(stack(n_dynamics).rowkey, stack(n_dynamics).value, group = stack(n_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, yscale = :log10, xlims = (0, Nyr))
