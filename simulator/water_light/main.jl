
##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- main.jl
##
## Contains exploratory analyses, accuracy checks for calculations,
## and miscellanious pieces that don't belong in a larger file.
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
Nspp::Int64 = 10; Nyr::Int64 = 1000; Ninit::Float64 = 1.0;

## ecological parameters
μ::Float64 = 0.1   ## mortality rate
E::Float64 = 0.5   ## evapotranspiration rate
l::Float64 = 1.5   ## leaf area allometric constant
b::Float64 = 3.0   ## biomass allometric constant
F::Float64 = 10.0  ## fecundity per unit biomass
W₀::Float64 = 0.4  ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity
P::Float64 = 10.0

## include function headers
include("utility_functions.jl")
include("simulation_functions.jl")
include("eq_functions.jl")
include("meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

##---------------------------------------------------------------
## Does single-species invasion criteria work?
##---------------------------------------------------------------

function calc_ir(g1, F, μ)
    (2 * F * g1^2) / (μ^3)
end

spp_data = generate_spp_data(1, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

## should successfully invade
μ = 0.1
calc_ir(spp_data.C₁[1], F, μ) ## invasion growth rate > 1
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 0.1, zeros(1,1), true, 3.0, 0.00005, false)
plot_simulation_dynamics(out)
## invasion is successful

## should not successfully invade
μ = 0.3
calc_ir(spp_data.C₁[1], F, μ)
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 0.1, zeros(1,1), true, 3.0, 0.00005, false)
plot_simulation_dynamics(out)
## invasion unsuccessful

out[10]

## check over range of parameters
μ_list = [0.1:0.01:0.4;]
ir_list = Vector{Float64}(undef, length(μ_list))
eq_n = Vector{Float64}(undef, length(μ_list))
for i in 1:length(μ_list)
    ir_list[i] = calc_ir(spp_data.C₁[1], F, μ_list[i])
    println(ir_list[i])
    out = sim_ppa(spp_data, 200, nrow(spp_data), Ninit, μ_list[i], F, 0.1, zeros(1,1), true, 3.0, 0.00005, false)
    eq_n[i] = out[2][200*10,2]
end

plot(log.(ir_list), eq_n, seriestype = :scatter, legend = :none, color = :black,
     frame = :box, ylab = "equilibrium pop. density", xlab = "log invasion growth rate")
vline!([0])
savefig("figures/inv_growth_ppa.pdf")

## invasion critera appears to work!

##---------------------------------------------------------------
## Is the canopy closure boundary stable?
##---------------------------------------------------------------
Nyr = 300
P = 10
Random.seed!(3)
spp_data = generate_spp_data(Nspp = 1, Wmax = 0.6, n_ht = 1, T = 1.0 / P,
                             F = F, μ = μ, b = 3.0, tradeoff_exp = 0.5,
                             tradeoff_sd = 0.0,
                             C₁max = 0.0001, C₂max = 0.00005);

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.5, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, 0.1)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.7, θ_fc, zeros(1,1),
                    true, 0.4, 3.0, 0.1, true, 0.8, 1.3)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.5, θ_fc, zeros(1,1),
                    true, 0.4, 3.0, 0.1, true, 0.8, 0.7)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.5, θ_fc, zeros(1,1),
                   true, 0.4, true, 0.8, 1.3)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 5.5, θ_fc, zeros(1,1),
                    true, 0.4, true, 0.8, 0.7)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.85, θ_fc)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.2, θ_fc)
plot_simulation_dynamics(out)
plot_canopy_cover(out)

maximum(Matrix(out[8][:,:]))


##---------------------------------------------------------------
## Do the strigul stability criteria work?
##---------------------------------------------------------------
Nyr = 200

## calculate x* by numerically
## solving equation D-7 from strigul 2008, modified to incorporate variation in T/P
function calc_xstar(g1, uf, P, μ, F)

    function xs(x_star)
        x = [1:1:1e4;]
        for i in 1:length(x)
            x[i] = exp(-μ * x[i] * (1.0/P)) * (x_star + g1*x[i]*(1.0/P))^(b-1)
        end
        F * (1.0/P) * exp(-(μ/(g1*uf)) * x_star) * sum(x) - 1
    end

    if xs(0.0) < 0.0
        return 0.0
    else
        return find_zero(xs, (0.0, 10.0))
    end
end

## calculate the stability criterion from equation 28 in strigul (doesn't appear to work??)
function stability_crit(spp_data, g1, understory_factor, P, μ, F)
    xstar = calc_xstar(spp_data.C₁[1], understory_factor, P, μ, F)
    (μ / (spp_data.C₁[1] * understory_factor)) * xstar
end


## do some exploration, varying parameter values to look at stability
spp_data = generate_spp_data(Nspp = 1, Wmax = 0.6, n_ht = 1,
                             T = 1.0 / P, F = F, μ = μ, b = 3.0,
                             tradeoff_exp = 0.5, tradeoff_sd = 0.0,
                             C₁max = 0.0001, C₂max = 0.00005)
spp_data.ht .= 0.5

F = 100.0
μ = 0.1
P = 10
Nyr = 200
uf = 0.1
μ = 0.1

stability_crit(spp_data, spp_data.C₁[1], uf, P, μ, F)
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, true, 0.5, 0.95)

plot(1:length(out[7]), out[7])
plot_simulation_dynamics(out)

xstar = calc_xstar(spp_data.C₁[1], uf, P, μ, F)
sqrt(xstar)

function ts(x)
    F * ((x^2 / μ) + ((2 * spp_data.C₁[1] * x) / μ^2) + ((2 * spp_data.C₁[1]^2) / μ^3)) -
        exp((μ / (spp_data.C₁[1] * uf)) * x)
end
sqrt(find_zero(ts, 1.0))

## stability collapses as the understory growth rate declines relative to overstory
## When this happens: understory growth doesn't compensate for overstory mortality, leading
## the canopy to occasionally open, crashing H* to zero, and then building up again.
##
## However, the system still "returns" to these oscillations after a perturbation, it doesnt
## just go to another equilibrium. So its all a bit unclear to me.

## stable
uf = 0.1
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)
plot(1:length(out[7]), out[7])

## stable
uf = 0.01
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)
plot(1:length(out[7]), out[7])

## unstable ?
uf = 0.0001
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)
plot(1:length(out[7]), out[7])

## try and get a sense of the pattern over a range of parameter values
uf_list = range(0.00001, stop = 0.025, length = 20)
var = Vector{Float64}(undef, length(uf_list))
ind = Vector{Float64}(undef, length(uf_list))

spp_data = generate_spp_data(1, 0.6, 1, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)
spp_data.ht .= 0.5
for i in 1:length(uf_list)
    ind[i] = stability_crit(spp_data, spp_data.C₁[1], uf_list[i], P, μ, F)
    out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 0.1, zeros(1,1), true, 3.0, uf_list[i], false)
    var[i] = std(out[7][Int(round(Nyr * 1.0/0.1))-100:Int(round(Nyr * 1.0/0.1))])
end

plot(uf_list, var, seriestype = :scatter)
plot(ind, var, seriestype = :scatter)

## Though I feel like I understand which parameters generate instability,
## I am not satisfied because I cannot get the threshold to work. Not sure
## if this is because the discrete nature of the simulator obscures things,
## or if something is very wrong...

##---------------------------------------------------------------
## Does the invasion criterion work?
##---------------------------------------------------------------

## first for light only:

F = 100.0
μ = 0.1
P = 10
Nyr = 300
uf = 0.1
μ = 0.1

Random.seed!(1)
spp_data = generate_spp_data(2, 0.6, 1, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

calc_xstar(spp_data[1,:C₁], uf, P, μ, F)
calc_xstar(spp_data[2,:C₁], uf, P, μ, F)
## species 1 should win

out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)

plot_simulation_dynamics(out)
plot(1:length(out[7]), out[7])
## and it does

## then for light and water:

F = 100.0
μ = 0.1
P = 10
Nyr = 600
uf = 0.1
μ = 0.1
mean_p = 2.0

Random.seed!(4) ## change seed to get iteration where both species have positive zstar
spp_data = generate_spp_data(2, 0.6, 1, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)
mono_zstar(spp_data, F, P, μ, mean_p, θ_fc, Nyr)
## species 1 should win

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, mean_p,
                    θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

spp_data

plot_simulation_dynamics(out)
plot(1:length(out[9]), out[9])
plot_canopy_cover(out)
## and species 1 does in fact win


## try for 1 species
Random.seed!(1)
spp_data = generate_spp_data(5, 0.6, 1, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

F = 100.0
μ = 0.1
P = 10
Nyr = 300
uf = 0.1
μ = 0.1
mean_p = 2.0

mono_zstar(spp_data, F, P, μ, mean_p, θ_fc, Nyr)
## species 2 should win (will differ for different runs)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, mean_p,
                    θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

plot_simulation_dynamics(out)
plot(1:length(out[9]), out[9])
plot_canopy_cover(out)
## and species 2 does appear to be on track to win (though it will take a long time for species 1 to go extinct)

##---------------------------------------------------------------
## Can multiple water strategies coexist with single height strategy?
##---------------------------------------------------------------

# illustrate with four species
Nyr = 500
P = 10
μ = 0.11
uf = 0.1
Random.seed!(4)
spp_data = generate_spp_data(4, 0.7, 1, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 11.0)
spp_data.ht .= 0.5

out_02 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.2, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_02)
plot_canopy_cover(out_02)

out_025 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.25, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_025)
plot_canopy_cover(out_025)

out_03 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_03)
plot_canopy_cover(out_03)

out_07 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.7, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_07)
plot_canopy_cover(out_07)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)

out_15 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_15)
plot_canopy_cover(out_15)

plot(plot_canopy_cover(out_02),
     plot_canopy_cover(out_025), plot_canopy_cover(out_03),
     plot_canopy_cover(out_07), plot_canopy_cover(out_1),
     plot_canopy_cover(out_15),
     legend = :none, layout = (3,2))
savefig("figures/sequence.pdf")


spp_data = generate_spp_data(15, 0.7, 1, 1.0 / P, F, μ, 3.0, 0.4, 0.0, 0.0001, 0.00005, 11.0)

p1 = plot(spp_data.Wᵢ[3:12], spp_data.C₁[3:12], marker_z = spp_data.C₁[3:12], xflip = true, seriestype = :scatter, legend = :none,
          markersize = 8, seriescolor = my_cgrad, frame = :box)

sd1 = DataFrame(spp_data[1,:])
sd1.spp .= 1
out_1 = sim_water_ppa(sd1, 400, nrow(sd1), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, false, 0.7, 0.9, true, 0.5)
plot_simulation_dynamics(out_1)
sd2 = DataFrame(spp_data[2,:])
sd2.spp .= 1
out_2 = sim_water_ppa(sd2, 400, nrow(sd2), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, false, 0.7, 0.9, true, 0.5)
plot_simulation_dynamics(out_2)

sd3 = DataFrame(spp_data[3,:])
sd3.spp .= 1
out_3 = sim_water_ppa(sd3, 400, nrow(sd3), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, false, 0.7, 0.9, true, 0.5)
plot_simulation_dynamics(out_3)

pdata = stack(out_1[2])
pdata.variable = parse.(Int64, pdata.variable)

pdata2 = stack(out_2[2])
pdata2.variable .= 2
pdata = vcat(pdata, pdata2)

pdata3 = stack(out_3[2])
pdata3.variable .= 3
pdata = vcat(pdata, pdata3)

pdata = pdata[pdata.rowkey .> 2000, :]
p1 = plot(pdata.rowkey ./ 10, pdata.value, group = pdata.variable, line_z = pdata.variable,
         xlim = [200, maximum(pdata.rowkey) / 10], ylim = [0, maximum(pdata.value)+10],
         seriescolor = my_cgrad, seriestype = :line, legend = :none, colorbar = :none,
         frame = :box, grid = false, linewidth = 4.5)



out_02 = sim_water_ppa(spp_data, 800, nrow(spp_data), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, false, 0.7, 0.9, true, 0.5, false, 0.7)
plot_simulation_dynamics(out_02)

pdata = stack(out_02[2])
pdata.variable = parse.(Int64, pdata.variable)
pdata = pdata[pdata.rowkey .> 4000, :]
p2 = plot(pdata.rowkey ./ 10, pdata.value, group = pdata.variable, line_z = pdata.variable,
         xlim = [400, maximum(pdata.rowkey) / 10], ylim = [0, maximum(pdata.value)+10],
         seriescolor = my_cgrad, seriestype = :line, legend = :none, colorbar = :none,
         frame = :box, grid = false, linewidth = 4.5)

plot(p1, p2)


uf = 0.1
spp_data = generate_spp_data(3, 0.7, 1, 1.0 / P, F, μ, 3.0, 0.4, 0.0, 0.0001, 0.00005, 11.0)

plot(spp_data.τ, spp_data.Wᵢ, marker_z = spp_data.C₁, seriestype = :scatter, legend = :none,
          markersize = 8, seriescolor = my_cgrad, frame = :box)

calc_eqN(spp_data, 1.0/P, 0.03, F, E, θ_fc, μ)
out_02 = sim_water_ppa(spp_data, 200, nrow(spp_data), Ninit, μ, F, P, 0.3, 0.7, zeros(1,1),
                       true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_02)
plot_canopy_cover(out_02)

pdata = stack(out_02[2])
pdata.variable = parse.(Int64, pdata.variable)
pdata
pdata = pdata[[8001*3:1:8001*12;], :]

p2 = plot(pdata.rowkey ./ 10, pdata.value, group = pdata.variable, line_z = pdata.variable,
         xlim = [0, maximum(pdata.rowkey) / 10],
         seriescolor = my_cgrad, seriestype = :line, legend = :none, colorbar = :none,
         frame = :box, grid = false, linewidth = 4.5)


## repeat but with perturbations to test stability
Nyr = 700
out_02 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.2, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_02)
plot_canopy_cover(out_02)

out_025 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.25, θ_fc, zeros(1,1),
                        true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_025)
plot_canopy_cover(out_025)

out_03 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.3, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_03)
plot_canopy_cover(out_03)

out_07 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.7, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_07)
plot_canopy_cover(out_07)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc, zeros(1,1),
                      true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)

out_15 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.5, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.6)
plot_simulation_dynamics(out_15)
plot_canopy_cover(out_15)

plot(plot_canopy_cover(out_02),
     plot_canopy_cover(out_025), plot_canopy_cover(out_03),
     plot_canopy_cover(out_07), plot_canopy_cover(out_1),
     plot_canopy_cover(out_15),
     legend = :none, layout = (3,2))
savefig("figures/sequence_perturb.pdf")


##---------------------------------------------------------------------
## Equilibrium population density - without transpiration in understory
##---------------------------------------------------------------------

Nyr = 900
μ = 0.31
F = 200.0
uf = 0.1

ρ = 0.015
P = 3
T = 1 / P
Random.seed!(4)
spp_data = generate_spp_data(10, 0.7, 1, 1.0 / P, F, μ, 3.0, 0.4, 0.0, 0.0001, 0.00005, 11.0, 0.3, 0.7)
spp_data = spp_data[1:3,:]
spp_data.ht .= 0.5

## first try it out on a small community
out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), [1.0, 1.0, 1.0], μ, F, P, ρ * P, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)
plot_canopy_cover(out)
plot_simulation_dynamics(out)

calc_eqN(spp_data, 1/P, ρ, F, E, θ_fc, μ, false, false)
calc_eq_biomass(DataFrame(spp_data), ρ, E, T, F, μ, uf)


##---------------------------------------------------------------
## calculating equilibrium biomass
##---------------------------------------------------------------

Nyr = 1500
μ = 0.31
F = 200.0
uf = 0.1

ρ = 0.04
P = 3
T = 1 / P
Random.seed!(4)
spp_data = generate_spp_data(10, 0.7, 1, 1.0 / P, F, μ, 3.0, 0.4, 0.0, 0.0001, 0.00005, 11.0, 0.3, 0.7)
spp_data = spp_data[1:5,:]
spp_data.ht .= 0.5

feas(spp_data, T, ρ, F, μ, θ_fc, uf, false, false)

## first try it out on a small community
out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), repeat([1.0], inner = nrow(spp_data)), μ, F, P, ρ * P, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)
plot_canopy_cover(out)
plot_simulation_dynamics(out)

calc_eqN(spp_data, 1/P, ρ, F, E, θ_fc, μ, false, false)
eqB = calc_eq_biomass(DataFrame(spp_data), ρ, E, T, F, μ, uf)

plot(eqB, vec(collect(out[12][4500, 2:6])), seriestype = :scatter)
plot!([0.0, maximum(eqB)], [0.0, maximum(eqB)], color = :black)



##---------------------------------------------------------------
## attempting to calculate invasion condition
##---------------------------------------------------------------

## first calculate the value of ρ at which the canopy first closes
## should be a function only of the characteristics of species 1

function calc_ρ_cc(T, spp_data, F, μ, b)
    E * spp_data.τ[1] + spp_data.Wᵢ[1] - spp_data.Wᵢ[nrow(spp_data)]
end

calc_ρ_cc(0.1, spp_data, F, μ, b)

function calc_T_cc(ρ, spp_data, F, μ, b)
    function st(T)
        [calc_τ(spp_data.C₁[1], spp_data.C₂[1], F, μ, T[1], b) - (spp_data.Wᵢ[nrow(spp_data)] + ρ - spp_data.Wᵢ[1]) / E]
    end

    nlsolve(st, [0.1])
end

calc_T_cc(0.01, spp_data, F, μ, b)

## when ρ or T are above their critical value, all species with an allometric height exponent
## shorter than the max are immediately excluded.

## after canopy closure, as ρ increases or T decreases, late season species will be sequentially
## excluded once they can no longer grow at g*

function calc_ρ_cx(spp_data, E, T)
    spp_data.Wᵢ[1] - spp_data.Wᵢ[nrow(spp_data)-1] + E * T *
        ((spp_data.C₁[nrow(spp_data)] + spp_data.C₂[nrow(spp_data)]) / (spp_data.C₁[1] + spp_data.C₂[1]))
end

function calc_T_cx(spp_data, E, ρ)
    ((spp_data.Wᵢ[nrow(spp_data)-1] + ρ - spp_data.Wᵢ[1]) / E) *
        ((spp_data.C₁[1] + spp_data.C₂[1]) / (spp_data.C₁[nrow(spp_data)] + spp_data.C₂[nrow(spp_data)]))
end


spp_data = generate_spp_data(10, 0.7, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 11.0)
spp_data.ht[[1,3,5,7,9]] .= 0.2
spp_data.ht[[2,4,6,8,10]] .= 0.8

feas(spp_data, 0.1, 0.04, F, μ, θ_fc, uf, false, false)

Nyr = 600
out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), repeat([1.0], nrow(spp_data)), μ, F, P, 0.4, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)
plot_canopy_cover(out)
plot_simulation_dynamics(out)

out = sim_water_ppa(DataFrame(spp_data[1,:]), Nyr, 1, [1.0], μ, F, P, 0.4, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)

out = sim_water_ppa(DataFrame(spp_data[3,:]), Nyr, 1, [1.0], μ, F, P, 0.4, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)

mono_zstar(spp_data, F, P, μ, 0.4, θ_fc, Nyr)

calc_xss(spp_data.C₁[4], spp_data.C₂[2], E, 0.1, 0.04, F, μ, uf, false)

##---------------------------------------------------------------
## Check if abundance calculations work
##---------------------------------------------------------------

ceq = check_eq_agreement(10, 4, 500, 0.4, 10.0, 0.4, μ, F, uf)

plot(1:maximum(ceq.eq_sim), 1:maximum(ceq.eq_sim))
plot!(ceq.eq_sim, ceq.eq_an, seriestype = :scatter)

ceq = check_eq_agreement(10, 4, 700, 0.4, 10.0, 0.4, μ, F, uf)

plot(1:700, 1:700, color = :lightgray, linewidth = 2.5,
     frame = :box, xlim = [0.0, 700.0], ylim = [0.0, 700.0], legend = :none)
plot!(ceq.eq_sim, ceq.eq_an, seriestype = :scatter, color = :black, xlab = "eq. abundance (simulated)",
      ylab = "eq. abundance (calculated)")
savefig("figures/eq_ab_check.pdf")

##---------------------------------------------------------------
## SCRATCH FROM FEASIBILITY
##---------------------------------------------------------------

Nyr = 400
μ = 0.31
F = 200.0
uf = 0.1

Random.seed!(4)
spp_data = generate_spp_data(4, 0.7, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 11.0)
spp_data = spp_data[1:3,:]
spp_data.ht .= 0.5

#spp_data = generate_spp_data(2, 0.6, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 20.0)
#spp_data.ht .= 0.5

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), [100.0, 100.0, 100.0], μ, F, P, 0.43, θ_fc,
                    zeros(1,1), true, 0.4, 3.0, uf, false)
plot_canopy_cover(out)
plot_simulation_dynamics(out)
spp_data

plot(1:length(out[9]), out[9], linewidth = 2.5)
plot(1:length(out[6]), out[6], linewidth = 2.5, ylim = [w3+0.0423, w2+0.043])

## There should be a constant growth rate for coexisting water strategies, if fecundity is the same.
T = (1.0 / P)
c1 = spp_data.C₁[1]
c2 = spp_data.C₂[1]
c22 = spp_data.C₁[2]
c3 = spp_data.C₁[3]
w0 = out[6][Int(round(Nyr*P))]
w1 = spp_data.Wᵢ[1]
w2 = spp_data.Wᵢ[2]
w3 = spp_data.Wᵢ[3]

w0 = w3 + 0.475 / P

w0 = w3 + 0.04424
w3 + 0.0443
w0 = w2 + 0.04425

function lrs3(gz)

    function int(x)
        F * exp(-μ * (gz[2] / (uf * gz[1]))) * exp(-μ * x) *
            (gz[2] + gz[1]*x)^2
    end

    xx = [0.0:gz[1]*uf*1e-4:gz[2];]
    s = sum(F .* 1e-4 * exp.(-μ .* (xx ./ (uf * gz[1]) )) .* xx .^ 2)
    [quadgk(int, 0.0, 1e4)[1] - 1.0,
     E * ((T * (gz[1] + c2)) / (c1 + c2)) *
         (s + 1.0) - (w0 - w1)]
end

sol3 = nlsolve(lrs3, [0.02, 0.0125], iterations = 5000)
sol3.zero[1]
sqrt(sol3.zero[2])

function calc_τ_cc(g, spp_data)
    (T .* (g .+ spp_data.C₂)) ./ (spp_data.C₁ .+ spp_data.C₂)
end

tau = calc_τ_cc(sol3.zero[1], spp_data)

xx = [0.0:sol3.zero[1]*uf*1e-4:sol3.zero[2];]
s = sum(F .* 1e-4 * exp.(-μ .* (xx ./ (uf * sol3.zero[1]) )) .* xx .^ 2)

(w2 - w3) / (E * (tau[3] - tau[2]))
(w1 - w2) / (E * (tau[2] - tau[1]))
(w0 - w1) / (E * (tau[1]))


## THE SOLUTION !!
ρ = 0.042;

function crit_ρ(r)
    g = ((w2 + r[1] - w1) / (E * T)) * (c1 + c2) - c2
    function lrs5(gz)
        function int(x)
            F * exp(-μ * (gz[1] / (uf * g))) * exp(-μ * x) *
                (gz[1] + g*x)^2
        end

        [quadgk(int, 0.0, 1e4)[1] - 1.0]
    end;
    sol5 = nlsolve(lrs5, [0.1])

    function int2(x)
        F * exp(-μ * (sol5.zero[1] / (uf * c3))) * exp(-μ * x) *
            (sol5.zero[1] + c3*x)^2
    end

    [quadgk(int2, 0.0, 1e4)[1] - 1.0]
end

sol = nlsolve(crit_ρ, [0.042])

g = ((w2 + 0.04249 - w1) / (E * T)) * (c1 + c2) - c2
function lrs5(gz)
    function int(x)
        F * exp(-μ * (gz[1] / (uf * g))) * exp(-μ * x) *
            (gz[1] + g*x)^2
    end

    [quadgk(int, 0.0, 1e4)[1] - 1.0]
end;
sol5 = nlsolve(lrs5, [0.1])

function int2(x)
    F * exp(-μ * (sol5.zero[1] / (uf * c3))) * exp(-μ * x) *
        (sol5.zero[1] + c3*x)^2
end

quadgk(int2, 0.0, 1e4)[1]
