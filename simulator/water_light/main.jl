## WATER_AND_LIGHT_SIMULATOR -- main.jl
##
## by: Jacob Levine - jacoblevine@princeton.edu
## November 2022
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. define includes and parameters
##---------------------------------------------------------------
##---------------------------------------------------------------

using Plots, QuadGK, DataFrames, Distributions,
    SpecialFunctions, NLsolve, AutoPreallocation,
    Profile, PProf, PlotThemes, CSV, GeoArrays, Random, Debugger,
    ForwardDiff, Polylogarithms, ProgressBars, Roots

## simulation parameters
Nspp::Int64 = 10; Nyr::Int64 = 1000; Ninit::Float64 = 1.0;

## ecological parameters
μ::Float64 = 0.1 ## mortality rate
E::Float64 = 0.5 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 3.0  ## biomass allometric constant
F::Float64 = 10.0  ## fecundity per unit biomass
W₀::Float64 = 0.4 ## initial water content (default)
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


## try for 5 species
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

savefig(p1, "~/Documents/Science/application_materials/f3.pdf")
savefig(p2, "~/Documents/Science/application_materials/f4.pdf")



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

##---------------------------------------------------------------
## Edaphic gradients -- CALCULATIONS
##---------------------------------------------------------------

meq = multi_eq(30, 100, 0.2, 0.05, 0.5, 40, 3, 10, 4, F, μ, 2, 0.1)
smeq = summarize_multi_eq(meq)
smeq.P .= round.(smeq.P, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.total ./ 10.0, sub.mean, group = sub.P, line_z = sub.P, ribbon = sub.sd, fill_z = sub.P, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [0.5, 5.0], ylims = [0.0, 30.0],
          xlab = "Monthly precipitation (volumetric soil water equivalent)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass(meq)
smeq_b.P .= round.(smeq_b.P, digits = 2)
bp = plot(smeq_b.total ./ 10.0, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P, fill_z = smeq_b.P,
          linewidth = 3.5, ribbon = smeq_b.sd,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [0.5, 5.0],
          xlab = "Monthly precipitation (volumetric soil water equivalent)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/total_gradient.pdf")

summarize_multi_eq_transpiration(meq)

plot(sub.mean, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P, fill_z = smeq_b.P,
     linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Species Richness", ylab = "Total biomass")
savefig("figures/div_es_total.pdf")


meq = multi_eq(30, 100, 0.4, 0.2, 0.5, 4, 3, 30, 40, F, μ, 2)
smeq = summarize_multi_eq(meq)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.P, sub.mean, group = sub.total ./ 10.0, line_z = sub.total, ribbon = sub.sd, fill_z = sub.total, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [3.0, 30.0], ylims = [0.0, 30.0],
          xlab = "Storm frequency (storms per month)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass(meq)
bp = plot(smeq_b.P, smeq_b.mean, group = smeq_b.total, line_z = smeq_b.total,
          fill_z = smeq_b.total, linewidth = 3.5, ribbon = smeq_b.sd,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [3.0, 30.0],
          xlab = "Monthly precipitation (volumetric soil water equivalent)", ylab = "Total biomass")

plot(np, bp, layout = [1,1])
savefig("figures/P_gradient.pdf")

plot(sub.mean, log.(smeq_b.mean), group = smeq_b.total, line_z = smeq_b.total, fill_z = smeq_b.total ,
     linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Species richness", ylab = "Total biomass")
savefig("figures/div_es_P.pdf")



meq = multi_eq(30, 100, 0.2, 0.05, 0.5, 40, 10, 10, 1, F, μ, 2, 0.1)
smeq = summarize_multi_eq(meq)
smeq.P .= round.(smeq.P, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.total ./ 100.0, sub.mean, group = sub.P, line_z = sub.P, ribbon = sub.sd, fill_z = sub.P, linewidth = 3.5, legend = :none,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [0.05, 0.5], ylims = [0.0, 30.0], ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass(meq)
smeq_b.P .= round.(smeq_b.P, digits = 2)
bp = plot(smeq_b.total ./ 100.0, smeq_b.mean .* 5, group = smeq_b.P, line_z = smeq_b.P, fill_z = smeq_b.P,
          linewidth = 3.5, ribbon = smeq_b.sd, legend = :none,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, xlims = [0.05, 0.5],
          xlab = "Mean storm magnitude (volumetric soil water equivalent)", ylab = "GPP (g C m⁻² day⁻¹)")

plot(np, bp, layout = [1,1])
savefig("~/Documents/Science/application_materials/f1.pdf")

meq = multi_eq(30, 100, 0.3, 0.25, 0.25, 1, 1, 20, 40, F, μ, 2, 0.1)
smeq = summarize_multi_eq(meq)
smeq.P .= round.(smeq.P, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.P, sub.mean, group = sub.total, line_z = sub.total, ribbon = sub.sd, fill_z = sub.total, linewidth = 3.5, legend = :none,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box, ylims = [0.0, 30.0], ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass(meq)
smeq_b.P .= round.(smeq_b.P, digits = 2)
bp = plot(smeq_b[40:-1:1, :P], smeq_b.mean .* 5, group = smeq_b.total, line_z = smeq_b.total, fill_z = smeq_b.total,
          linewidth = 3.5, legend = :none, ribbon = smeq_b.sd,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Mean storm frequency (storms per month)", ylab = "GPP (g C m⁻² day⁻¹)")
plot(np, bp, layout = [1,1])
savefig("~/Documents/Science/application_materials/f2.pdf")



spp_data = generate_spp_data(15, 0.7, 1, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.0001, 11.0)
calc_eqN(spp_data, 0.01, 0.03)
calc_eq_biomass(spp_data, 0.03, E, 10, F, μ, 0.1, 0.4)
T = 0.1
ρ = 0.5


##---------------------------------------------------------------
## Edaphic gradients -- SIMULATIONS
##---------------------------------------------------------------

# first try with 2 different light strategies
spp_data = generate_spp_data(20, 0.6, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 20.0)
spp_data[[1:2:nrow(spp_data);], :ht] .= 0.5
spp_data[[2:2:nrow(spp_data);], :ht] .= 0.65
spp_data

uf = 0.1
Nyr = 800
P = 5.0

precip_list = range(0.1, stop = 1.5, length = 100)
results = DataFrame(total_precip = precip_list,
                    n_coex = Vector{Int64}(undef, length(precip_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(precip_list)),
                    canopy_cover = Vector{Float64}(undef, length(precip_list)),
                    starting_swc = Vector{Float64}(undef, length(precip_list)),
                    ending_swc = Vector{Float64}(undef, length(precip_list)),
                    avg_height = Vector{Float64}(undef, length(precip_list)))

mt = mortality_table(Nyr*Int(round(P)), μ, repeat([1.0/P], inner = Nyr*Int(round(P))))

Threads.@threads for i in 1:length(precip_list)
    res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, results[i, :total_precip],
                        0.6, mt, false, 0.4, 3.0, uf, false)
    results[i, :n_coex] = sum(Matrix(res[2])[Nyr*P, 2:nrow(spp_data)+1] .> 0.01) ## set threshold for persistance
    results[i, :spp_id] = findall(Matrix(res[2])[Nyr*P, 2:nrow(spp_data)+1] .> 0.01)
    results[i, :canopy_cover] = sum(Matrix(res[8])[Nyr*P, 2:nrow(spp_data)+1])
    results[i, :starting_swc] = res[6][Nyr*P]
    results[i, :ending_swc] = res[5][Nyr*P]
    results[i, :avg_height] = mean(Matrix(res[7])[Nyr*P, 2:nrow(spp_data)+1])
    println("completed iteration: " * string(i) * " of " * string(length(precip_list)))
end
results

CSV.write("gradient_results_map.csv", results)

results = CSV.read("gradient_results_map.csv", DataFrame)
results.total_precip .= results.total_precip .* 1000

p1 = plot(results.total_precip, results.n_coex, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :orange, ylab = "spp richness", ylim = (0, 20))
p2 = plot(results.total_precip, results.canopy_cover, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover")
p3 = plot(results.total_precip, results.avg_height, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "mean annual precipitation")
p4 = plot(results.total_precip, results.ending_swc, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :royalblue, ylab = "pre-storm SWC", xlab = "mean annual precipitation")
p5 = plot(results.total_precip, results.starting_swc, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :lightblue, ylab = "SWC", xlab = "mean annual precipitation")
plot(p1, p2, p4, layout = (3,1))
savefig("figures/gradient.pdf")


res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.2, 0.75,
                        0.6, mt, false, 0.4, 3.0, uf, false)
plot_simulation_dynamics(res)

## FOR P, maintaining storm size
Nyr = 500
p_list = range(4, stop = 12, length = 100)
p_results = DataFrame(storm_frequency = p_list,
                    n_coex = Vector{Int64}(undef, length(p_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(p_list)),
                    canopy_cover = Vector{Float64}(undef, length(p_list)),
                    starting_swc = Vector{Float64}(undef, length(p_list)),
                    ending_swc = Vector{Float64}(undef, length(p_list)),
                    avg_height = Vector{Float64}(undef, length(p_list)))

Threads.@threads for i in 1:length(p_list)
    P = p_results[i, :storm_frequency]
    res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.05 * P,
                        0.6, zeros(1,1), false, 0.4, 3.0, uf, false)
    p_results[i, :n_coex] = sum(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01) ## set threshold for persistance
    p_results[i, :spp_id] = findall(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01)
    p_results[i, :canopy_cover] = sum(Matrix(res[8])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    p_results[i, :starting_swc] = res[6][Nyr*Int(round(P))]
    p_results[i, :ending_swc] = res[5][Nyr*Int(round(P))]
    p_results[i, :avg_height] = mean(Matrix(res[7])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    println("completed iteration: " * string(i) * " of " * string(length(p_list)))
end
p_results
CSV.write("gradient_results_freq.csv", p_results)

Nyr * Int(round(p_results[i,:storm_frequency]))

p1 = plot(p_results.storm_frequency, p_results.n_coex, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, ylab = "# spp", ylim = (0, nrow(spp_data)))
p2 = plot(p_results.storm_frequency, p_results.canopy_cover, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover", ylim = (0.2, 1.0))
p3 = plot(p_results.storm_frequency, p_results.avg_height, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "storm frequency")
p4 = plot(p_results.storm_frequency[1:97], p_results.ending_swc[1:97], seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :red, ylab = "SWC", xlab = "storm frequency")
#p4 = plot!(p_results.storm_frequency[1:97], p_results.starting_swc[1:97], seriestype = :line, linewidth = 3,
#          frame = :box, legend = :none, color = :orange, ylab = "SWC", xlab = "storm frequency")
plot(p1, p2, p4, layout = (3,1))

savefig("figures/gradient_freq.pdf")


## FOR P, variable storm size
Nyr = 500
p_list = (range(0.75, stop = 25, length = 100))
pv_results = DataFrame(storm_frequency = p_list,
                    n_coex = Vector{Int64}(undef, length(p_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(p_list)),
                    canopy_cover = Vector{Float64}(undef, length(p_list)),
                    starting_swc = Vector{Float64}(undef, length(p_list)),
                    ending_swc = Vector{Float64}(undef, length(p_list)),
                    avg_height = Vector{Float64}(undef, length(p_list)))

Threads.@threads for i in 1:length(p_list)
    P = pv_results[i, :storm_frequency]
    res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.2,
                        0.2, zeros(1,1), false, 0.4, 3.0, uf, false)
    pv_results[i, :n_coex] = sum(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01) ## set threshold for persistance
    pv_results[i, :spp_id] = findall(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01)
    pv_results[i, :canopy_cover] = sum(Matrix(res[8])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    pv_results[i, :starting_swc] = res[6][Nyr*Int(round(P))]
    pv_results[i, :ending_swc] = res[5][Nyr*Int(round(P))]
    pv_results[i, :avg_height] = mean(Matrix(res[7])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    println("completed iteration: " * string(i) * " of " * string(length(p_list)))
end

pv_results
CSV.write("gradient_results_freq_novar.csv", pv_results)
plot_simulation_dynamics(res)
pv_results.ending_swc .= 0.086984

pv_results.storm_time .= 30 ./ pv_results.storm_frequency

p1 = plot(pv_results.storm_frequency, pv_results.n_coex, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :orange, ylab = "# spp", ylim = (0, nrow(spp_data)))
p2 = plot(pv_results.storm_frequency, pv_results.canopy_cover, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover", ylim = (0.2, 1.0))
p3 = plot(pv_results.storm_frequency, pv_results.avg_height, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "total precip.")
p4 = plot(pv_results.storm_frequency, pv_results.ending_swc, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :royalblue, ylab = "SWC", xlab = "storm frequency", ylim = (0.085, 0.105))
plot(p1, p2, p4, layout = (3,1))

pv_results

savefig("figures/gradient_freq_novar.pdf")


## FOR P, variable storm size
Nyr = 500
p_list = (range(0.75, stop = 25, length = 100))
pv_results = DataFrame(storm_frequency = p_list,
                    n_coex = Vector{Int64}(undef, length(p_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(p_list)),
                    canopy_cover = Vector{Float64}(undef, length(p_list)),
                    starting_swc = Vector{Float64}(undef, length(p_list)),
                    ending_swc = Vector{Float64}(undef, length(p_list)),
                    avg_height = Vector{Float64}(undef, length(p_list)))

P = 2

Threads.@threads for i in 1:length(p_list)
    P = pv_results[i, :storm_frequency]
    res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.5,
                        0.2, zeros(1,1), false, 0.4, 3.0, uf, false)
    pv_results[i, :n_coex] = sum(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01) ## set threshold for persistance
    pv_results[i, :spp_id] = findall(Matrix(res[2])[Nyr*Int(round(P)), 2:nrow(spp_data)+1] .> 0.01)
    pv_results[i, :canopy_cover] = sum(Matrix(res[8])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    pv_results[i, :starting_swc] = res[6][Nyr*Int(round(P))]
    pv_results[i, :ending_swc] = res[5][Nyr*Int(round(P))]
    pv_results[i, :avg_height] = mean(Matrix(res[7])[Nyr*Int(round(P)), 2:nrow(spp_data)+1])
    println("completed iteration: " * string(i) * " of " * string(length(p_list)))
end

pv_results
CSV.write("gradient_results_freq_novar_high.csv", pv_results)

plot_simulation_dynamics(res)
pv_results.ending_swc .= 0.086984

pv_results.storm_time .= 30 ./ pv_results.storm_frequency

p1 = plot(pv_results.storm_frequency, pv_results.n_coex, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :orange, ylab = "# spp", ylim = (0, nrow(spp_data)))
p2 = plot(pv_results.storm_frequency, pv_results.canopy_cover, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover", ylim = (0.2, 1.01))
p3 = plot(pv_results.storm_frequency, pv_results.avg_height, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "total precip.")
p4 = plot(pv_results.storm_frequency, pv_results.ending_swc, seriestype = :line, linewidth = 3, xflip = true,
          frame = :box, legend = :none, color = :royalblue, ylab = "SWC", xlab = "storm frequency", ylim = (0.085, 0.11))
plot(p1, p2, p4, layout = (3,1))

pv_results

savefig("figures/gradient_freq_novar_high.pdf")

# #---------------------------------------------------------------
## Stochastic simulations
##---------------------------------------------------------------

Nyr = 400
μ = 0.31
F = 200.0
uf = 0.1

Random.seed!(4)
spp_data = generate_spp_data(4, 0.7, 1, 1.0 / P, F, μ, b, 0.4, 0.0, 0.0001, 0.00005, 11.0)

P = 5.0
mp = 0.2

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

## start dialing up variation in mean annual precip
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.1, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1000.0, mp, 0.2, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1000.0, mp, 0.4, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


##---------------------------------------------------------------
## For low mean map (only latest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 5.0
mp = 0.05

r = generate_rainfall_regime(400, P, 1000.0, mp, 0.0, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])

## sd = 0.033
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.033, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


## sd = 0.1
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.1, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


## sd = 0.4
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.45, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)

plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


##---------------------------------------------------------------
## For high mean map (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 5.0
mp = 0.5

r = generate_rainfall_regime(400, P, 1000.0, mp, 0.0, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])

## sd = 0.13
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.13, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


## sd = 0.15
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.15, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


## sd = 0.5
r = generate_rainfall_regime(400, P, 1000.0, mp, 0.5, true, false)
rplot = plot_rainfall_regime(r)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)

plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])



##---------------------------------------------------------------
## For moderate mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 10.0
mp = 0.3

r = generate_rainfall_regime(400, P, 1000.0, mp, 0.0, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1e5, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1.0, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 0.01, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1e-3, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])



##---------------------------------------------------------------
## For small mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 3.0
mp = 0.6

r = generate_rainfall_regime(400, P, 1000.0, mp, 0.0, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1e5, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1.0, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 0.01, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(400, P, 1e-3, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])



##---------------------------------------------------------------
## For high mean P (only earliest species persists)
##---------------------------------------------------------------

## start dialing up variation in mean annual precip
P = 20.0
mp = 0.6

r = generate_rainfall_regime(200, P, 1000.0, mp, 0.0, true, false)
rplot = plot_rainfall_regime(r, 1:200)
length(r[1])

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(200, P, 1e5, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(200, P, 1.0, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(200, P, 0.01, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


r = generate_rainfall_regime(200, P, 1e-3, mp, 0.0, false, true, false)
rplot = plot_rainfall_regime(r, 1:200)

out = sim_water_ppa_stochastic(spp_data, length(r[1]), nrow(spp_data), 1.0, r, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot = plot_simulation_dynamics(out)
plot(plot(rplot, xlab = ""), plot(dynplot, colorbar = :none), layout = [1,1])


##---------------------------------------------------------------
## multi_eq methods
##
##---------------------------------------------------------------


out = multi_eq_variable_map(20, 10, 400, 0.4, 0.05, 0.4, 10, 0.0, 1.0, 3, 5.0, 4.0, 2, F, μ, false)
out[1]

smeq = summarize_multi_eq_variable_map(out)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.mapmean, sub.mean, group = sub.mapsd, line_z = sub.mapsd, ribbon = sub.sd, fill_z = sub.mapsd, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Storm frequency (storms per month)", ylab = "Species richness")

##---------------------------------------------------------------
## Alternative forcings
##---------------------------------------------------------------

##---------------------------------------------------------------
## carbon and MAP
##---------------------------------------------------------------

calc_aₘ(500, 1000)

x = collect(range, )
plot(x, calc_aₘ.(500, calc_vpd.(x, 30.0)))

x = collect(range(100.0, 1500.0, 100))
plot(x, calc_aₘ.(x, calc_vpd(20.0, 30.0)))
plot!(x, calc_aₘ.(x, calc_vpd(20.0, 30.0), 1.04545, 1000.0, 50.0, 10.0))

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
## temp/VPD
##---------------------------------------------------------------

meq = multi_eq_temp_map(50, 5, 0.2, 10.0, 45.0, 100, 0.2, 0.4, 5, 10.0, F, μ, 2, 0.1)
smeq = summarize_multi_eq_temp_map(meq)
smeq.t .= round.(smeq.t, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.t, sub.mean, group = sub.map ./ 10.0, line_z = sub.map ./ 10.0,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.map ./ 10.0, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_map(meq)
smeq_b.t .= round.(smeq_b.t, digits = 2)
bp = plot(smeq_b.t, (smeq_b.mean), group = smeq_b.map ./ 10.0, line_z = smeq_b.map ./ 10.0,
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
## temp and storm freq
##---------------------------------------------------------------

meq = multi_eq_temp_P(50, 30, 0.2, 10.0, 45.0, 100, 7.0, 15.0, 5, 0.3, F, μ, 2, 0.1)
smeq = summarize_multi_eq_temp_P(meq)
smeq.t .= round.(smeq.t, digits = 2)
sub = smeq[smeq.var .== "n", :]
np = plot(sub.t, sub.mean, group = sub.P, line_z = sub.P,
          ribbon = (sub.mean .- sub.lower, sub.upper .- sub.mean),
          fill_z = sub.P, linewidth = 3.5,
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Species richness")

smeq_b = summarize_multi_eq_biomass_temp_P(meq)
smeq_b.t .= round.(smeq_b.t, digits = 2)
bp = plot(smeq_b.t, (smeq_b.mean), group = smeq_b.P, line_z = smeq_b.P,
          fill_z = smeq_b.P,
          linewidth = 3.5,
          ribbon = (smeq_b.mean .- smeq_b.lower, smeq_b.upper .- smeq_b.mean),
          seriescolor = my_cgrad, fillalpha = 0.3, colorbar = false, frame = :box,
          xlab = "Average temperature (C)", ylab = "Log ecosystem biomass")

plot(np, bp, layout = [1,1])
savefig("figures/map_freq.pdf")
