## WATER_ONLY_SIMULATOR -- main.jl
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
P::Int64 = 10

## include function headers
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
μ = 0.2
calc_ir(spp_data.C₁[1], F, μ)
out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 0.1, zeros(1,1), true, 3.0, 0.00005, false)
plot_simulation_dynamics(out)
## invasion unsuccessful

## check over range of parameters
μ_list = [0.1:0.01:0.3;]
ir_list = Vector{Float64}(undef, length(μ_list))
eq_n = Vector{Float64}(undef, length(μ_list))
for i in 1:length(μ_list)
    ir_list[i] = calc_ir(spp_data.C₁[1], F, μ_list[i])
    println(ir_list[i])
    out = sim_ppa(spp_data, 200, nrow(spp_data), Ninit, μ_list[i], F, 0.1, zeros(1,1), true, 3.0, 0.00005, false)
    eq_n[i] = out[2][200*10,2]
end

plot(ir_list, eq_n, seriestype = :scatter)

## invasion critera appears to work!

##---------------------------------------------------------------
## Is the canopy closure boundary stable?
##---------------------------------------------------------------
Nyr = 300
P = 10
spp_data = generate_spp_data(1, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

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
spp_data = generate_spp_data(1, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)
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
uf_list = range(0.00001, stop = 0.1, length = 10)
var = Vector{Float64}(undef, length(uf_list))
ind = Vector{Float64}(undef, length(uf_list))

spp_data = generate_spp_data(1, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)
spp_data.ht .= 0.5
for i in 1:length(uf_list)
    ind[i] = stability_crit(spp_data, spp_data.C₁[1], uf_list[i], P, μ, F)
    out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 0.1, zeros(1,1), true, 3.0, uf_list[i], false)
    var[i] = std(out[7][Int(round(Nyr * 1.0/0.1))-100:Int(round(Nyr * 1.0/0.1))])
end

plot(uf_list, var, seriestype = :scatter)
plot(ind, var, seriestype = :scatter)

## Though I feel like I understand which parameters generate instability,
## I am not very satisfied because I cannot get the threshold to work. Not sure
## if this is because the discrete nature of the simulator obscures things,
## or if something is very wrong...

##---------------------------------------------------------------
## Does the invasion criterion work?
##---------------------------------------------------------------

## first for light only:
spp_data = generate_spp_data(2, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

F = 100.0
μ = 0.1
P = 10
Nyr = 300
uf = 0.1
μ = 0.1

calc_xstar(spp_data[1,:C₁], uf, P, μ, F)
calc_xstar(spp_data[2,:C₁], uf, P, μ, F)
## species 1 should win

out = sim_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1.0/P, zeros(1,1), true, 3.0, uf, false)

plot_simulation_dynamics(out)
plot(1:length(out[7]), out[7])
## and it does


## then for light and water:
spp_data = generate_spp_data(2, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

F = 100.0
μ = 0.1
P = 10
Nyr = 300
uf = 0.1
μ = 0.1
mean_p = 2.0

mono_zstar(spp_data, F, P, μ, mean_p, θ_fc, Nyr)
## species 1 should win (will differ for different runs)

out = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, mean_p,
                    θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

plot_simulation_dynamics(out)
plot(1:length(out[9]), out[9])
plot_canopy_cover(out)
## and species 1 does in fact win


## try for 5 species

## then for light and water:
spp_data = generate_spp_data(5, 0.6, 1.0 / P, F, μ, 3.0, 0.5, 0.0, 0.0001, 0.00005)

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
spp_data = generate_spp_data(4, 0.6, 1.0 / P, F, μ, 2.5, 0.5, 0.0, 0.0001, 0.00005)
spp_data.ht .= 0.5

out_05 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_05)
plot_canopy_cover(out_05)

out_07 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.7, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_07)
plot_canopy_cover(out_07)

out_085 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.85, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_085)
plot_canopy_cover(out_085)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)

out_12 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.2, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_12)
plot_canopy_cover(out_12)

out_15 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_15)
plot_canopy_cover(out_15)

out_175 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_175)
plot_canopy_cover(out_175)

out_20 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 2.0, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_20)
plot_canopy_cover(out_20)

out_25 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 2.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_25)
plot_canopy_cover(out_25)

plot(plot_canopy_cover(out_05), plot_canopy_cover(out_1),
     plot_canopy_cover(out_175), plot_canopy_cover(out_25), legend = :none)
savefig("figures/sequence.pdf")


## repeat but with perturbations to test stability

out_05 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.5, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.8)
plot_simulation_dynamics(out_05)
plot_canopy_cover(out_05)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc, zeros(1,1),
                      true, 0.4, 3.0, uf, true, 0.7, 0.8)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)

out_175 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.75, θ_fc, zeros(1,1),
                        true, 0.4, 3.0, uf, true, 0.7, 0.8)
plot_simulation_dynamics(out_175)
plot_canopy_cover(out_175)

out_25 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 2.5, θ_fc, zeros(1,1),
                       true, 0.4, 3.0, uf, true, 0.7, 0.8)
plot_simulation_dynamics(out_25)
plot_canopy_cover(out_25)

x1 = plot(plot_canopy_cover(out_05), plot_canopy_cover(out_1),
          plot_canopy_cover(out_175), plot_canopy_cover(out_25), legend = :none)
x2 = plot(plot_simulation_dynamics(out_05), plot_simulation_dynamics(out_1),
          plot_simulation_dynamics(out_175), plot_simulation_dynamics(out_25), legend = :none,
          colorbar = :none)
plot(x1, x2)

savefig("figures/sequence.pdf")


##---------------------------------------------------------------
## Edaphic gradients
##---------------------------------------------------------------

# first try with 2 different light strategies

spp_data = generate_spp_data(10, 0.6, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 20.0)
spp_data[[1:2:nrow(spp_data);], :ht] .= 0.5
spp_data[[2:2:nrow(spp_data);], :ht] .= 0.65
spp_data

uf = 0.1
Nyr = 800
P = 5

precip_list = range(0.1, stop = 1.5, length = 100)
results = DataFrame(total_precip = precip_list,
                    n_coex = Vector{Int64}(undef, length(precip_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(precip_list)),
                    canopy_cover = Vector{Float64}(undef, length(precip_list)),
                    starting_swc = Vector{Float64}(undef, length(precip_list)),
                    ending_swc = Vector{Float64}(undef, length(precip_list)),
                    avg_height = Vector{Float64}(undef, length(precip_list)))

mt = mortality_table(Nyr*P, μ, repeat([1.0/P], inner = Nyr*P))

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

p1 = plot(results.total_precip, results.n_coex, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, ylab = "# spp", ylim = (0, nrow(spp_data)))
p2 = plot(results.total_precip, results.canopy_cover, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover")
p3 = plot(results.total_precip, results.avg_height, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "total precip.")
plot(p1, p2, p3, layout = (3,1))
savefig("figures/gradient.pdf")

plot(results.total_precip, results.ending_swc, seriestype = :line)
plot(results.total_precip, results.starting_swc, seriestype = :line)


res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.75,
                        0.6, mt, false, 0.4, 3.0, uf, false)
plot_simulation_dynamics(res)


## FOR P
Nyr = 500
p_list = range(2, stop = 20, length = 10)
p_results = DataFrame(storm_frequency = p_list,
                    n_coex = Vector{Int64}(undef, length(p_list)),
                    spp_id = Vector{Vector{Int64}}(undef, length(p_list)),
                    canopy_cover = Vector{Float64}(undef, length(p_list)),
                    starting_swc = Vector{Float64}(undef, length(p_list)),
                    ending_swc = Vector{Float64}(undef, length(p_list)),
                    avg_height = Vector{Float64}(undef, length(p_list)))
p_results.storm_frequency .= Int.(round.(p_results.storm_frequency))

Threads.@threads for i in 1:length(p_list)
    P = p_results[i, :storm_frequency]
    res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.3,
                        0.6, zeros(1,1), false, 0.4, 3.0, uf, false)
    p_results[i, :n_coex] = sum(Matrix(res[2])[Nyr*P, 2:nrow(spp_data)+1] .> 0.01) ## set threshold for persistance
    p_results[i, :spp_id] = findall(Matrix(res[2])[Nyr*P, 2:nrow(spp_data)+1] .> 0.01)
    p_results[i, :canopy_cover] = sum(Matrix(res[8])[Nyr*P, 2:nrow(spp_data)+1])
    p_results[i, :starting_swc] = res[6][Nyr*P]
    p_results[i, :ending_swc] = res[5][Nyr*P]
    p_results[i, :avg_height] = mean(Matrix(res[7])[Nyr*P, 2:nrow(spp_data)+1])
    println("completed iteration: " * string(i) * " of " * string(length(p_list)))
end
p_results

p1 = plot(p_results.storm_frequency, p_results.n_coex, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, ylab = "# spp", ylim = (0, nrow(spp_data)))
p2 = plot(p_results.storm_frequency, p_results.canopy_cover, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :green, ylab = "% canopy cover")
p3 = plot(p_results.storm_frequency, p_results.avg_height, seriestype = :line, linewidth = 3,
          frame = :box, legend = :none, color = :red, ylab = "avg. height", xlab = "total precip.")
plot(p1, p2, p3, layout = (3,1))


res = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1, 1.0,
                        0.6, zeros(1,1), false, 0.4, 3.0, uf, false)
plot_simulation_dynamics(res)
plot_canopy_cover(res)



## first gradient in total precip

mono_zstar(spp_data, F, P, μ, 0.5, θ_fc, Nyr)

out_05 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 0.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_05)
plot_canopy_cover(out_05)


mono_zstar(spp_data, F, P, μ, 1.0, θ_fc, Nyr)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.0, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)


mono_zstar(spp_data, F, P, μ, 1.75, θ_fc, Nyr)

out_175 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_175)
plot_canopy_cover(out_175)


mono_zstar(spp_data, F, P, μ, 2.5, θ_fc, Nyr)

out_25 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 2.5, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_25)
plot_canopy_cover(out_25)


mono_zstar(spp_data, F, P, μ, 5.0, θ_fc, Nyr)

out_5 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, 5.0, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_5)
plot_canopy_cover(out_5)


## next gradient in P

mono_zstar(spp_data, F, 20, μ, 1.75, θ_fc, Nyr)

out_20 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 20, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_20)
plot_canopy_cover(out_20)

out_20 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 20, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_20)
plot_canopy_cover(out_20)

out_15 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 15, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_15)
plot_canopy_cover(out_15)

out_10 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 10, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_10)
plot_canopy_cover(out_10)

mono_zstar(spp_data, F, 5, μ, 1.75, θ_fc, Nyr)

out_5 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 5, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_5)
plot_canopy_cover(out_5)

out_3 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 3, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_3)
plot_canopy_cover(out_3)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1, 1.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)


## vary p, for lower total rainfall

out_20 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 20, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_20)
plot_canopy_cover(out_20)

out_15 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 15, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_15)
plot_canopy_cover(out_15)

out_10 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 10, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_10)
plot_canopy_cover(out_10)

mono_zstar(spp_data, F, 5, μ, 1.75, θ_fc, Nyr)

out_5 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 5, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_5)
plot_canopy_cover(out_5)

out_3 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 3, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_3)
plot_canopy_cover(out_3)

out_1 = sim_water_ppa(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, 1, 0.75, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
plot_simulation_dynamics(out_1)
plot_canopy_cover(out_1)


## do full gradient, make plots of species traits and canopy closure








##---------------------------------------------------------------
## multi-equil
##---------------------------------------------------------------

mintotal = 0.1
maxtotal = 5
lengthtotal = 2
minP = 2
maxP = 5
lengthP = 2

multi_1 = multi_eq(15, 5, 400, 0.4,
                   0.01*10, 0.5*10, 4, 2, 20, 8,
                   F, μ)
summary_1 = summarize_multi_eq(multi_1)
sub = summary_1[summary_1.var .== "n",:]

plot(sub.T, sub.mean, group = sub.total, line_z = sub.total,
     seriescolor = my_cgrad,
     seriestype = :line,
     xlim = [minimum(sub.T), maximum(sub.T)],
     legend = :none, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
     colorbar = true, colorbar_title = "mean annual precip.")

sub = summary_1[summary_1.var .== "avg",:]

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
