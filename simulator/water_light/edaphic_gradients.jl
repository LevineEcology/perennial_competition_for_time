##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- edaphic_gradients.jl
##
## Diversity effects of precipitation gradients.
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
