##---------------------------------------------------------------
## WATER COMPETITION - ANNUALS & PERENNIALS
## By: Jacob Levine
## December 2022
##---------------------------------------------------------------

import NLsolve, DataFrames, Roots, Distributions, Plots
using JuMP, HiGHS, Ipopt

## function to calculate (numerically solve) break-even time from spp traits (using Roots.jl)
function calc_tau!(spp_data::DataFrames.DataFrame, F::Float64, μ::Float64,
                   T₀::Float64, B::Float64 = 2.5)
    spp_data.tau = zeros(DataFrames.nrow(spp_data))
    for i ∈ 1:DataFrames.nrow(spp_data)
        f(t) = F * spp_data.G[i]^B * t^B - μ * (T₀ - t) - 1
        spp_data[i, :tau] = Roots.find_zero(f, (0, 1e4))
    end
    return nothing
end

## function to solve for equilibrium population
function solve_eqN(spp_data::DataFrames.DataFrame, init::Float64 = 10.0)

    Q = DataFrames.nrow(spp_data)

    model = Model(Ipopt.Optimizer)

    @variables(model, begin
                   ## population densities
                   x[1:Q] ≥ 0.0, (start = init)
                   end);

    ## empty objective -- no optimization
    @objective(model, Max, 1.0)

    ## define expressions for growing season length
    Σᵢ = Array{NonlinearExpression}(undef, Q)
    for i ∈ reverse(1:Q)
        if i == Q
            Σᵢ[i] = @NLexpression(model, x[i] * spp_data.G[i]^(B-1))
        else
            Σᵢ[i] = @NLexpression(model, Σᵢ[i+1] + x[i] * spp_data.G[i]^(B-1))
        end
    end

    τᵢ = Array{NonlinearExpression}(undef, Q)
    for i ∈ 1:Q
        if i == 1
            τᵢ[i] = @NLexpression(model, (W₀ - spp_data.Wᵢ[i]) / (E * Σᵢ[i]))
        else
            τᵢ[i] = @NLexpression(model, τᵢ[i-1] + (spp_data.Wᵢ[i-1] - spp_data.Wᵢ[i]) / (E * Σᵢ[i]))
        end
    end

    ## define constraints for solver
    for i ∈ 1:Q
        @NLconstraint(model, 0.0 ≤ τᵢ[i] ≤ T₀^B)
        @NLconstraint(model, 1.0 == F * spp_data.G[i]^B * τᵢ[i] - μ * (T₀ - τᵢ[i]^(1/B)))
    end

    ## solve and print solution
    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    """)
    return [value.(x), "$(termination_status(model))"]

end

F = 100.0
μ = 0.0001
T₀ = 100.0
B = 2.5
W₀ = 0.6
E = 0.1

function feasibility(niter::Int64 = 30, nspp::Int64 = 10)

    feas_tradeoff = Vector{Bool}(undef, niter)
    feas_notradeoff = Vector{Bool}(undef, niter)

    ## for both tradeoff and non-tradeoff, iterate and determine feasability of a 10 spp eq.
    for i ∈ 1:niter

        ## for tradoeff case
        spp_data = DataFrames.DataFrame(G = rand(Distributions.Uniform(0.001, 0.01), 10))
        calc_tau!(spp_data, F, μ, T₀)
        spp_data.Wᵢ = 0.6 .* exp.(-0.03 .* spp_data.tau)
        DataFrames.sort!(spp_data, :Wᵢ, rev = true)
        spp_data = spp_data[spp_data.tau .< T₀, :]

        out = solve_eqN(spp_data, 1.0)
        if out[2] == "LOCALLY_INFEASIBLE"
            feas_tradeoff[i] = false
        else
            feas_tradeoff[i] = true
        end

        ## for no tradoeff case
        spp_data = DataFrames.DataFrame(G = rand(Distributions.Uniform(0.001, 0.01),10))
        calc_tau!(spp_data, F, μ, T₀)
        spp_data.Wᵢ = 0.6 .* exp.(-0.03 .* spp_data.tau) +
            rand(Distributions.Normal(0, 0.05), DataFrames.nrow(spp_data))
        DataFrames.sort!(spp_data, :Wᵢ, rev = true)
        spp_data = spp_data[spp_data.tau .< T₀, :]

        out = solve_eqN(spp_data, 1.0)
        if out[2] == "LOCALLY_INFEASIBLE"
            feas_notradeoff[i] = false
        else
            feas_notradeoff[i] = true
        end

    end

    return [feas_tradeoff, feas_notradeoff]

end

fs_results = feasibility(30, 10)

sum(fs_results[1]) / length(fs_results[1])
sum(fs_results[2]) / length(fs_results[2])



spp_data = DataFrames.DataFrame(G = rand(Distributions.Uniform(0.001, 0.01), 10))
calc_tau!(spp_data, F, μ, T₀)
spp_data.Wᵢ = 0.6 .* exp.(-0.03 .* spp_data.tau) +
     rand(Distributions.Normal(0, 0.05), DataFrames.nrow(spp_data))
DataFrames.sort!(spp_data, :Wᵢ, rev = true)
spp_data = spp_data[spp_data.tau .< T₀, :]

Plots.plot(vcat(0.0, spp_data.tau), vcat(W₀, spp_data.Wᵢ), seriestype = :scatter, framestyle = :box,
           grid = false, color = :black, legend = false)

out = solve_eqN(spp_data, 10.0)
if all(out[1] .> 0.0)
    notdone = false
else
    notdone = true
end

while notdone
    spp_data = spp_data[out[1] .> 0.0, :]
    out = solve_eqN(spp_data, 10.0)
    if all(out[1] .> 0.0)
        notdone = false
    end
end
length(out[1])

spp_data

out = solve_eqN(spp_data, 10.0)
spp_data = (spp_data[out[1] .> 0.0, :])
out = solve_eqN(spp_data[out[1] .> 0.0, :], 10.0)
solve_eqN(spp_data[out[1] .> 0.0, :], 10.0)
