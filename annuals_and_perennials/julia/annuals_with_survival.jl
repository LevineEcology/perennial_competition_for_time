
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
