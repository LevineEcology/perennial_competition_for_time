using Plots: _plots_plotly_defaults
##---------------------------------------------------------------
## WATER_ONLY_SIMULATOR -- eq_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## November 2022
##
## This script contains functions related to calculated equilibrium values
## and determining coexistence outcomes. All functions defined here are primarily
## for finding analytical solutions, rather than performing simulations
##---------------------------------------------------------------

"""
    minfeas(μ::Float64 = 0.1, F::Float64 = 1.0, b::Float64 = 2.5, analytical = false)

Calculates the minimum feasible growth rate, `g`, given a mortality rate `μ`, fecundity rate `F`, and
allometric exponent `b`. If `analytical` is `true`, the denominator of the growth rate expression is
calculated using a gamma function (as in the continuous time model). If `analytical` is `false`, the
denominator is calculated using the polylogarithm function (as in the discrete time model).
"""
function minfeas(μ::Float64 = 0.1, F::Float64 = 1.0, b::Float64 = 2.5, analytical = false)

    if analytical
        s = gamma(b, 0.0) * μ^-(b)
    else
        s = sigma(μ, b)
    end

   (1/(F * s))^(1/(b-1))

end


"""
    calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)

Calculates the equilibrium leaf area given a vector of equilibrium abundances `eqN.`
"""
function calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)
    eqN ./ (μ .* F)
end


"""
    μ_feas(spp_data::DataFrame, F::Float64, μ::Float64, analytical::Bool = false)

Determine which species in `spp_data` are feasible in a non-competitive context. Operates by calling
the `minfeas()` function, which calculates the minimum feasible growth rate given the mortality rate, `μ`,
fecundity rate, `F`, and allometric exponent `b`. This function returns a list of indices indicating which
species are feasible. Indices refer to species ordering in `spp_data`.
"""
function μ_feas(spp_data::DataFrame, F::Float64, μ::Float64, analytical::Bool = false)
    mf = minfeas(μ, F, b, analytical)
    findall(spp_data.C₁ .> mf)
end



"""
    calc_ρ_cc(T::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

Calculate the value of `ρ` (the minimum storm size at which the canopy closes) as a function of
the time between storm events `T`, species characteristics `spp_data`, fecundity `F`, mortality rate `μ`,
and allometric exponent `b`.

"""
function calc_ρ_cc(T::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)
    E * spp_data.τ[1] + spp_data.Wᵢ[1] - spp_data.Wᵢ[nrow(spp_data)]
end



"""
    calc_T_cc(ρ::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

Calculate the value of `T` (the time between storm events) as a function of
storm size `ρ`, species characteristics `spp_data`, fecundity `F`, mortality rate `μ`,
and allometric exponent `b`.

"""
function calc_T_cc(ρ::Float64, spp_data::DataFrame, F::Float64, μ::Float64, b::Float64)

    function st(T)
        [calc_τ(spp_data.C₁[1], spp_data.C₂[1], F, μ, T[1], b) - (spp_data.Wᵢ[nrow(spp_data)] + ρ - spp_data.Wᵢ[1]) / E]
    end

    nlsolve(st, [0.1])
end



"""
    calc_ρ_cx(spp_data::DataFrame, E::Float64, T::Float64)

Calculate the minimum value of `ρ` (storm size) at which the most drought tolerant species in
`spp_data` is excluded from the community due to capped transpiration.

"""
function calc_ρ_cx(spp_data::DataFrame, E::Float64, T::Float64)
    spp_data.Wᵢ[1] - spp_data.Wᵢ[nrow(spp_data)-1] + E * T *
        ((spp_data.C₁[nrow(spp_data)] + spp_data.C₂[nrow(spp_data)]) / (spp_data.C₁[1] + spp_data.C₂[1]))
end



"""
    calc_T_cx(spp_data::DataFrame, E::Float64, T::Float64)

Calculate the maximum value of `T` (storm interval length) at which themost drought tolerant species in
`spp_data` is excluded from the community due to capped transpiration.

"""
function calc_T_cx(spp_data::DataFrame, E::Float64, ρ::Float64)
    ((spp_data.Wᵢ[nrow(spp_data)-1] + ρ - spp_data.Wᵢ[1]) / E) *
        ((spp_data.C₁[1] + spp_data.C₂[1]) / (spp_data.C₁[nrow(spp_data)] + spp_data.C₂[nrow(spp_data)]))
end



"""
    calc_xss(C₁::Float64, C₂::Float64, E::Float64, T::Float64, ρ::Float64,
             F::Float64, μ::Float64, uf::Float64, understory_transpiration::Bool = false)

Calculate the equilibrium canopy closure size, `xss` or ``x^{**}`` in the manuscript. This value is calculated
from the characteristics of the least drought tolerant species (the species with the
highest value of ``W_{i}^{*}``) with a given height exponent, `h`.

"""
function calc_xss(C₁::Float64, C₂::Float64, E::Float64, T::Float64,
                  ρ::Float64, F::Float64, μ::Float64, uf::Float64, understory_transpiration::Bool = false)

    ## calculation method changes if understory species transpire soil water
    if understory_transpiration

    else

        ## determine if canopy closes
        if ρ / (E * calc_τ(C₁, C₂, F, μ, T, b)) < 1.0

            return 0.0

        else

            ## calculate eq growth rate
            gss = (ρ / (E * T)) * (C₁ + C₂) - C₂

            ## create function representing expression which must be solved numerically
            function cc(xs)

                ## integral within the expression
                function intcc(x)
                    exp(-μ * x) * (xs[1] + gss * x)^2
                end

                ## full expression
                [F * exp(-μ * (xs[1] / (uf * gss))) * quadgk(intcc, 0.0, 1e4)[1] - 1.0]

            end

            ## solve root numerically
            return nlsolve(cc, [0.01]).zero[1]

        end

    end

end



"""
    feas(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
    F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
    analytical::Bool = false, understory_transpiration::Bool = false)

Use convex hull algorithm to determine the feasibility of each species in `sd`. `sd` should be a `DataFrame`
as generated by `generate_spp_data()`. Returns a vector of indices giving the subset of species identifiers
in `sd` which are feasible given their break-even time and critical water content.
"""
function feas(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04,
              F::Float64 = 10.0, μ::Float64 = 0.1, θ_fc::Float64 = 0.4, uf::Float64 = 0.1,
              analytical::Bool = false, understory_transpiration::Bool = false)

    W₀ = minimum([minimum(sd[:, :Wᵢ]) + ρ, θ_fc])

    ## first test for feasibility in water-only scenario
    infeas = μ_feas(sd, F, μ, analytical)
    first_feas = minimum(findall(sd.Wᵢ .< W₀), init = 9999)
    if first_feas == 9999 || length(infeas) == 0 || maximum(infeas) < first_feas ## cases where no species are feasible
        return []
    else

        ## create dataframe of species characteristics (including w₀) for convex hull algorithm
        points = DataFrame(τ = vcat(sd.τ[infeas], 0.0), Wᵢ = vcat(sd.Wᵢ[infeas], W₀),
                           spp = vcat(sd.spp[infeas], 0))
        ## remove all infeasible species from points dataframe
        points = points[first_feas:nrow(points), :]
        sort!(points, :Wᵢ, rev = true) ## sort in reverse drought tolerance order


        if nrow(points) == 1 ## case where all species infeasible
            return []
        else

            coex = [] ## create empty object to populate with coexisting species indices
            e = 1
            ## loop over species starting with least drought tolerant
            while e < nrow(points)-1

                ## calculate slope between first point and all others
                points.slope .= (points.Wᵢ[e] .- points.Wᵢ) ./ (points.τ[e] .- points.τ)
                cand = points.slope[e+1:nrow(points)]

                ## find the species with the steepest slope (must be negative)
                k = findall(points.slope .== minimum(cand[cand .< 0.0]))[1]
                if k == nrow(points) ## if k is the most drought tolerant species
                    break
                else
                    push!(coex, points[k,:spp]) ## add species to list of coexisting species
                    e = k ## start next iteration from this species
                end
            end
            push!(coex, points[nrow(points),:spp]) ## add latest feasible species to list (automatically coexists)

            ## Now determine whether the species which coexist on water alone close the canopy
            if ρ > calc_ρ_cc(T, sd[coex, :], F, μ, b)

                ## determine best height strategy and eliminate all others
                ht_list = unique(sd[coex, :ht])
                xss_list = Vector{Float64}(undef, length(ht_list))
                sd_sub = sd[coex, :]
                ## loop through height strategies and calculate equilibrium canopy closure size
                for i in 1:length(ht_list)
                    sd_sub_sub = sd_sub[sd_sub.ht .== ht_list[i], :]
                    xss_list[i] = calc_xss(sd_sub_sub[sd_sub_sub.Wᵢ .== maximum(sd_sub_sub.Wᵢ), :C₁][1],
                                           sd_sub_sub[sd_sub_sub.Wᵢ .== maximum(sd_sub_sub.Wᵢ), :C₂][1],
                                           E, T, ρ, F, μ, uf, understory_transpiration)

                end

                if all(xss_list .== 0.0) ## indicates no strategy closes canopy
                else

                    ## identify the best height strategy
                    ht_id = findall(xss_list .^ ht_list .== maximum(xss_list .^ ht_list))

                    ## thin list of coexisting species
                    coex = coex[findall(sd[coex, :ht] .== ht_list[ht_id])]

                end

                ## next check if capped transpiration limits growth period of late species
                if length(coex) == 1 ## if only one species remains in coex list, we are good
                    return coex[sd[coex, :Wᵢ] .< W₀] ## provided that species' Wᵢ is lower than W₀
                ## if expression true, latest species is capped
                elseif ρ > calc_ρ_cx(sd[coex, :], E, T)
                    ## get minimum critical water content from coexisting species
                    Wᵢ_min = sd[coex, :Wᵢ][length(coex)]

                    if length(coex) > 2

                        ## loop through species, eliminating those which are limited by capped transpiration
                        for i in [length(coex)-1:-1:2;]
                            if ρ > calc_ρ_cx(sd[coex[1:i], :], E, T)
                                Wᵢ_min = sd[coex, :Wᵢ][i]
                            else
                                break
                            end
                        end

                    end

                    return coex[sd[coex,:Wᵢ] .< W₀ .&& sd[coex, :Wᵢ] .> Wᵢ_min]
                ## otherwise return current list
                else
                    return coex[sd[coex,:Wᵢ] .< W₀]
                end

            ## if canopy never closes, return list of species which coexist on water alone
            else
                return coex[sd[coex,:Wᵢ] .< W₀]
            end

        end

    end

end;



"""
    calc_gstar(spp_data::DataFrame, ρ::Float64, E::Float64, T::Float64, θ_fc::Float64 = 0.6)

Calculates the equilibrium growth rate (common to all coexisting species) given species' characteristics,
the storm size, `ρ`, transpiration constant `E`, interstorm interval length `T`, and field capacity `θ_fc`.
"""
function calc_gstar(spp_data::DataFrame, ρ::Float64, E::Float64, T::Float64, θ_fc::Float64 = 0.6)
    ((minimum([θ_fc, spp_data[nrow(spp_data), :Wᵢ] + ρ]) - spp_data[1, :Wᵢ]) / (E * T)) *
        (spp_data.C₁[1] + spp_data.C₂[1]) - spp_data.C₂[1]
end



"""
    calc_τ_cc(gstar::Float64, T::Float64, C₁::Float64, C₂::Float64)

Calculates the equilibrium between-storm growth period length (the break-even time, `τ`) when the canopy is closed
as a function of the equilibrium growth rate `gstar`, the interstorm interval length `T`, non-water-limited
growth rate `C₁` and water-limited growth rate `C₂`. This function calculates `τ` for an individual species, rather
than for each species in the community as some other functions do.
"""
function calc_τ_cc(gstar, T, C₁, C₂)
    ((T * (gstar + C₂)) / (C₁ + C₂))
end



"""
    calc_Λ(ρ::Float64, w::Vector{Float64}, τ::Float64, E::Float64, θ_fc::Float64 = 0.6)

Calculates the equilibrium canopy area for each species in a community given the storm size `ρ`,
a vector of species' critical water contents `w`, a vector of species' break-even times `τ`, the transpiration
constant `E`, and the field capactity `θ_fc`.

Note that this function does not do a feasibility screen. The supplied vectors of species characteristics
must be pre-screened.
"""
function calc_Λ(ρ::Float64, w::Vector{Float64}, τ::Vector{Float64}, E::Float64, θ_fc::Float64 = 0.6)
    Λ = Vector{Float64}(undef, length(τ))
    w = vcat(minimum([θ_fc, w[length(w)]+ρ]), w)
    τ = vcat(0.0, τ)

    ## loop through all species except most drought tolerant species
    for i in 2:length(w)-1
        Λ[i-1] = (1/E) * ((w[i-1] - w[i]) / (τ[i] - τ[i-1]) - (w[i] - w[i+1]) / (τ[i+1] - τ[i]))
    end

    ## calculate Λ for most drought tolerant species
    Λ[length(w)-1] = (1 / E) * ((w[length(w)-1] - w[length(w)]) / (τ[length(τ)] - τ[length(τ)-1]))

    return Λ

end



"""
    calc_eqN_cc(spp_data::DataFrame, ρ::Float64, E::Float64,
                T::Float64, F::Float64, μ::Float64, θ_fc::Float64 = 0.6)

Calculates equilibrium abundance of species in `spp_data` provided those species close the canopy.
This function does not itself check whether the species close the canopy, rather it is intended as a
helper function for the more general function `calc_eqN()`.

"""
function calc_eqN_cc(spp_data, ρ, E, T, F, μ, θ_fc::Float64 = 0.6)
    gstar = calc_gstar(spp_data, ρ, E, T, θ_fc)
    τ = Vector{Float64}(undef, nrow(spp_data))
    for i in 1:nrow(spp_data)
        τ[i] = calc_τ_cc(gstar, T, spp_data.C₁[i], spp_data.C₂[1])
    end
    Λ = calc_Λ(ρ, spp_data.Wᵢ, τ, E, θ_fc)
    return F .* Λ .* T .* (exp(μ * T) / (exp(μ * T) - 1))
end



"""
    calc_eqN(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
             θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
             uf::Float64 = 0.1)

Calculate equilibrium population density for species in `spp_data`. Determines whether species close the canopy, and then
calculates resulting equilibria accordingly.

"""
function calc_eqN(sd::DataFrame, T::Float64 = 0.1, ρ::Float64 = 0.04, F::Float64 = F, E::Float64 = E,
                  θ_fc::Float64 = 0.4, μ::Float64 = 0.1, analytical::Bool = false, understory_transpiration::Bool = false,
                  uf::Float64 = 0.1)

    output = Vector{Int64}(undef, nrow(sd))
    eqN = Vector{Float64}(undef, nrow(sd))
    sort!(sd, :Wᵢ, rev = true) ## ensure sd is ordered correctly

    ## check feasibility of each species in sd
    fs = feas(sd, T, ρ, F, μ, θ_fc, uf, analytical, understory_transpiration)

    if isempty(fs) ## if no species are feasible set each pop. density to 0
        sd.eqN .= 0.0
    else

        ## calculate starting soil water content, capped at θ_fc (field capacity)
        W₀ = minimum([θ_fc, sd[fs[length(fs)], :Wᵢ] + ρ])

        ## set equilibrium pop. density of all species not in fs (feasible spp) to 0
        eqN[(!in).(sd.spp, Ref(fs))] .= 0.0

        ## check if canopy closed
        if ρ > calc_ρ_cc(T, sd[fs,:], F, μ, b)

            ## calculate equilibrium population density for closed canopy system
            eqN[fs] .= calc_eqN_cc(sd[fs, :], ρ, E, T, F, μ, θ_fc)

            ## and set to output
            sd.eqN = eqN

        ## if canopy is open
        else

            ## if analytical is set to true, treat system as continuously reproducting
            if analytical

                ## calculate constant that will be reused in subsequent expressions
                cnst = F * T / E

                ## if only one species has a feasible equilibrium, calculate equilibrium reproduction
                if length(fs) == 1
                    eqN[fs[1]] = cnst * ((W₀ - sd[fs[1], :Wᵢ]) / (sd[fs[1], :τ]))
                else
                    ## loop through feasible species and calculate equilibrium reproduction, expression depends on relative phenology
                    for s in [1:1:length(fs);]
                        if s == 1
                            eqN[fs[s]] = cnst * ((W₀ - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ]) -
                                                 (sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ]))
                        elseif s == length(fs)
                            eqN[fs[s]] = cnst * ((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ]))
                        else
                            eqN[fs[s]] = cnst * (((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ])) -
                                                 ((sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ])))
                        end
                    end
                end

            ## if analytical is set to false, treat system as if reproduction occurs once per storm interval
            else

                ## calculate constant that will be reused in subsequent expressions
                cnst = F * T / E

                if length(fs) == 1
                    eqN[fs[1]] = cnst * ((W₀ - sd[fs[1], :Wᵢ]) / (sd[fs[1], :τ]))
                else
                    for s in [1:1:length(fs);]
                        if s == 1
                            eqN[fs[s]] = cnst * ((W₀ - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ]) -
                                (sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ]))
                        elseif s == length(fs)
                            eqN[fs[s]] = cnst * ((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ]))
                        else
                            eqN[fs[s]] = cnst * (((sd[fs[s-1], :Wᵢ] - sd[fs[s], :Wᵢ]) / (sd[fs[s], :τ] - sd[fs[s-1], :τ])) -
                                    ((sd[fs[s], :Wᵢ] - sd[fs[s+1], :Wᵢ]) / (sd[fs[s+1], :τ] - sd[fs[s], :τ])))
                        end
                    end
                end
            end
            ## multiply equilibrium density by following constant to account for mortality
            sd.eqN = eqN .* (exp(μ * T) / (exp(μ * T) - 1))
        end
    end

    ## sort output and return
    output = sort(sd, :spp)
    return output

end;



"""
    calc_eq_biomass(spp_data::DataFrame, ρ::Float64, E::Float64, T::Float64, F::Float64,
                    μ::Float64, uf::Float64, θ_fc::Float64 = 0.6)

Calculates equilibrium biomass of species in `spp_data.`
"""
function calc_eq_biomass(spp_data, ρ, E, T, F, μ, uf, θ_fc::Float64 = 0.6, analytical::Bool = false,
                         understory_transpiration::Bool = false ## no support for this yet
                         )

    ## check feasibility of each species in spp_data
    fs = feas(spp_data, T, ρ, F, μ, θ_fc, uf, analytical, understory_transpiration)

    ## create empty vector of eq. biomass values, to be populated later
    eqB = zeros(nrow(spp_data))

    if isempty(fs) ## if no species are feasible set each pop. density to 0
        return eqB
    else

        ## remove infeasible species from dataframe
        spp_data = spp_data[fs, :]

        ## If canopy is open
        if ρ < calc_ρ_cc(T, spp_data, F, μ, b)

            Λ = calc_Λ(ρ, spp_data.Wᵢ, spp_data.τ, E, θ_fc)

            r = F .* Λ

            function int(t)
                exp(-μ * t) * t^b * (1 / (F * T^b * sigma(μ*T, b-1))) ^ (b / (b-1))
            end

            eqB[fs] .= r .* quadgk(int, 0.0, 1e3)[1]
            return eqB

        else

            ## calculate equilibrium growth rate (common to all species)
            gstar = calc_gstar(spp_data, ρ, E, T, θ_fc)

            ## use gstar to calculate break-even time for each species
            τ = Vector{Float64}(undef, nrow(spp_data))
            for i in 1:nrow(spp_data)
                τ[i] = calc_τ_cc(gstar, T, spp_data.C₁[i], spp_data.C₂[1])
            end

            ## calculate equilibrium canopy area of each species
            Λ = calc_Λ(ρ, spp_data.Wᵢ, τ, E, θ_fc)
            ## calculate initial density of each cohort at equilibrium (aka equilibrium reproduction)
            r = F .* Λ

            ## create new function to solve for canopy closure size
            function cc(xs)

                function intcc(x)
                    exp(-μ * x) * (xs[1] + gstar * x)^2
                end

                [F * exp(-μ * (xs[1] / (uf * gstar))) * quadgk(intcc, 0.0, 1e4)[1] - 1.0]

            end

            ## numerically solve for canopy closure size
            xstar = nlsolve(cc, [0.01]).zero[1]

            ## first integral sums individual biomass values of all cohorts with x < xstar (understory cohorts),
            ## discounted for mortality
            function int1(t)
                exp(-μ * t) * (uf * gstar * t) ^ b
            end

            ## second integral sums individual biomass values of all cohorts in overstory, discounted for
            ## mortality
            function int2(t)
                exp(-μ * t) * (xstar + gstar * (t - (xstar / (uf * gstar)))) ^ b
            end

            ## multiply sum of the individual biomass integrals by the starting population density of each cohort, r
            eqB[fs] .= r .* (quadgk(int1, 0.0, (xstar / (uf * gstar)))[1] +
                         quadgk(int2, (xstar / (uf * gstar)), 1e3)[1])
            return eqB
        end
    end
end



"""
    check_eq(bd, tol = 1e-15)

Checks whether system has reached equilibrium given a DataFrame of density dynamics, as
generated by `sim_water_ppa()` and a tolerance `tol.`
"""
function check_eq(bd, tol = 1e-15)
    check = (Array(bd[nrow(bd),[2:1:ncol(bd);]]) .-
        Array(bd[nrow(bd)-1,[2:1:ncol(bd);]])) ./
        Array(bd[nrow(bd)-1,[2:1:ncol(bd);]])
    if all(check[.!isnan.(check)] .< tol)
        return true
    else return false
    end
end;


"""
    check_eq_agreement(iter::Int64, Nspp::Int64, Nyr::Int64 = 2000,
                            W0::Float64 = 0.6)

Checks agreement between equilibrium population density calculations and simulations.
"""
function check_eq_agreement(iter::Int64, Nspp::Int64, Nyr::Int64 = 300, θ_fc::Float64 = 0.4,
                            P::Float64 = 10, mean_p::Float64 = 0.4,
                            μ::Float64 = 0.15, F::Float64 = 10.0, uf::Float64 = 0.1, n_hts::Int64 = 1)
    eq_sim = Vector{Float64}(undef, iter*Nspp)
    eq_an = Vector{Float64}(undef, iter*Nspp)
    mt = mortality_table(Int(Nyr * P), μ, repeat([1.0 / P], inner = Int(Nyr * P)))
    Threads.@threads for i in 1:iter
        spp_data = generate_spp_data(Nspp, 0.7, n_hts, 1.0 / P, F, μ, 2.5, 0.4, 0.0, 0.0001, 0.00005, 11.0, 0.3, 0.6)
        eq_sim[[((i-1)*Nspp)+1:1:i*Nspp;]] =
            Vector(sim_water_ppa(spp_data, Nyr, Nspp, repeat([1.0], Nspp), μ, F, P, mean_p, θ_fc,
                    mt, false, 0.4, 3.0, uf, false)[2][Int(round(Nyr*P)),2:Nspp+1])
        eq_an[[((i-1)*Nspp)+1:1:i*Nspp;]] =
           calc_eqN(spp_data, 1.0 / P, mean_p / P, F, E, θ_fc, μ, false, false)[:,:eqN]
    end
    return DataFrame(eq_sim = eq_sim, eq_an = eq_an)
end;



"""
    plot_eq_agreement(data::DataFrame, save::Bool = false, filename::String = "")

Generates a plot of calculated vs. simulated equilibrium population density.
"""
function plot_eq_agreement(data::DataFrame, save::Bool = false, filename::String = "")
   p = plot(framestyle = :box, grid = false,
            legend = :none, ylab = "Equilibrium density (simulated)", xlab = "Equilibrium density (predicted)",
            ylim = [0, maximum(data.eq_sim)], xlim = [0, maximum(data.eq_an)])
   p = Plots.abline!(1,0)
   p = plot!(data.eq_an, data.eq_sim, seriestype = "scatter")
   if save
       savefig(filename)
   end
   return p
end
