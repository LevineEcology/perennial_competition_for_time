using DataFrames: DataFrameColumns, PrettyTablesConf
## define parameters


"""
    canopy_proportion(height_data::DataFrame, n_data::DataFrame,
                      biomass_data::DataFrame, zstar::Float64)

Helper function for running simulations.

Calculates the proportion of total canopy area occupied by each species in a community,
given `height_data`, `n_data`, `biomass_data`, and the canopy closure height, `zstar`.

"""
function canopy_proportion(height_data::DataFrame, n_data::DataFrame, biomass_data::DataFrame, zstar::Float64)
    ht_vec = vec(reshape(Matrix(height_data[:,2:ncol(height_data)]), 1, :))
    n_vec = vec(reshape(Matrix(n_data[:,2:ncol(n_data)]), 1, :))

    t_ind = findall(ht_vec .> zstar)
    total_canopy = sum(n_vec[t_ind])

    canopy_p = Vector{Float64}(undef, ncol(height_data)-1)
    for i in 2:ncol(height_data)
        ind = findall(height_data[:, i] .> zstar)
        canopy_p[i-1] = sum(n_data[ind, i] .* biomass_data[ind, i] .^ ((b-1)/b))
    end
    return canopy_p

end



"""
    calc_t!(t::Vector{Float64}, biomass_data::DataFrame,
            n_data::DataFrame, v::Vector{Float64},
            Wᵢ::Vector{Float64}, W₀::Float64,
            T::Float64, ht_data::DataFrame, zstar::Float64)

Helper function for running simulations.

Calculates the growing period length of each species, modifying the vector
`t`, using information on biomass `biomass_data`, population density `n_data`,
and the species' characteristics.
"""
function calc_t!(t::Vector{Float64}, biomass_data::DataFrame,
                 n_data::DataFrame, v::Vector{Float64},
                 Wᵢ::Vector{Float64}, W₀::Float64,
                 T::Float64, ht_data::DataFrame, zstar::Float64)

    ## determine the proportion of the canopy occupied by each species
    x = canopy_proportion(ht_data, n_data, biomass_data, zstar)

    ## only perform calculations for spp with Wᵢ < W₀
    g0 = findall(W₀ .> Wᵢ)
    t[1:g0[1]-1] .= 0.0
    t[g0] .= T
    if length(g0) != 0
        for s in g0[1]:length(x)
            if s == g0[1] ## special case for earliest species
                t[s] = (W₀ - Wᵢ[s]) / (E * sum(x[g0]))
            else
                t[s] = ((Wᵢ[s-1] - Wᵢ[s]) / (E * sum(x[s:length(v)]))) + t[s-1]
            end
            if t[s] > T
                t[s] = T
                break
            end
        end
    end
    nothing
end

"""
    calc_w(w_cur::Float64, t::Vector{Float64}, T::Float64,
           biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64},
           Wᵢ::Vector{Float64}, ht_data::DataFrame, zstar::Float64)

Helper function for running simulations.

Calculates the amount of soil water left after all species have ceased growth.
`w_cur` is the starting soil volumetric water content.
"""
function calc_w(w_cur::Float64, t::Vector{Float64}, T::Float64,
                biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64},
                Wᵢ::Vector{Float64}, ht_data::DataFrame, zstar::Float64)

    ## determine the proportion of the canopy occupied by each species
    v = canopy_proportion(ht_data, n_data, biomass_data, zstar)
    ## determine which species have critical water contents below starting water content
    sb = findall(Wᵢ .< w_cur)
    if isempty(sb) ## if no species grow, return current soil water content
        return w_cur
    elseif all(t[sb] .< T) ## if all species grow, return critical water content of latest species
        return Wᵢ[length(t)]
    elseif all(t[sb] .== T) ## if all species grow, but do not shut off
         return w_cur - E * sum(v[sb]) * T
    else ## if only some species grow
        fs = minimum(findall(t .== T))
        if fs == sb[1]
            return w_cur - E * sum(v[fs:length(t)]) * (T - t[fs-1])
        else
            return Wᵢ[fs-1] - E * sum(v[fs:length(t)]) * (T - t[fs-1])
        end
    end
end



"""
    grow(biomass::Vector{Float64}, height::Vector{Float64}, zstar::Float64,
         b::Float64, g::Float64, gᵤ::Float64, T::Float64)

Helper function for simulations.

Increments biomass of all species in a simulation given their growth rate in the overstory,
`g` and understory `gᵤ.`
"""
function grow(biomass::Vector{Float64}, height::Vector{Float64}, zstar::Float64, b::Float64, g::Float64, gᵤ::Float64, T::Float64)

    for i in 1:length(biomass)
        if height[i] >= zstar
            biomass[i] = maximum([biomass[i] .^ (1/b) .+ (g * T), 0.0])^b
        else
            biomass[i] = maximum([biomass[i] .^ (1/b) .+ (gᵤ * T), 0.0])^b
        end
    end

    return biomass

end;

"""
    mortality_table(Nyr::Int64, μ::Float64, Tvec::Vector{Float64})

Helper function for simulations.

Generates a table of mortality discounting constants to facilitate rapid simulations. Has a
large upfront computational cost, but saves time overall for long simulation runs or large communities.
"""
function mortality_table(Nyr::Int64, μ::Float64, Tvec::Vector{Float64})
    ## generate empty array
    mt = zeros(Nyr+1, Nyr+1)

    ## run computations in parallel
    Threads.@threads for i in 1:Nyr
        for j in 1:Nyr
            if j >= i
                mt[i,j] = exp(-μ * sum(Tvec[i:j]))
            end
        end
    end
    mt
end

"""
    die!(nd::DataFrame, rd::DataFrame, yr::Int64, mt::Matrix{Float64})

Helper function for simulations.

Reduces population densities for cohorts in simulation given mortality table `mt.`
"""
function die!(nd::DataFrame, rd::DataFrame, yr::Int64, mt::Matrix{Float64})
    for j in 2:ncol(nd)
        nd[:,j] .= rd[:,j] .* mt[:,yr]
    end
    nothing
end;

"""
    ΣB(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})

Helper function for simulations.

Calculates the total biomass of each species in a simulation.
"""
function ΣB(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})
    for i in [2:1:ncol(biomass_data);]
        v[i-1] = sum(biomass_data[:,i] .* n_data[:,i])
    end
    v
end;

"""
    ΣL(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})

Helper function for simulations.

Calculates the total canopy area of each species in a simulation.
"""
function ΣL(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})
    for i in 2:ncol(biomass_data)
        v[i-1] = sum(biomass_data[:,i] .^ ((b-1)/b) .* n_data[:,i])
    end
    v
end;


"""
    calc_ht(biomass_data::DataFrame, height_data::DataFrame, spp_data::DataFrame)

Helper function for simulations.

Calculates the height of each cohort of each species in a simulation.
"""
function calc_ht(biomass_data::DataFrame, height_data::DataFrame, spp_data::DataFrame)
    for i in 2:ncol(biomass_data)
        height_data[:,i] .= biomass_data[:,i] .^ (spp_data[i-1,:ht]/b)
    end

    return height_data
end;



"""
    calc_zstar(biomass_data::DataFrame, height_data::DataFrame, n_data::DataFrame)

Helper function for simulations.

Determines the canopy closure height for a community given data on their biomass, height and population
density.
"""
function calc_zstar(biomass_data::DataFrame, height_data::DataFrame, n_data::DataFrame)


    ht_totals = Vector{Float64}(undef, nrow(biomass_data)*(ncol(biomass_data)-1))

    ht_vec = vec(reshape(Matrix(height_data[:,2:ncol(height_data)]), 1, :))
    biom_vec = vec(reshape(Matrix(biomass_data[:,2:ncol(biomass_data)]), 1, :))
    n_vec = vec(reshape(Matrix(n_data[:,2:ncol(n_data)]), 1, :))

    ord = sortperm(ht_vec, rev = true)

    zstarind = 0

    ht_totals[1] = biom_vec[ord[1]]^((b-1)/b) * n_vec[ord[1]]
    for i in 2:length(ht_vec)
        ht_totals[i] = ht_totals[i-1] + biom_vec[ord[i]]^((b-1)/b) * n_vec[ord[i]]
        if ht_totals[i] >= 1.0
            zstarind = i
            break
        end
    end

    if zstarind == 0
        return 0.0
    else
        return ht_vec[ord[zstarind]]
    end

end;

"""
    Σn!(n_data::DataFrame, n_dynamics::DataFrame, yr::Int64)

Helper function for simulations.

Sums population densities by species.
"""
function Σn!(n_data::DataFrame, n_dynamics::DataFrame, yr::Int64)
    for i in [1:1:ncol(n_data)-1;]
            n_dynamics[yr,i+1] = sum(n_data[:,i+1])
    end
    nothing
end;


"""
    birth(n_data::DataFrame, canopy_dynamics::DataFrame,
          yr::Int64, F::Float64, T::Float64)

Helper function for simulations.

Creates new cohort for each species based on current canopy occupancy. Assumes fecundity proportional
to canopy occupancy.
"""
function birth(n_data::DataFrame, canopy_dynamics::DataFrame,
                yr::Int64, F::Float64, T::Float64)

    ## loop through species and generate new cohorts.
    for i in 2:ncol(n_data)
        n_data[yr+1,i] = canopy_dynamics[yr,i] * F * T
    end

    return n_data
end;


"""
    birth_b(n_data::DataFrame, height_data::DataFrame, biomass_data::DataFrame,
            zstar::Float64, v::Vector{Float64},
            yr::Int64, F::Float64, T::Float64)

Helper function for simulations.

Creates new cohort for each species based on total biomass. Assumes fecundity proportional
to biomass rather than canopy occupancy.

"""
function birth_b(n_data::DataFrame, height_data::DataFrame, biomass_data::DataFrame,
                 zstar::Float64, v::Vector{Float64},
                 yr::Int64, F::Float64, T::Float64)
    v .= 0.0
    for i in 1:nrow(n_data)
        for j in 2:ncol(n_data)
            if height_data[i,j] > zstar
                v[j-1] = v[j-1] + biomass_data[i,j] * n_data[i,j]
            end
        end
    end

    for i in 2:ncol(n_data)
        n_data[yr+1, i] = v[i-1] * F * T
    end

    return n_data
end;

######## SIMULATION FUNCTIONS ########

"""
    iterate_water_ppa!(yr::Int64, biomass_data::DataFrame, biomass_dynamics::DataFrame,
                       n_data::DataFrame, n_dynamics::DataFrame, r_data::DataFrame,
                       w_data::DataFrame, height_data::DataFrame, canopy_dynamics::DataFrame,
                       g::Vector{Float64}, t::Vector{Float64}, v::Vector{Float64}, Nspp::Int64,
                       μ::Float64, F::Float64, mt::Matrix{Float64}, W₀vec::Vector{Float64},
                       Tvec::Vector{Float64}, understory_factor::Float64, θ_fc::Float64, b::Float64)

Iterates through years of a simulation and returns results. Called by simulation wrappers.
"""
function iterate_water_ppa(Nyr::Int64, spp_data::DataFrame,
                           biomass_data::DataFrame, biomass_dynamics::DataFrame,
                           n_data::DataFrame, n_dynamics::DataFrame,
                           r_data::DataFrame, w_data::Vector{Float64},
                           height_data::DataFrame, canopy_dynamics::DataFrame,
                           g::Vector{Float64}, t::Vector{Float64},
                           v::Vector{Float64}, Nspp::Int64,
                           μ::Float64 = 0.1, F::Float64 = 10.0, mt::Matrix{Float64} = zeros(1,1),
                           W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                           Tvec::Vector{Float64} = repeat([40.0], inner = Nyr),
                           understory_factor::Float64 = 0.1,
                           θ_fc::Float64 = 0.4, pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5,
                           perturb::Bool = false, perturb_iter::Int64 = 200, perturb_fac::Float64 = 0.9)

    ## initialize soil water content
    w = minimum([θ_fc, w_init])
    w_in_data = copy(w_data)
    zstar_data = copy(w_data)

    ## begin iterating
    if pb
        prog = ProgressBar(total = Nyr)
    end

    zstar = 0.0

    if perturb

        for yr in 1:perturb_iter-1

            #println(zstar)

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            ## calculate initial water content
            w = minimum([w + W₀vec[yr], θ_fc])
            w_in_data[yr] = w

            ## calculate growth time for each species
            if any(spp_data.Wᵢ .< w)
                calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ,
                        w, Tvec[yr], height_data, zstar)
            else
                t = repeat([0.0], inner = Nspp)
            end

            ## recalculate soil water content
            w = maximum([calc_w(w, t, Tvec[yr], biomass_data, n_data, v, spp_data.Wᵢ, height_data, zstar), 0.0])
            w_data[yr] = w

            ## calculate average growth rate for each species
            g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, g[s],
                         g[s] * understory_factor, Tvec[yr])

            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, Tvec[yr])
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

        r_data[:,2:ncol(r_data)] .= r_data[:,2:ncol(r_data)] .* perturb_fac
        n_data[:,2:ncol(r_data)] .= n_data[:,2:ncol(r_data)] .* perturb_fac

        for yr in perturb_iter:Nyr

            #println(zstar)

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            ## calculate initial water content
            w = minimum([w + W₀vec[yr], θ_fc])
            w_in_data[yr] = w

            ## calculate growth time for each species
            if any(spp_data.Wᵢ .< w)
                calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ,
                        w, Tvec[yr], height_data, zstar)
            else
                t = repeat([0.0], inner = Nspp)
            end

            ## recalculate soil water content
            w = maximum([calc_w(w, t, Tvec[yr], biomass_data, n_data, v, spp_data.Wᵢ, height_data, zstar), 0.0])
            w_data[yr] = w

            ## calculate average growth rate for each species
            g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, g[s],
                         g[s] * understory_factor, Tvec[yr])

            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, Tvec[yr])
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

    else

        for yr in 1:Nyr

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            ## calculate initial water content
            w = minimum([w + W₀vec[yr], θ_fc])
            w_in_data[yr] = w

            ## calculate growth time for each species
            if any(spp_data.Wᵢ .< w)
                calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ,
                        w, Tvec[yr], height_data, zstar)
            else
                t = repeat([0.0], inner = Nspp)
            end

            ## recalculate soil water content
            w = maximum([calc_w(w, t, Tvec[yr], biomass_data, n_data, v, spp_data.Wᵢ, height_data, zstar), 0.0])
            w_data[yr] = w

            ## calculate average growth rate for each species
            g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, g[s],
                         g[s] * understory_factor, Tvec[yr])

            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, Tvec[yr])
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

    end

    return [biomass_data, n_dynamics, n_data, r_data, w_data,
            w_in_data, height_data, canopy_dynamics, zstar_data, g, t, biomass_dynamics]

end;


"""
    sim_water_ppa(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                       Ninit::Any, μ::Float64 = 0.15, F::Float64 = 10.0,
                       P::Float64 = 40.0, mean_p::Float64 = 16.0, θ_fc = 0.4, mt::Matrix{Float64} = zeros(1,1),
                       pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5, understory_factor::Float64 = 0.1,
                       perturb::Bool = false, perturb_frac::Float64 = 0.5, perturb_factor::Float64 = 0.9,
                       perturb_water::Bool = false, perturb_factor_water::Float64 = 0.7,
                       perturb_water_return::Bool = false, perturb_water_return_frac::Float64 = 0.7)

Perform simulation of a community of species limited by both water and (potentially) light. The community
of species must be described by `spp_data`, a data frame generated by the `generate_spp_data()` function.
"""
function sim_water_ppa(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                       Ninit::Any, μ::Float64 = 0.15, F::Float64 = 10.0,
                       P::Float64 = 40.0, mean_p::Float64 = 16.0, θ_fc = 0.4, mt::Matrix{Float64} = zeros(1,1),
                       pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5, understory_factor::Float64 = 0.1,
                       perturb::Bool = false, perturb_frac::Float64 = 0.5, perturb_factor::Float64 = 0.9,
                       perturb_water::Bool = false, perturb_factor_water::Float64 = 0.7,
                       perturb_water_return::Bool = false, perturb_water_return_frac::Float64 = 0.7)

    Nyr = Nyr * Int(round(P))

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    r_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    height_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                            spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                            N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))

    biomass_data = unstack(biomass_data, :spp, :B)
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]])
    n_data = unstack(n_data, :spp, :N)
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]])
    r_data = unstack(r_data, :spp, :N)
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]])
    height_data = unstack(height_data, :spp, :N)
    height_data[!,[2:ncol(height_data);]] .= convert.(Float64,height_data[!,[2:ncol(height_data);]])

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[1:Nyr,:])
    n_dynamics = copy(n_data[1:Nyr,:])
    canopy_dynamics = copy(n_data[1:Nyr,:])

    ## inital population:
    if typeof(Ninit) == Float64
        n_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
        r_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
    elseif typeof(Ninit) == Vector{Float64}
        n_data[1, 2:Nspp+1] = Ninit
        r_data[1, 2:Nspp+1] = Ninit
    else
        println("Please supply Ninit value that is a Float64 or Vector{Float64}")
    end

    g = Vector{Float64}(undef, Nspp)
    t = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    w_data = Vector{Float64}(undef, Nyr)

    T = 1.0 / P
    W₀ = mean_p / P

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, repeat([T], inner = Nyr))
    end

    ρ_list = repeat([W₀], inner = Nyr)
    if perturb_water
        ρ_list[Int(round(perturb_frac * Nyr)):Nyr] .= perturb_factor_water * W₀
        if perturb_water_return
            ρ_list[Int(round(perturb_water_return_frac * Nyr)):Nyr] .= W₀
        end
    end

    iterate_water_ppa(Nyr, spp_data,
                      biomass_data, biomass_dynamics,
                      n_data, n_dynamics,
                      r_data, w_data,
                      height_data, canopy_dynamics,
                      g, t, v, Nspp, μ, F, mt,
                      ρ_list,
                      repeat([T], inner = Nyr), understory_factor, θ_fc,
                      pb, w_init, b, perturb, Int(round(perturb_frac * Nyr)),
                      perturb_factor)

end;


"""
    mono_zstar(spp_data::DataFrame, F::Float64, P::Float64, μ::Float64,
               mean_p::Float64 = 4.0, θ_fc::Float64 = 0.4, Nyr::Int64 = 200)

Uses simulation approach to determine canopy closure height in equilibrium, `zstar`,
for each species in `spp_data`. Does so by looping through each species and simulating
a monoculture to equilibrium. Computationally expensive and slow for large assemblages.
"""
function mono_zstar(spp_data::DataFrame, F::Float64, P::Float64, μ::Float64,
                    mean_p::Float64 = 4.0, θ_fc::Float64 = 0.4, Nyr::Int64 = 200)

    zs = Vector{Float64}(undef, nrow(spp_data))
    for i in 1:nrow(spp_data)
        sd = DataFrame(spp_data[i,:])
        sd[:,:spp] = [1]
        out = sim_water_ppa(sd, Nyr, 1, [Ninit], μ, F, P, mean_p, θ_fc,
                            zeros(1,1), true, 0.4, 3.0, uf, false)
        zs[i] = out[9][Nyr]
        println(out[2][Nyr,:])
    end

    return zs

end;


"""
    water_only(spp_data::DataFrame, Nyr::Int64, Ninit, μ::Float64,
                    F::Float64, P::Float64, mean_p::Float64, θ_fc::Float64)

Simulates community of species competing for water only, calls simulation function defined
in separate script.
"""
function water_only(spp_data::DataFrame, Nyr::Int64, Ninit, μ::Float64,
                    F::Float64, P::Float64, mean_p::Float64, θ_fc::Float64)

    Include("../water_only/simulation_functions.jl")
    sim_water_only(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, mean_p, θ_fc)

end


##---------------------------------------------------------------
## PPA ONLY
##---------------------------------------------------------------


"""
    iterate_ppa(Nyr::Int64, spp_data::DataFrame,
                biomass_data::DataFrame, biomass_dynamics::DataFrame,
                n_data::DataFrame, n_dynamics::DataFrame,
                r_data::DataFrame, zstar_data::Vector{Float64},
                height_data::DataFrame, canopy_dynamics::DataFrame,
                g::Vector{Float64},
                v::Vector{Float64}, Nspp::Int64,
                μ::Float64 = 0.1, F::Float64 = 10.0, T::Float64 = 0.1,
                mt::Matrix{Float64} = zeros(1,1),
                understory_factor::Float64 = 0.1, pb::Bool = true,
                b::Float64 = 2.5,
                perturb::Bool = false, perturb_iter::Int64 = 200, perturb_fac::Float64 = 0.9)


"""
function iterate_ppa(Nyr::Int64, spp_data::DataFrame,
                     biomass_data::DataFrame, biomass_dynamics::DataFrame,
                     n_data::DataFrame, n_dynamics::DataFrame,
                     r_data::DataFrame, zstar_data::Vector{Float64},
                     height_data::DataFrame, canopy_dynamics::DataFrame,
                     g::Vector{Float64},
                     v::Vector{Float64}, Nspp::Int64,
                     μ::Float64 = 0.1, F::Float64 = 10.0, T::Float64 = 0.1,
                     mt::Matrix{Float64} = zeros(1,1),
                     understory_factor::Float64 = 0.1, pb::Bool = true,
                     b::Float64 = 2.5,
                     perturb::Bool = false, perturb_iter::Int64 = 200, perturb_fac::Float64 = 0.9)

    ## begin iterating
    if pb
        prog = ProgressBar(total = Nyr)
    end

    zstar = 0.0

    if perturb

        for yr in 1:perturb_iter-1

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, spp_data[s,:C₁],
                         spp_data[s,:C₁] * understory_factor, T)
            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, T)
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

        r_data[:,2:ncol(r_data)] .= r_data[:,2:ncol(r_data)] .* perturb_fac
        n_data[:,2:ncol(r_data)] .= n_data[:,2:ncol(r_data)] .* perturb_fac

        for yr in perturb_iter:Nyr

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, spp_data[s,:C₁],
                         spp_data[s,:C₁] * understory_factor, T)
            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, T)
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

    else

        for yr in 1:Nyr

            ## calculate canopy closure height
            zstar = calc_zstar(biomass_data, height_data, n_data)
            zstar_data[yr] = zstar

            for s in spp_data.spp
                biomass_data[1:yr, s+1] =
                    grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, spp_data[s,:C₁],
                         spp_data[s,:C₁] * understory_factor, T)
            end

            die!(n_data, r_data, yr, mt)
            height_data = calc_ht(biomass_data, height_data, spp_data)

            v = ΣB(biomass_data, n_data, v)
            biomass_dynamics[yr, [2:1:Nspp+1;]] = v

            canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

            ## generate new cohort and record population dynamics
            n_data = birth(n_data, canopy_dynamics, yr, F, T)
            for i in 2:ncol(r_data)
                r_data[yr+1, i] = n_data[yr+1, i]
            end
            Σn!(n_data, n_dynamics, yr)

            ## extinction cutoff
            r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

            if pb
                update(prog)
            end

            #println(zstar)

        end

    end

    return [biomass_data, n_dynamics, n_data, r_data,
            height_data, canopy_dynamics, zstar_data]

end;

## simulate water only
function sim_ppa(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                 Ninit::Float64, μ::Float64 = 0.15, F::Float64 = 10.0, T::Float64 = 0.1,
                 mt::Matrix{Float64} = zeros(1,1),
                 pb::Bool = true, b::Float64 = 2.5, understory_factor::Float64 = 0.1,
                 perturb::Bool = false, perturb_frac::Float64 = 0.5, perturb_factor::Float64 = 0.9)

    Nyr = Int(round(Nyr * 1.0/T))

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    r_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    height_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                            spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                            N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))

    biomass_data = unstack(biomass_data, :spp, :B)
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]])
    n_data = unstack(n_data, :spp, :N)
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]])
    r_data = unstack(r_data, :spp, :N)
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]])
    height_data = unstack(height_data, :spp, :N)
    height_data[!,[2:ncol(height_data);]] .= convert.(Float64,height_data[!,[2:ncol(height_data);]])

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[1:Nyr,:])
    n_dynamics = copy(n_data[1:Nyr,:])
    canopy_dynamics = copy(n_data[1:Nyr,:])

    ## inital population:
    n_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
    r_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)

    g = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    zstar_data = Vector{Float64}(undef, Nyr)

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, repeat([T], inner = Nyr))
    end

    iterate_ppa(Nyr, spp_data,
                biomass_data, biomass_dynamics,
                n_data, n_dynamics,
                r_data, zstar_data,
                height_data, canopy_dynamics,
                g, v, Nspp, μ, F, T, mt, understory_factor,
                pb, b, perturb, Int(round(perturb_frac * Nyr)),
                perturb_factor)


end;


##---------------------------------------------------------------
## STOCHASTIC SIMULATIONS
##---------------------------------------------------------------


####### HELPER FUNCTIONS ########

"""
    generate_intervals(Pmean::Float64, cluster::Bool = false)

Helper function for simulations.

Creates a vector of stochastic interstorm interval lengths. Lengths are
drawn from an exponential distribution given mean frequency `Pmean`. If `cluster`
is `true`, cluster all rainfall events at start of the season by sorting by length
in increasing order.
"""
function generate_intervals(Pmean::Int64, cluster::Bool = false)

    ## draw interstorm intervals from exponential distribution with λ = Pmean
    inter = rand(Exponential(1/Pmean), Pmean)

    ## if the total length of the interstorm intervals exceeds 1 year, truncate
    if sum(inter) > 1.0
        rs = 0.0
        for i in 1:length(inter)
            rs = rs + inter[i]
            if rs > 1.0
                inter = inter[1:i]
                inter[i] = 1.0 - sum(inter[1:i-1])
                inter = inter[1:i]
                break
            end
        end
    end
    if cluster
        sort!(vec(inter))
    end
    return inter
end


"""
    generate_rainfall_regime(Nyr::Int64, Pmean::Float64, Pdisp::Float64,
                             map_mean::Float64, map_var::Float64, cluster::Bool = false)

Creates a complete, stochastic rainfall regime for use in stochastic simulations as implemented
in `sim_water_ppa_stochastic()`.

"""
function generate_rainfall_regime(Nyr::Int64, Pmean::Float64, Pdisp::Float64,
                                  map_mean::Float64, map_sd::Float64,
                                  constant_P::Bool = false, cluster::Bool = false)

    if constant_P
        Preal = Int.(repeat([Pmean], Nyr))
        Tlist = repeat([1.0 / Pmean], Nyr * Int.(Pmean))
    else

        ## reparameterize for draws from negative binomial distribution
        var = Pmean + 1 / Pdisp * Pmean ^ 2
        p = (var - Pmean) / var

        ## draw storm frequencies from negative binomial distribution
        Plist = rand(NegativeBinomial(Pdisp, 1-p), Nyr)
        Plist[Plist .== 0] .= 1.0 ## don't allow 0 values

        ## true values of P are determined after calling generate_intervals()
        Preal = Float64[]
        Tlist = Int64[]
        for yr in 1:Nyr
            new_int = generate_intervals(Plist[yr]) ## get random interstorm interval lengths
            Tlist = vcat(Tlist, new_int) ## add to full list of interval lengths
            Preal = vcat(Preal, length(new_int)) ## record true P value as length of new_int
        end
        Preal = Int.(Preal) ## convert to integer

    end

    ## draw precipitation totals from (truncated) normal distribution
    precip_list = rand(Normal(map_mean, map_sd), Nyr)
    precip_list[precip_list .< 0.0] .= 0.0

    ## now calculate storm sizes from yearly precip totals and interval lengths
    ss_list = Vector{Float64}(undef, length(Tlist))
    i = 1
    for yr in 1:Nyr
        ss_list[i:i+Preal[yr]-1] .= Tlist[i:i+Preal[yr]-1] .* precip_list[yr]
        i = i+Preal[yr]
    end
    ss_list[ss_list .< 0.0] .= 0.0 ## no negative storm totals

    return [ss_list, Tlist, Preal]

end

####### SIMULATION FUNCTIONS ########


"""
    sim_water_ppa_stochastic(spp_data::DataFrame, Nyr::Int64, Nspp::Int64, Ninit::Any,
                             rainfall_regime::Vector{Vector{Float64}},
                             F::Float64 = 10.0, θ_fc::Float64 = 0.4, μ::Float64 = 0.05,
                             mt::Matrix{Float64} = zeros(1,1), w_init = 0.4, b::Float64 = 2.5,
                             understory_factor::Float64 = 0.1, pb::Bool = true)

Simulate community dynamics under stochastic rainfall regime.

"""
function sim_water_ppa_stochastic(spp_data::DataFrame, Nyr::Int64, Nspp::Int64, Ninit::Any,
                                  rainfall_regime::Vector{Vector{Float64}},
                                  F::Float64 = 10.0, θ_fc::Float64 = 0.4, μ::Float64 = 0.05,
                                  mt::Matrix{Float64} = zeros(1,1), w_init = 0.4, b::Float64 = 2.5,
                                  understory_factor::Float64 = 0.1, pb::Bool = true)

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    r_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    height_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                            spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                            N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))

    biomass_data = unstack(biomass_data, :spp, :B)
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]])
    n_data = unstack(n_data, :spp, :N)
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]])
    r_data = unstack(r_data, :spp, :N)
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]])
    height_data = unstack(height_data, :spp, :N)
    height_data[!,[2:ncol(height_data);]] .= convert.(Float64,height_data[!,[2:ncol(height_data);]])

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[1:Nyr,:])
    n_dynamics = copy(n_data[1:Nyr,:])
    canopy_dynamics = copy(n_data[1:Nyr,:])

    ## inital population:
    if typeof(Ninit) == Float64
        n_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
        r_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
    elseif typeof(Ninit) == Vector{Float64}
        n_data[1, 2:Nspp+1] = Ninit
        r_data[1, 2:Nspp+1] = Ninit
    else
        println("Please supply Ninit value that is a Float64 or Vector{Float64}")
    end

    g = Vector{Float64}(undef, Nspp)
    t = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    w_data = Vector{Float64}(undef, Nyr)

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, rainfall_regime[2])
    end

    iterate_water_ppa(Nyr, spp_data,
                      biomass_data, biomass_dynamics,
                      n_data, n_dynamics,
                      r_data, w_data,
                      height_data, canopy_dynamics,
                      g, t, v, Nspp, μ, F, mt,
                      rainfall_regime[1],
                      rainfall_regime[2], understory_factor, θ_fc,
                      pb, w_init, b, false)

end;



## simulate water only
function sim_water_only_stochastic(spp_data::DataFrame, Nspp::Int64,
                                   Ninit::Float64, rainfall_regime::Vector{Vector{Float64}},
                                   θ_fc::Float64 = 0.4, μ::Float64 = 0.05, F::Float64 = 10.0,
                                   mt::Matrix{Float64} = zeros(1,1), pb::Bool = true)

    Nyr = length(rainfall_regime[1])

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
    r_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));

    biomass_data = unstack(biomass_data, :spp, :B);
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]]);
    n_data = unstack(n_data, :spp, :N);
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]]);
    r_data = unstack(r_data, :spp, :N);
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]]);

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
    n_dynamics = copy(n_data[[1:1:Nyr;],:]);

    ## inital population:
    n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);
    r_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

    g = Vector{Float64}(undef, Nspp);
    t = Vector{Float64}(undef, Nspp);
    v = Vector{Float64}(undef, Nspp);

    w_data = Vector{Float64}(undef, Nyr);

    W₀vec = copy(rainfall_regime[1])
    W₀vec[W₀vec .> θ_fc] .= θ_fc

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, rainfall_regime[2])
    end

    iterate_water_only_sim(Nyr, spp_data,
                           biomass_data, biomass_dynamics,
                           n_data, n_dynamics, r_data, w_data,
                           g, t, v, Nspp, mt, μ, F,
                           vec(W₀vec), vec(rainfall_regime[2]), θ_fc,
                           pb)

end;

