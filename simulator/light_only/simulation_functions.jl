## define parameters

function sigma(μ, ex)
    Float64(polylog(-ex, exp(-μ)))
end;

## calculate average growth rates
function calc_g(t::Float64, C₁::Float64, C₂::Float64, T::Float64)
    (t * C₁ - (T-t) * C₂)/T
end;

function calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)
    eqN ./ (μ .* F)
end

## grow plant biomass
function grow(biomass::Vector{Float64}, height::Vector{Float64}, zstar::Float64, b::Float64, g::Float64, gᵤ::Float64)

    for i in 1:length(biomass)
        if height[i] > zstar
            biomass[i] = maximum([biomass[i] .^ (1/b) .+ (g), 0.0])^b
        else
            biomass[i] = maximum([biomass[i] .^ (1/b) .+ (gᵤ), 0.0])^b
        end
    end

    return biomass

end;

function mortality_table(Nyr, μ, Tvec)
    mt = zeros(Nyr+1, Nyr+1)
    for i in 1:Nyr
        for j in 1:Nyr
            if j >= i
                mt[i,j] = exp(-μ * sum(Tvec[i:j]))
            end
        end
    end
    mt
end

## kill plants
function die!(nd::DataFrame, rd::DataFrame, yr::Int64, mt::Matrix{Float64})
    for j in 2:ncol(nd)
        nd[:,j] .= rd[:,j] .* mt[:,yr]
    end
    nothing
end;

## aggregate biomass by species
function ΣB(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})
    for i in [2:1:ncol(biomass_data);]
        v[i-1] = sum(biomass_data[:,i] .* n_data[:,i])
    end
    v
end;

## aggregate leaf area by species
function ΣL(biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64})
    for i in 2:ncol(biomass_data)
        v[i-1] = sum(biomass_data[:,i] .^ ((b-1)/b) .* n_data[:,i])
    end
    v
end;

function calc_ht(biomass_data::DataFrame, height_data::DataFrame, spp_data::DataFrame)
    for i in 2:ncol(biomass_data)
            height_data[:,i] .= biomass_data[:,i] .^ (spp_data[i-1,:ht]/b)
    end

    return height_data
end;

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

## aggregate population by species and add to dynamics data
function Σn!(n_data::DataFrame, n_dynamics::DataFrame, yr::Int64)
    for i in [1:1:ncol(n_data)-1;]
            n_dynamics[yr,i+1] = sum(n_data[:,i+1])
    end
    nothing
end;

## aggregate population by species and add to dynamics data
function birth(n_data::DataFrame, height_data::DataFrame, biomass_data::DataFrame, zstar::Float64, v::Vector{Float64},
                yr::Int64, F::Float64)
    v .= 0.0
    for i in 1:nrow(n_data)
        for j in 2:ncol(n_data)
            if height_data[i,j] > zstar
                v[j-1] = v[j-1] + biomass_data[i,j] * n_data[i,j]
            end
        end
    end

    for i in 2:ncol(n_data)
        n_data[yr+1, i] = v[i-1] * F
    end

    return n_data
end;

function generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.07,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.001)

    spp_data = DataFrame(C₁ = rand(Uniform(0.0000005, C₁max), Nspp) .* 365.0,
                         C₂ = repeat([C₂max], Nspp) .* 365.0);
    spp_data.spp = Int[1:1:Nspp;];
    spp_data.ht = rand(Uniform(0.4, 0.6), nrow(spp_data))
    return spp_data
end;

## perform iterations and return output
function iterate_water_ppa(Nyr::Int64, spp_data::DataFrame,
                                biomass_data::DataFrame, biomass_dynamics::DataFrame,
                                n_data::DataFrame, n_dynamics::DataFrame,
                                r_data::DataFrame,
                                height_data::DataFrame, canopy_dynamics::DataFrame,
                                g::Vector{Float64},
                                v::Vector{Float64}, Nspp::Int64,
                                μ::Float64 = 0.1, F::Float64 = 10.0, mt::Matrix{Float64} = zeros(1,1),
                                understory_factor::Float64 = 0.1, pb::Bool = true)

    ## begin iterating
    if pb
        prog = ProgressBar(total = Nyr)
    end

    zstar = 0.0

    for yr in 1:Nyr

        ## calculate canopy closure height
        zstar = calc_zstar(biomass_data, height_data, n_data)

        for s in spp_data.spp
            biomass_data[1:yr, s+1] =
                grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, spp_data.C₁[s],
                     spp_data.C₁[s]*understory_factor)
        end

        die!(n_data, r_data, yr, mt)

        height_data = calc_ht(biomass_data, height_data, spp_data)

        v = ΣB(biomass_data, n_data, v)
        biomass_dynamics[yr, [2:1:Nspp+1;]] = v

        canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

        ## generate new cohort and record population dynamics
        n_data = birth(n_data, biomass_data, height_data, zstar, v, yr, F)
        for i in 2:ncol(r_data)
            r_data[yr+1, i] = n_data[yr+1, i]
        end
        Σn!(n_data, n_dynamics, yr)

        ## extinction cutoff
        r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

        if pb
            update(prog)
        end

    end

    return [biomass_data, n_dynamics, n_data, r_data, height_data, canopy_dynamics, zstar]

end;


## simulate water only
function sim_water_ppa(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                       Ninit::Float64, μ::Float64 = 0.15, F::Float64 = 10.0,
                       mt::Matrix{Float64} = zeros(1,1),
                       pb::Bool = true)

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
    t = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, repeat([1.0], inner = Nyr))
    end

    iterate_water_ppa(Nyr, spp_data,
                      biomass_data, biomass_dynamics,
                      n_data, n_dynamics,
                      r_data,
                      height_data, canopy_dynamics,
                      g, v, Nspp, μ, F, mt,
                      0.05,
                      pb)

end;

function mono_zstar(spp_data, F, P, μ, Nyr = 200)

    zs = Vector{Float64}(undef, nrow(spp_data))
    for i in 1:nrow(spp_data)
        sd = DataFrame(spp_data[i,:])
        sd[:,:spp] = [1]
        out = sim_water_ppa(sd, 200, 1, Ninit, μ, F)
        zs[i] = out[7]
    end

    return zs

end

function water_only(spp_data, Nyr, Ninit, μ, F, P, mean_p, θ_fc)

    include("../water_only/simulation_functions.jl")
    sim_water_only(spp_data, Nyr, nrow(spp_data), Ninit, μ, F, P, mean_p, θ_fc)

end
