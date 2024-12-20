## define parameters
##calculate growth time within season
function calc_t!(t::Vector{Float64}, biomass_data::DataFrame,
                 n_data::DataFrame, v::Vector{Float64},
                 Wᵢ::Vector{Float64}, W₀::Float64 = 0.6, T::Float64 = 40.0)
    v = ΣL(biomass_data, n_data, v)

    ## only perform calculations for spp with Wᵢ < W₀
    g0 = findall(W₀ .> Wᵢ)
    t[1:g0[1]-1] .= 0.0
    t[g0] .= T
    if length(g0) == 0
        g0 = [1]
    end

    for s in g0[1]:length(v)
        if s == 1
            t[s] = (W₀ - Wᵢ[s]) / (E * sum(v))
        else
            t[s] = ((Wᵢ[s-1] - Wᵢ[s]) / (E * sum(v[[s:1:length(v);]]))) + t[s-1]
        end
        if t[s] > T
            t[s] = T
            break
        end
    end
    nothing
end

## calculate average growth rates
function calc_g(t::Float64, C₁::Float64, C₂::Float64, T::Float64 = 40.0)
    (t * (C₁ + C₂) - T * C₂)/T
end;

## doesnt give exact result, though it is correct to at least 3 decimal places
function calc_τ(C₁::Float64, C₂::Float64, F::Float64, mu::Float64, T::Float64, expon::Float64)
    x = Vector{Int64}([0:1:100000;])
    y = (1-mu).^x .* x.^expon
    (1/(C₁ + C₂)) * (((1/F) * 1/(sum(y)))^(1/expon) + T*C₂)
end

## grow plant biomass
function grow(biomass::Vector{Float64}, b::Float64, g::Float64, T::Float64)
    replace(x -> isless(x, 0) ? 0 : x, (biomass .^ (1/b) .+ (g .* T))) .^ b
end;

## kill plants
function die!(nd::DataFrame, μ::Float64, yr::Int64)
    for i in [2:1:ncol(nd);]
        nd[[1:1:yr+1;],i] = nd[[1:1:yr+1;],i] .* (1.0-μ)
    end
    nd
end;

## aggregate biomass by species
function ΣB(bd::DataFrame, nd::DataFrame, v::Vector{Float64})
    for i in [2:1:ncol(bd);]
        v[i-1] = sum(bd[:,i] .* nd[:,i])
    end
    v
end;

## aggregate leaf area by species
function ΣL(bd::DataFrame, nd::DataFrame, v::Vector{Float64})
    for i in [2:1:ncol(bd);]
        v[i-1] = sum(bd[:,i] .^ (l/b) .* nd[:,i])
    end
    v
end;

## aggregate population by species and add to dynamics data
function Σn!(nd::DataFrame, ndy::DataFrame, yr::Int64)
    for i in [1:1:ncol(nd)-1;]
            ndy[yr,i+1] = sum(nd[:,i+1])
    end
    ndy
end;

## aggregate population by species and add to data
function birth!(nd::DataFrame, v::Vector{Float64}, yr::Int64, F::Float64)
    for i in [2:1:ncol(nd);]
            nd[yr+1, i] = v[i-1] .* F
    end
    nd
end;

function generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.07,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.0005)

    spp_data = DataFrame(C₁ = rand(Uniform(0, C₁max), Nspp),
                         C₂ = rand(Uniform(0, C₂max), Nspp));
    spp_data.τ = broadcast(calc_τ, spp_data.C₁, spp_data.C₂, F, μ, T, b);
    spp_data.Wᵢ = Wmax .* exp.(-tradeoff_exp * spp_data.τ) .+
        rand(Normal(0, tradeoff_sd), Nspp)
    spp_data = sort(spp_data, :Wᵢ, rev = true)
    replace(x -> isless(0, x) ? 0 : x, spp_data.Wᵢ)
    replace(x -> isless(Wmax, x) ? Wmax : x, spp_data.Wᵢ)
    spp_data.spp = Int[1:1:Nspp;];
    return spp_data
end;


## perform iterations and return output
function iterate_water_only_sim(Nyr::Int64, spp_data::DataFrame,
                                biomass_data::DataFrame, biomass_dynamics::DataFrame,
                                n_data::DataFrame, n_dynamics::DataFrame,
                                t::Vector{Float64},
                                v::Vector{Float64}, Nspp::Int64,
                                μ::Float64 = 0.1,
                                W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                                Tvec::Vector{Float64} = repeat([40.0], inner = Nyr))

    ## begin iterating
    for yr in [1:1:Nyr;]
        ## calculate season length for each species
        calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ, W₀vec[yr], Tvec[yr])
        println(t)

        ## calculate average growth rate for each species
        g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

        for s in spp_data.spp
            biomass_data[[1:1:yr;],s+1] = grow(biomass_data[[1:1:yr;],s+1], b, g[s], Tvec[yr])
        end

        ## mortality
        if yr < Nyr
            n_data = die!(n_data, μ, yr)
        else
            n_data = die!(n_data, μ, yr-1)
        end

        ## record biomass in dynamic data
        v = ΣB(biomass_data, n_data, v)
        biomass_dynamics[yr, [2:1:Nspp+1;]] = v

        ## generate new cohort and record population dynamics
        n_data = birth!(n_data, v, yr, F)
        n_dynamics = Σn!(n_data, n_dynamics, yr)

    end

    return [biomass_dynamics, biomass_data, n_dynamics, n_data]

end;


## simulate water only
function sim_water_only(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                        Ninit::Float64, W₀::Float64 = 0.6, T::Float64 = 40.0,
                        rand_W₀::Bool = false,
                        W₀mean::Float64 = 0.6, W₀sd::Float64 = 0.1,
                        rand_T::Bool = false,
                        Tmean::Float64 = 40.0, Tsd::Float64 = 5.0)

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));

    biomass_data = unstack(biomass_data, :spp, :B);
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]]);
    n_data = unstack(n_data, :spp, :N);
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64, n_data[!,[2:ncol(n_data);]]);

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
    n_dynamics = copy(n_data[[1:1:Nyr;],:]);

    ## inital population:
    n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

    t = Vector{Float64}(undef, Nspp);
    v = Vector{Float64}(undef, Nspp);

    if rand_W₀
        W₀vec = rand(Normal(W₀mean, W₀sd), Nyr)
    else
        W₀vec = repeat([W₀], inner = Nyr)
    end

    if rand_T
        Tvec = rand(Normal(Tmean, Tsd), Nyr)
    else
        Tvec = repeat([T], inner = Nyr)
    end

    iterate_water_only_sim(Nyr, spp_data,
                           biomass_data, biomass_dynamics,
                           n_data, n_dynamics,
                           t, v, Nspp, μ, Tvec, W₀vec)

end;

function sigma(μ, l, b)
    x = [0:1:100000;]
    y1 = sum((1-μ).^x .* x.^l)
    y2 = sum((1-μ).^x .* x.^b)
    y2^(l/b) / y1
end;
