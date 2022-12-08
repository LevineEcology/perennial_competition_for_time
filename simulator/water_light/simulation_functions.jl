## define parameters
##calculate growth time within season
function calc_t!(t::Vector{Float64}, L::DataFrame,
                 Wᵢ::Vector{Float64}, W₀::Float64 = 0.6, T::Float64 = 40.0)

    ## only perform calculations for spp with Wᵢ < W₀
    g0 = findall(W₀ .> Wᵢ)
    t[1:g0[1]-1] .= 0.0
    t[g0] .= T
    if length(g0) != 0
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
    for i in [2:1:ncol(biomass_data);]
        v[i-1] = sum(biomass_data[:,i] .^ (l/b) .* n_data[:,i])
    end
    v
end;

## aggregate population by species and add to dynamics data
function Σn!(n_data::DataFrame, n_dynamics::DataFrame, yr::Int64)
    for i in [1:1:ncol(n_data)-1;]
            n_dynamics[yr,i+1] = sum(n_data[:,i+1])
    end
    nothing
end;

## aggregate population by species and add to dynamics data
function birth!(n_data::DataFrame, v::Vector{Float64}, yr::Int64, F::Float64)
    for i in [2:1:ncol(n_data);]
            n_data[yr+1, i] = v[i-1] .* F
    end
    nothing
end;


## number of species must be even
function generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.07,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.0005)

    spp_data = DataFrame(C₁ = rand(Uniform(0, C₁max), Int64(Nspp/2)),
                         C₂ = rand(Uniform(0, C₂max), Int64(Nspp/2)),
                         h = rand(Uniform(0.4, 0.6), Int64(Nspp/2)));
    spp_data.τ = broadcast(calc_τ, spp_data.C₁, spp_data.C₂, F, μ, T, b);
    spp_data.Wᵢ = Wmax .* exp.(-tradeoff_exp * spp_data.τ) .+
        rand(Normal(0, tradeoff_sd), Int64(Nspp/2))
    spp_data = vcat(spp_data, spp_data)
    spp_data.k = repeat([1, 0.85], inner = Int64(Nspp/2))
    spp_data = sort(spp_data, :Wᵢ, rev = true)
    spp_data.Wᵢ = replace(x -> isless(x, 0) ? 0 : x, spp_data.Wᵢ)
    spp_data.Wᵢ = replace(x -> isless(Wmax, x) ? Wmax : x, spp_data.Wᵢ)
    spp_data.spp = Int[1:1:Nspp;];
    return spp_data
end;

spp_data = generate_spp_data(4)

Nspp = 4
Nyr = 10
Ninit = 1.0

#########
#########
#########
#########
#########

## outline
## 1. loop through disturbance intervals
## 2. at the beginning of each disturbance interval, use accumulated fecundity from previous interval
##    to determine initial population values
## 3. loop through interrain periods.
##    a. calculate total leaf area of each size class
##    b. determine which size classes remain in canopy, if size class partially in canopy, modify leaf area accordingly
##    c. calculate the growing time for each water strategy. assign t = 0, g = 0
##    d. grow species according to growing time
##    e. calculate reproduction

function ΣL(biomass_data::Vector{Float64}, n_data::Vector{Float64}, L::DataFrame)
    for i in 1:length(biomass_data)
        L[i,:L] = sum(biomass_data[i] .^ (l/b) .* n_data[i])
    end
    L
end;

function ΣLhc(L::DataFrame, hcL::Vector{Float64})
    for i in 1:length(unique(L.k))
        hcL[i] = sum(L[L.k .== unique(L.k)[i], :L])
    end
    return hcL
end

function ΣLwc(L::DataFrame, wcL::Vector{Float64})
    for i in 1:length(unique(L.Wᵢ))
        wcL[i] = sum(L[L.Wᵢ .== unique(L.Wᵢ)[i], :L])
    end
    return wcL
end


t₀::Int64 = 10 * 40 ## disturbance interval

## loop through years (disturbance intervals)

n_data = repeat([Ninit], Nspp)
biomass_data = repeat([0.0], Nspp)
L = spp_data[:,[:k, :Wᵢ, :spp]]
L.L = repeat([0.0], inner = Nspp)
hcL = Vector{Float64}(undef, length(unique(spp_data.k)))
wcL = Vector{Float64}(undef, length(unique(spp_data.Wᵢ)))
shortest = 1

g = Vector{Float64}(undef, Nspp);
t = Vector{Float64}(undef, Nspp);
v = Vector{Float64}(undef, Nspp);
Tvec = repeat([40.0], inner = Nyr);
W₀vec = repeat([0.6], inner = Nyr);

for yr in 1:Nyr
    yr = 1

    fc = zeros(Nspp)

    for r in 1:Int64(t₀/T)

        ## sort L by water class
        L = ΣL(biomass_data, n_data, L)
        L = sort(L, :spp)
        wcL = ΣLwc(L, wcL) ## sum leaf area by water strategy



        ## calculate season length for each species
        if any(spp_data.Wᵢ .< W₀vec[yr])
            calc_t!(t, L, W₀vec[yr], Tvec[yr])
        else
            t = repeat([0.0], inner = Nspp)
        end

        ## calculate average growth rate for each species
        g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

        for s in spp_data.spp
            biomass_data[s] = grow(biomass_data[s], b, g[s], Tvec[yr])
        end

        ## determine total leaf area for each height class
        L = ΣL(biomass_data, n_data, L)
        L = sort(L,:k)
        hcL = ΣLhc(L, hcL)

        ## check if canopy closed
        moveon == sum(hcL[shortest:Nspp] > 1.0)
        while moveon == false
            if (sum(hcL[shortest+1:Nspp]) > 1.0)
                shortest = shortest + 1
                moveon = sum(hcL[shortest:Nspp] > 1.0)
            else
                ## calculate fecundity for partially overtopped species
                fc[L[L.k == unique(L.k)[shortest], :spp]] .= fc[L[L.k == unique(L.k)[shortest], :spp]] .+
                    (1 - sum(hcL[shortest+1:Nspp])) .* L[L.k == unique(L.k)[shortest], :L] ./ sum(L[L.k == unique(L.k)[shortest], :L]) .^ (b/l) .* F
                moveon = true
            end
        end

        ## calculate fecundity for canopy species
        fc[L[L.k ∉ unique(L.k)[1:shortest], :spp]] .= fc[L[L.k ∉ unique(L.k)[1:shortest], :spp]] .+
            L[L.k ∉ unique(L.k)[shortest], :L] .^ (b/l) .* F






        sub = L[L.k .== unique(L.k)[2],:]
        Lrat = sub[1,:sum_L] / sum(sub.sum_L)
        out = false
        if sum(L[[2,4], 2]) > 1
            out = true
        elseif sum(L[:,2]) > 1
            L[[1,3],:sum_L] .= [Lrat, 1-Lrat] .* (1-sum(L[[1,2],:sum_L]))
        end
        L[L.sum_L .<= 0, :sum_L] .= 0
        Lc = combine(groupby(L, :Wi), :sum_L => sum => :sum2)
        t = repeat(calc_t_hr(Lc[:,:sum2], Lc.Wi), inner = 2)

        ## calculate average growth rate for each species
        g = Vector{Float64}(undef, Nspp);
        for s in spp_data.spp
            g[s] = calc_g(t[s], spp_data[s, :C1], spp_data[s, :C2]) * spp_data[s,:k]
        end

        if out
            g[[1,3]] .= 0
        end

        ## grow species
        for s in spp_data.spp
            new_b = biomass_data[1,s+1]^(1/b) + g[s]*T
            if new_b < 0
                biomass_data[1,s+1] = 0
            else
                biomass_data[1,s+1] = new_b^b
            end
        end

    end

    ## reproduction
    tL = sum_L_hr(biomass_data, n_data, yr)
    fecund = tL.sum_L .* F
    n_data[yr+1,[2:1:Nspp+1;]] = fecund

end




## perform iterations and return output
function iterate_water_only_sim(Nyr::Int64, spp_data::DataFrame,
                                biomass_data::DataFrame, biomass_dynamics::DataFrame,
                                n_data::DataFrame, n_dynamics::DataFrame,
                                g::Vector{Float64}, t::Vector{Float64},
                                v::Vector{Float64}, Nspp::Int64,
                                μ::Float64 = 0.1,
                                W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                                Tvec::Vector{Float64} = repeat([40.0], inner = Nyr))

    ## begin iterating
    for yr in [1:1:Nyr;]


        ## calculate season length for each species
        if any(spp_data.Wᵢ .< W₀vec[yr])
            calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ, W₀vec[yr], Tvec[yr])
        else
            t = repeat([0.0], inner = Nspp)
        end

        ## calculate average growth rate for each species
        g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

        for s in spp_data.spp
            biomass_data[[1:1:yr;],s+1] = grow(biomass_data[[1:1:yr;],s+1], b, g[s], Tvec[yr])
        end

        ## mortality
        if yr < Nyr
            die!(n_data, μ, yr)
        else
            die!(n_data, μ, yr-1)
        end

        ## record biomass in dynamic data
        v = ΣB(biomass_data, n_data, v)
        biomass_dynamics[yr, [2:1:Nspp+1;]] = v

        ## generate new cohort and record population dynamics
        birth!(n_data, v, yr, F)
        Σn!(n_data, n_dynamics, yr)

    end

    return [biomass_dynamics, n_dynamics, n_data]

end;

## simulate water only
function sim_water_only(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                        Ninit::Float64, μ::Float64 = 0.1, W₀::Float64 = 0.6,
                        T::Float64 = 40.0)

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
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]]);

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
    n_dynamics = copy(n_data[[1:1:Nyr;],:]);

    ## inital population:
    n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

    g = Vector{Float64}(undef, Nspp);
    t = Vector{Float64}(undef, Nspp);
    v = Vector{Float64}(undef, Nspp);

    iterate_water_only_sim(Nyr, spp_data,
                           biomass_data, biomass_dynamics,
                           n_data, n_dynamics,
                           g, t, v, Nspp, μ,
                           repeat([W₀], inner = Nyr),
                           repeat([T], inner = Nyr))

end;

function sigma(μ, l, b)
    x = [0:1:100000;]
    y1 = sum((1-μ).^x .* x.^l)
    y2 = sum((1-μ).^x .* x.^b)
    y2^(l/b) / y1
end;

##---------------------------------------------------------------
## STOCHASTIC SIMULATIONS
##---------------------------------------------------------------



## simulate water only
function sim_water_only_stochastic(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                                   Ninit::Float64, rand_W₀::Bool = true,
                                   rand_T::Bool = false,
                                   W₀mean::Float64 = 0.6, W₀sd::Float64 = 0.05,
                                   Tmean::Float64 = 40.0, Tsd::Float64 = 7.0)

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
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]]);

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
    n_dynamics = copy(n_data[[1:1:Nyr;],:]);

    ## inital population:
    n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

    g = Vector{Float64}(undef, Nspp);
    t = Vector{Float64}(undef, Nspp);
    v = Vector{Float64}(undef, Nspp);

    if rand_W₀
        W₀vec = rand(Normal(W₀mean, W₀sd), Nyr)
        W₀vec = replace(x -> isless(x, 0) ? 1e-10 : x, W₀vec)
    else
        W₀vec = repeat([W₀mean], inner = Nyr)
    end

    if rand_T
        Tvec = rand(Normal(Tmean, Tsd), Nyr)
        Tvec = replace(x -> isless(x, 0) ? 1e-10 : x, Tvec)
    else
        Tvec = repeat([Tmean], inner = Nyr)
    end

    iterate_water_only_sim(Nyr, spp_data,
                           biomass_data, biomass_dynamics,
                           n_data, n_dynamics,
                           g, t, v, Nspp, μ,
                           W₀vec, Tvec)

end;

function plot_simulation_dynamics(results, save::Bool = false)
    pdata = stack(results[2])
    pdata.variable = parse.(Int64, pdata.variable)

    p = plot(pdata.rowkey, pdata.value, group = pdata.variable,line_z = pdata.variable,
             ylim = [0, round(maximum(pdata.value) + 100)], xlim = [0, maximum(pdata.rowkey)],
             seriescolor = my_cgrad, seriestype = :line,
             legend = :topleft, frame = :box, grid = false, linewidth = 1.5)

    if save
        savefig(p, filename)
    end

    return p
end
