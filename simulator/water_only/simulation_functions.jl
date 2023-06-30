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

function calc_w(w_cur, t, T, biomass_data, n_data, v, Wᵢ)

    v = ΣL(biomass_data, n_data, v)

    sb = findall(Wᵢ .< w_cur)
    if isempty(sb)
        return w_cur
    elseif all(t[sb] .< T)
        return Wᵢ[length(t)]
    elseif all(t[sb] .== T)
         return w_cur - E * sum(v) * T
    else
        fs = minimum(findall(t .== T))
        return Wᵢ[fs-1] - E * sum(v[fs:length(t)]) * T
    end
end

function sigma(μ, ex)
    Float64(polylog(-ex, exp(-μ)))
end;

## calculate average growth rates
function calc_g(t::Float64, C₁::Float64, C₂::Float64, T::Float64 = 40.0)
    (t * C₁ - (T-t) * C₂)/T
end;

function calc_eq_leaf_area(eqN::Vector{Float64}, F::Float64, μ::Float64)
    eqN ./ (μ .* F)
end

function calc_τ(C₁::Float64, C₂::Float64, F::Float64, μ::Float64, T::Float64, b::Float64 = 2.5)
    s = sigma(μ*T, b)
    (T/(C₁ + C₂)) * ((1/(F * T^(b+1) * s))^(1/b) + C₂)
end

## grow plant biomass
function grow(biomass::Vector{Float64}, b::Float64, g::Float64, T::Float64)
    replace(x -> isless(x, 0) ? 0 : x, (biomass .^ (1/b) .+ (g .* T))) .^ b
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
function birth!(n_data::DataFrame, v::Vector{Float64},
                yr::Int64, F::Float64, T::Float64)
    for i in 2:ncol(n_data)
            n_data[yr+1, i] = v[i-1] .* F .* T
    end
    nothing
end;

function generate_spp_data(Nspp::Int64, Wmax::Float64 = 0.8, T::Float64 = 40.0,
                           F::Float64 = 100.0, μ::Float64 = 0.1,
                           b::Float64 = 2.5, tradeoff_exp::Float64 = 0.4,
                           tradeoff_sd::Float64 = 0.07,
                           C₁max::Float64 = 0.01, C₂max::Float64 = 0.0005)

    spp_data = DataFrame(C₁ = rand(Uniform(0.0000005, 0.01), Nspp) .* 365.0,
                         C₂ = repeat([0.001], Nspp) .* 365.0);
    spp_data.τ = broadcast(calc_τ, spp_data.C₁, spp_data.C₂, F, μ, T, b);
    spp_data.Wᵢ = Wmax .* exp.(-0.4 .* (5 .- spp_data.C₁)) .+
        rand(Normal(0, tradeoff_sd), Nspp)
    spp_data = sort(spp_data, :Wᵢ, rev = true)
    spp_data.Wᵢ = replace(x -> isless(x, 0) ? 1e-4 : x, spp_data.Wᵢ)
    spp_data.Wᵢ = replace(x -> isless(Wmax, x) ? Wmax : x, spp_data.Wᵢ)
    spp_data.spp = Int[1:1:Nspp;];
    return spp_data
end;

## perform iterations and return output
function iterate_water_only_sim(Nyr::Int64, spp_data::DataFrame,
                                biomass_data::DataFrame, biomass_dynamics::DataFrame,
                                n_data::DataFrame, n_dynamics::DataFrame,
                                r_data::DataFrame, w_data::Vector{Float64},
                                g::Vector{Float64}, t::Vector{Float64},
                                v::Vector{Float64}, Nspp::Int64, mt::Matrix{Float64},
                                μ::Float64 = 0.1, F::Float64 = 10.0,
                                W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                                Tvec::Vector{Float64} = repeat([40.0], inner = Nyr),
                                θ_fc::Float64 = 0.4, pb::Bool = true)

    ## initialize soil water content at field capacity
    w = θ_fc
    w_in_data = copy(w_data)

    ## begin iterating
    if pb
        prog = ProgressBar(total = Nyr)
    end
    for yr in 1:Nyr

        ## calculate initial water content
        w = minimum([w + W₀vec[yr], θ_fc])
        w_in_data[yr] = w

        ## calculate growth time for each species
        if any(spp_data.Wᵢ .< w)
            calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ,
                    w, Tvec[yr])
        else
            t = repeat([0.0], inner = Nspp)
        end

        ## calculate average growth rate for each species
        g = calc_g.(t, spp_data[:, :C₁], spp_data[:, :C₂], Tvec[yr])

        for s in spp_data.spp
            biomass_data[[1:1:yr;], s+1] =
                grow(biomass_data[[1:1:yr;], s+1], b, g[s], Tvec[yr])
        end

        die!(n_data, r_data, yr, mt)

        ## recalculate soil water content
        w = maximum([calc_w(w, t, Tvec[yr], biomass_data, n_data, v, spp_data.Wᵢ), 0.0])
        w_data[yr] = w

        v = ΣB(biomass_data, n_data, v)
        biomass_dynamics[yr, [2:1:Nspp+1;]] = v

        ## generate new cohort and record population dynamics
        birth!(n_data, v, yr, F, Tvec[yr])
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

    return [biomass_dynamics, n_dynamics, n_data, r_data, w_data, w_in_data]

end;


## simulate water only
function sim_water_only(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                        Ninit::Float64, μ::Float64 = 0.15, F::Float64 = 10.0,
                        P::Int64 = 40, mean_p::Float64 = 16.0, θ_fc = 0.4,
                        mt::Matrix{Float64} = zeros(1,1), pb::Bool = true)

    Nyr = Nyr * P

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

    biomass_data = unstack(biomass_data, :spp, :B)
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]])
    n_data = unstack(n_data, :spp, :N)
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]])
    r_data = unstack(r_data, :spp, :N)
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]])

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[1:Nyr,:])
    n_dynamics = copy(n_data[1:Nyr,:])

    ## inital population:
    n_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
    r_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)

    g = Vector{Float64}(undef, Nspp)
    t = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    w_data = Vector{Float64}(undef, Nyr)

    T = 1.0 / P
    W₀ = mean_p / P

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, repeat([T], inner = Nyr))
    end

    iterate_water_only_sim(Nyr, spp_data,
                           biomass_data, biomass_dynamics,
                           n_data, n_dynamics,
                           r_data, w_data,
                           g, t, v, Nspp, mt, μ, F,
                           repeat([minimum([θ_fc, W₀])], inner = Nyr),
                           repeat([T], inner = Nyr), θ_fc,
                           pb)

end;

##---------------------------------------------------------------
## STOCHASTIC SIMULATIONS
##---------------------------------------------------------------
##

## generate a year's worth of interval lengths given P
function generate_intervals(Pmean, cluster = false)
    inter = rand(Exponential(1/Pmean), Pmean)
    sum(inter)
    if sum(inter) > 1.0
        rs = 0.0
        for i in 1:length(inter)
            rs = rs + inter[i]
            if rs > 1.0
                inter = inter[1:i]
                inter[i] = 1.0 - sum(inter[1:i-1])
            end
            break
        end
    end
    if cluster
        sort!(vec(inter))
    end
    return inter
end

function generate_rainfall_regime(Nyr::Int64, Pmean::Int64, Pdisp::Float64, map_mean::Float64, map_var::Float64, cluster::Bool = false)

    var = Pmean + 1 / Pdisp * Pmean ^ 2
    p = (var - Pmean) / var
    Plist = rand(NegativeBinomial(Pdisp, 1-p), Nyr)
    Plist[Plist .== 0] .= 1

    precip_list = rand(Normal(map_mean, map_var), Nyr)

    Preal = Float64[]
    Tlist = Int64[]
    for yr in 1:Nyr
        new_int = generate_intervals(Plist[yr])
        Tlist = vcat(Tlist, new_int)
        Preal = vcat(Preal, length(new_int))
    end
    Tlist
    Preal = Int.(Preal)

    W₀list = Vector{Float64}(undef, length(Tlist))
    i = 1
    for yr in 1:Nyr
        W₀list[i:i+Preal[yr]-1] .= Tlist[i:i+Preal[yr]-1] .* precip_list[yr]
        i = i+Preal[yr]
    end

    return [W₀list, Tlist, Preal]

end

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


function plot_rainfall_regime(rainfall_regime)

    t = Vector{Float64}(undef, length(rainfall_regime[1]))
    t[1] = rainfall_regime[2][1]
    for i in 2:length(t)
        t[i] = t[i-1] + rainfall_regime[2][i]
    end

    w = copy(rainfall_regime[1])

    p = plot(t, w, seriestype = :scatter, color = :black, markersize = 0,
             frame = :box, legend = :none)

    for i in 1:length(t)
       p = plot!(vcat(t[i], t[i]), vcat(0.0, w[i]), color = :black, linewidth = 2)
    end

    return p

end


##---------------------------------------------------------------
## PLOTTING UTILITY
##---------------------------------------------------------------

function plot_simulation_dynamics(results, save::Bool = false, filename = "")
    pdata = stack(results[2])
    pdata.variable = parse.(Int64, pdata.variable)

    p = plot(pdata.rowkey, pdata.value, group = pdata.variable,line_z = pdata.variable,
             ylim = [0, round(maximum(pdata.value))+1.0], xlim = [0, maximum(pdata.rowkey)],
             seriescolor = my_cgrad, seriestype = :line,
             legend = :none, frame = :box, grid = false, linewidth = 2.5)

    pdata.value

    if save
        savefig(p, filename)
    end

    return p
end

function plot_simulation_dynamics_stochastic(results, rainfall_regime, save::Bool = false, filename = "")

    t = Vector{Float64}(undef, length(rainfall_regime[1]))
    t[1] = rainfall_regime[2][1]
    for i in 2:length(t)
        t[i] = t[i-1] + rainfall_regime[2][i]
    end

    pdata = stack(results[2])
    pdata.variable = parse.(Int64, pdata.variable)

    p = plot(repeat(t, outer = length(unique(pdata.variable))), pdata.value, group = pdata.variable,line_z = pdata.variable,
             ylim = [0, round(maximum(pdata.value))+1.0], xlim = [0, maximum(t)],
             seriescolor = my_cgrad, seriestype = :line,
             legend = :none, frame = :box, grid = false, linewidth = 2.5)
    p2 = plot_rainfall_regime(rainfall_regime)
    p2 = plot!(t, out[5], color = :blue, ylim = [0, 0.4], xlim = [0, maximum(t)],)

    p3 = plot(p, p2, layout = (2,1))

    if save
        savefig(p, filename)
    end

    return p3
end
