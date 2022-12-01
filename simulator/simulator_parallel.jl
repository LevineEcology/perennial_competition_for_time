##
using Distributed
ncores = length(Sys.cpu_info())
addprocs(ncores)
@everywhere using Plots, QuadGK, DataFrames, Distributions, Distributed

E = 0.001
G = 0.03
l = 1.5
@everywhere b = 2.5
mu = 0.02
F = 1
T = 40

function calc_g(t, C1, C2, T = 40)
    (t * (C1 + C2) - T * C2)/T
end

##calculate growth time within season
function calc_t(L, ## total leaf area of each spp
                Wi, W0 = 0.8, T = 40)
    df = DataFrame(L = L, Wi = Wi, t = repeat(Float64[T], inner = length(L)), rowkey = [1:1:length(L);])
    df = sort(df, :Wi, rev = true) ## sort in increasing order of drought tolerance
    for s in [1:1:nrow(df);]
        if s == 1
            ti = (W0 - df[s, :Wi]) / (E * sum(df.L))
        else
            ti = (df[s-1, :Wi] - df[s, :Wi]) / (E * sum(df.L[[s:1:nrow(df);]])) + sum(df.t[[1:1:s-1;]])
        end
        if ti < T
            df[s,:t] = ti
        else
            break ## if ti > T then all subsequent times will be as well and output can be returned
        end
    end

    ## return dataframe to original order so output matches input
    df = sort(df, :rowkey)
    return(df.t) ## return times
end

function calc_g_from_L(L, T = 40)
    calc_g(calc_t(L))
end

function calc_g_from_tau(tau)
    calc_g(calc_t(calc_L(tau)))
end

function calc_L(tau)
    L = 1 - tau^3 * G^3
    if L <= 0
        return 0
    else return L
    end
end

function sum_B(biomass_data, n_data)
    bstack = stack(biomass_data, [2:1:ncol(biomass_data);])
    bstack[:, :total_B] = bstack[:,3] .* stack(n_data, [2:1:ncol(biomass_data);])[:,3]
    combine(groupby(bstack, :variable), :total_B => sum => :sum_B)
end

function B2L(B)
    B^(l/b)
end

function calc_tau(C1, C2, F, mu, T, b)
    (1/(C1 + C2)) * (((1/F)^(1/b)) * (((1-exp(-mu*T/b))^2/(exp(-mu*T/b)))) + T*C2)
end

@everywhere function grow(df, g, cohort = 1)
    for s in [1:1:ncol(df)-1;]
        new_b = df[1,s+1]^(1/b) + (g[s])
        if new_b < 0
            df[1,s+1] = 0
        else
            df[1,s+1] = new_b^b
        end
    end
    df[1,:rowkey] = cohort
    df
end

function check_eq(bd, tol = 1e-15)
    check = (Array(bd[nrow(bd),[2:1:ncol(bd);]]) .-
        Array(bd[nrow(bd),[2:1:ncol(bd);]])) ./
        Array(bd[nrow(bd),[2:1:ncol(bd);]])
    if all(check[.!isnan.(check)] .< tol)
        return true
    else return false
    end
end

## start with 2 species
@everywhere Nspp = 7
Nyr = 1000
Ninit = 10
W0 = 0.6
spp = Int[1:1:Nspp;];
spp_data = DataFrame(spp = spp,
                     C1 = rand(Uniform(0, 0.08), Nspp),
                     C2 = rand(Uniform(0, 0.005), Nspp));
spp_data = DataFrame(spp = spp,
                     C1 = [0.02:0.01:0.08;],
                     C2 = repeat([0.005], length([0.02:0.01:0.08;])));
spp_data.tau = broadcast(calc_tau, spp_data.C1, spp_data.C2, 1, 0.05, 40, 3);
spp_data.Wi = 0.6 .* exp.(-0.2*spp_data.tau) .+ rand(Normal(0,0.02), Nspp) .+ 0.2;
plot(spp_data.tau, spp_data.Wi, seriestypes = :scatter, group = spp_data.spp)

## generate biomass and population data frames
biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                         spp = repeat(string.(spp), outer = Nyr+1),
                         B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                         spp = repeat(string.(spp), outer = Nyr+1),
                         N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
biomass_data = unstack(biomass_data, :spp, :B);
n_data = unstack(n_data, :spp, :N);

## dynamics data frames
biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
n_dynamics = copy(n_data[[1:1:Nyr;],:]);

## inital population:
n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);


## loop through years
for yr in [1:1:Nyr;]

    live_cohorts = n_data[any.(>(0), eachrow(n_data[:,[2:1:Nspp+1;]])),:rowkey]
    total_biomass = sum_B(biomass_data, n_data)

    ## calculate season length for each species
    t = calc_t(B2L.(total_biomass.sum_B), spp_data.Wi)

    ## calculate average growth rate for each species
    g = Vector{Float64}(undef, Nspp);
    for spp in spp_data.spp
        g[spp] = calc_g(t[spp], spp_data[spp, :C1], spp_data[spp, :C2])
    end

    ## grow species
    for cohort in live_cohorts
        biomass_data[cohort,:] = grow(DataFrame(biomass_data[cohort,:]), g, cohort)
    end

    ## mortality
    n_data[[1:1:yr;], [2:1:Nspp+1;]] = n_data[[1:1:yr;], [2:1:Nspp+1;]] .* exp(-mu)

    ## reproduction
    total_biomass = sum_B(biomass_data, n_data)
    f = total_biomass.sum_B .* F
    n_data[yr+1,[2:1:Nspp+1;]] = f

    ## record population and biomass in dynamic data
    biomass_dynamics[yr, [2:1:Nspp+1;]] = total_biomass.sum_B
    nstack = stack(n_data, [2:1:Nspp+1;])
    n_dynamics[yr, [2:1:Nspp+1;]] = combine(groupby(nstack, :variable), :value => sum => :sum_n).sum_n

end

plot(spp_data.tau, spp_data.



## loop through years
for yr in [1:1:Nyr;]

    live_cohorts = n_data[any.(>(0), eachrow(n_data[:,[2:1:Nspp+1;]])),:rowkey]
    @everywhere live_cohorts = local(live_cohorts)
    total_biomass = sum_B(biomass_data, n_data)

    ## calculate season length for each species
    t = calc_t(B2L.(total_biomass.sum_B), spp_data.Wi)

    ## calculate average growth rate for each species
    g = Vector{Float64}(undef, Nspp);
    for spp in spp_data.spp
        g[spp] = calc_g(t[spp], spp_data[spp, :C1], spp_data[spp, :C2])
    end
    @everywhere g = local(g)

    ## grow species
    biomass_data[live_cohorts,:] = @distributed (vcat) for cohort in live_cohorts
        grow(DataFrame(biomass_data[cohort,:]), g, cohort)
    end

    @everywhere biomass_data = local(biomass_data)

    ## mortality
    n_data[[1:1:yr;], [2:1:Nspp+1;]] = n_data[[1:1:yr;], [2:1:Nspp+1;]] .* exp(-mu)

    ## reproduction
    total_biomass = sum_B(biomass_data, n_data)
    f = total_biomass.sum_B .* F
    n_data[yr+1,[2:1:Nspp+1;]] = f

    ## record population and biomass in dynamic data
    biomass_dynamics[yr, [2:1:Nspp+1;]] = total_biomass.sum_B
    nstack = stack(n_data, [2:1:Nspp+1;])
    n_dynamics[yr, [2:1:Nspp+1;]] = combine(groupby(nstack, :variable), :value => sum => :sum_n).sum_n

end

plot(spp_data.tau, spp_data.Wi, seriestypes = :scatter, group = spp_data.spp)
plot(stack(biomass_dynamics).rowkey, stack(biomass_dynamics).value, group = stack(biomass_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, xlims = (0,Nyr))

plot(stack(n_dynamics).rowkey, stack(n_dynamics).value, group = stack(n_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, yscale = :log10, xlims = (0, Nyr))

check_eq(biomass_dynamics)

function calc_eqN(sd, F = 1, E = 0.001, W0 = 0.6)
    eqN = Vector{Float64}(undef, nrow(sd))
    sd = sort(sd, :Wi, rev = true)
    for s in sd.spp
        if s == 1
            eqN[s] = (F^((b-1)/b) / E) * ((W0 - sd[s,:Wi])/(sd[s,:tau]) - (sd[s,:Wi] - sd[s+1,:Wi])/(sd[s+1,:tau]-sd[s,:tau]))
        elseif s == nrow(sd)
            eqN[s] = (F^((b-1)/b) / E) * ((W0 - sd[s,:Wi])/(sd[s,:tau]))
        else
            eqN[s] = (F^((b-1)/b) / E) * ((sd[s-1,:Wi] - sd[s,:Wi])/(sd[s,:tau]-sd[s-1,:tau]) - (sd[s,:Wi] - sd[s+1,:Wi])/(sd[s+1,:tau]-sd[s,:tau]))
        end
    endad
    sd.eqN = aeqN
    sort(sd, :spp)
end

calc_eqN(spp_data)
sd = copy(spp_data)
    ((sd[s-1,:Wi] - sd[s,:Wi])/(sd[s,:tau]-sd[s-1,:tau]) )
    (sd[s,:Wi]  - sd[s+1,:Wi])/(sd[s+1,:tau]-sd[s,:tau])

    sd[s, :Wi]
    (sd[s,:Wi] - sd[s+1,:Wi])
    sd[5,:]
    s = 2
    sd[s,:]
