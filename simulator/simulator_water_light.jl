
using Plots, QuadGK, DataFrames, Distributions, SpecialFunctions, LazySets, NLsolve
include("utility_functions.jl")

## start with 2 species
Nspp = 1
Nyr = 1000
Ninit = 1
W0 = 0.6
F = 500
mu = 0.03
mu_d = 0.036

spp_data = DataFrame(C1 = rand(Uniform(0, 0.001), Nspp),
                     C2 = rand(Uniform(0, 0.0005), Nspp),
                     h = rand(Uniform(0.4, 0.6), Nspp));
spp_data.tau = broadcast(calc_tau, spp_data.C1, spp_data.C2, F, mu, 40, b)

spp_data.Wi = 0.6 .* exp.(-0.2 * spp_data.tau)
##.+ rand(Normal(0, 0.015), Nspp) ## add random noise to tradeoff
spp_data = sort(spp_data, :Wi, rev = true)
spp_data = spp_data[spp_data.Wi .< W0,:]
Nspp = nrow(spp_data)
spp_data.spp = Int[1:1:Nspp;];
plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

## generate biomass and population data frames
biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                         spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                         B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                         spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                         N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
c_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                   spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                   N = repeat(repeat(Float64[1], inner = Nspp), inner = Nyr+1));
biomass_data = unstack(biomass_data, :spp, :B);
n_data = unstack(n_data, :spp, :N);
c_data = unstack(c_data, :spp, :N);

## dynamics data frames
biomass_dynamics = copy(biomass_data[[1:1:Nyr;],:]);
n_dynamics = copy(n_data[[1:1:Nyr;],:]);

## inital population:
n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

z_data = repeat([0.0], inner = Nyr)

## loop through years
for yr in [1:1:Nyr;]

    live_cohorts = n_data[any.(>(0), eachrow(n_data[:,[2:1:Nspp+1;]])),:rowkey]

    cm = canopy_mem(biomass_data, c_data, n_data) ## determine canopy membership
    c_data = cm[1]
    z_data[yr] = cm[2]

    ## calculate season length for each species
    t = calc_t(sum_L(biomass_data, n_data)[:,2], spp_data.Wi)
    ## calculate average growth rate for each species
    g = Vector{Float64}(undef, Nspp);
    g_d = Vector{Float64}(undef, Nspp);
    for s in spp_data.spp
        g[s] = calc_g(t[s], spp_data[s, :C1], spp_data[s, :C2])
    end
    for s in spp_data.spp
        g_d[s] = calc_g(t[s], spp_data[s, :C1]*0.6, spp_data[s, :C2])
    end

    ## grow species
    for c in live_cohorts
        for s in spp_data.spp
            if c_data[c,s+1] == 1
                new_b = biomass_data[c,s+1]^(1/b) + g[s]*T
            else
                new_b = biomass_data[c,s+1]^(1/b) + g_d[s]*T
            end
            if new_b < 0
                biomass_data[c,s+1] = 0
            else
                biomass_data[c,s+1] = new_b^b
            end
        end
    end

    ## mortality
    if yr < Nyr
        for s in spp_data.spp
            n_data[n_data.rowkey .<= Nyr+1 .&& c_data[:,s+1] .== 1, s+1] = n_data[n_data.rowkey .<= Nyr+1 .&& c_data[:,s+1] .== 1, s+1] .* (1-mu)
            n_data[n_data.rowkey .<= Nyr+1 .&& c_data[:,s+1] .== 0, s+1] = n_data[n_data.rowkey .<= Nyr+1 .&& c_data[:,s+1] .== 0, s+1] .* (1-mu_d)
        end
    else
         for s in spp_data.spp
             n_data[n_data.rowkey .<= Nyr .&& c_data[:,s+1] .== 1, s+1]  = n_data[n_data.rowkey .<= Nyr .&& c_data[:,s+1] .== 1, s+1] .* (1-mu)
             n_data[n_data.rowkey .<= Nyr .&& c_data[:,s+1] .== 0, s+1]  = n_data[n_data.rowkey .<= Nyr .&& c_data[:,s+1] .== 0, s+1] .* (1-mu_d)
         end
    end

    ## reproduction
    total_biomass = sum_B(biomass_data, n_data)
    fecund = total_biomass.sum_B .* F
    n_data[yr+1,[2:1:Nspp+1;]] = fecund

    ## record population and biomass in dynamic data
    biomass_dynamics[yr, [2:1:Nspp+1;]] = total_biomass.sum_B
    nstack = stack(n_data, [2:1:Nspp+1;])
    n_dynamics[yr, [2:1:Nspp+1;]] = combine(groupby(nstack, :variable), :value => sum => :sum_n).sum_n

end

calc_eqN(spp_data)
n_data[nrow(n_data),:]

plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

plot(stack(biomass_dynamics).rowkey, stack(biomass_dynamics).value.+1, group = stack(biomass_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, xlims = (0,Nyr), yscale = :log10)

plot(stack(n_dynamics).rowkey, stack(n_dynamics).value, group = stack(n_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, yscale = :log10, xlims = (0, Nyr))

plot([1:1:Nyr;], z_data, seriestype = :line, legend = :bottomright)

check_eq(biomass_dynamics, 1e-12)
