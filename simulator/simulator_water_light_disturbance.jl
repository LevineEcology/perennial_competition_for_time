
using Plots, QuadGK, DataFrames, Distributions, SpecialFunctions, LazySets, NLsolve
include("utility_functions.jl")

## start with 2 species
Nspp = 4
t_0 = 40*35
Nyr = 500
Ninit = 1
W0 = 0.6
F = 10
E = 0.02


spp_data = DataFrame(C1 = rand(Uniform(0, 0.001), Int64(Nspp/2)),
                     C2 = rand(Uniform(0, 0.0005), Int64(Nspp/2)),
                     h = rand(Uniform(0.4, 0.6), Int64(Nspp/2)));
spp_data = vcat(spp_data, spp_data)
spp_data.k = repeat([1,0.85], inner = Int64(Nspp/2))

spp_data.tau = broadcast(calc_tau, spp_data.C1, spp_data.C2, F, mu, 40, b)
spp_data.Wi = 0.6 .* exp.(-0.1 * spp_data.tau)


##.+ rand(Normal(0, 0.015), Nspp) ## add random noise to tradeoff
spp_data = sort(spp_data, :Wi, rev = true)
spp_data = spp_data[spp_data.Wi .< W0,:]
Nspp = nrow(spp_data)
spp_data.spp = Int[1:1:Nspp;];
plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                         spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                         N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1));
n_data = unstack(n_data, :spp, :N);

## inital population:
n_data[1,[2:1:Nspp+1;]] = repeat([Ninit], inner = Nspp);

## loop through years
for yr in [1:1:Nyr;]

    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1], inner = Nspp),
                         spp = repeat(string.(spp_data.spp), outer = 1),
                         B = repeat(repeat(Float64[0], inner = Nspp), inner = 1));
    biomass_data = unstack(biomass_data, :spp, :B);

    for r in [1:1:t_0/40;]

        ## calculate season length for each species
        L = sum_L_hr(biomass_data, n_data, yr)
        L.k = spp_data.k
        L.Wi = spp_data.Wi
        L = sort(L,:k)
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
n_data
p_data = stack(n_data);
plot(p_data.rowkey, p_data.value, group = p_data.variable)


spp_data

calc_eqN(spp_data)
n_data[nrow(n_data),:]

plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

plot(stack(biomass_dynamics).rowkey, stack(biomass_dynamics).value.+1, group = stack(biomass_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, xlims = (0,Nyr), yscale = :log10)

plot(stack(n_dynamics).rowkey, stack(n_dynamics).value, group = stack(n_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, yscale = :log10, xlims = (0, Nyr))

plot([1:1:Nyr;], z_data, seriestype = :line, legend = :bottomright)

check_eq(biomass_dynamics, 1e-12)



x1 = [1:1:50;]
x2 = [51:1:100;]

g1 = 10
g2 = 5

y1 = x1 .* g1
y2 = y1[50] .+ (x2 .- 51) .* g2

plot([x1, x2 .- 1], [y1, y2])
plot([x1, x2 .- 1], [y1 .^ 1.5, y2 .^ 1.5])

y3 = y1[50]^1.5 .+ ((x2 .- 51) .* g2) .^ 1.5


plot!([x1, x2 .- 1], [y1 .^ 1.5, y3])
