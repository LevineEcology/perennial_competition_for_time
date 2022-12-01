
using Plots, QuadGK, DataFrames, Distributions, Random,
    SpecialFunctions, LazySets, NLsolve, CSV
include("utility_functions.jl")

mu = 0.1
W0 = 0.6
iter = 200
Nspp = 5

results = DataFrame(prop = Array{Float64}(undef, 3*iter*length([0:.003:0.05;])),
                    a = repeat([0:0.003:0.05;], inner = iter*3),
                    Nspp = repeat([100, 150, 200], outer = iter * length([0:.003:0.05;])))

results = DataFrame(prop = Array{Float64}(undef, iter),
                    a = repeat([0], inner = iter),
                    Nspp = repeat([40], inner = iter))

for it in [1:1:nrow(results);]

    Nspp = results[it, :Nspp]

    ##spp_data = DataFrame(G = rand(Uniform(0.1, 0.5), Nspp));
    ##spp_data.tau .= 1 ./ (F .* spp_data.G.^b)
    spp_data = DataFrame(tau = [1:100/Nspp:100;]);
    spp_data.G .= (1 ./ (F .* spp_data.tau))

    spp_data = sort(spp_data, :tau)
    spp_data.Wi = Array{Float64}(undef, nrow(spp_data))

    for s in [1:1:nrow(spp_data);]
        if s == 1
            spp_data.Wi[s] = 0.6 - rand(Exponential(0.6 / Nspp))
        else
            spp_data.Wi[s] = spp_data.Wi[s-1] - rand(Exponential(0.6 / Nspp))
        end
    end

    spp_data = sort(spp_data, :Wi, rev = true)
    spp_data = spp_data[spp_data.Wi .< W0 .&& spp_data.Wi .> 0,:]
    Nspp = nrow(spp_data)
    spp_data.spp = Int[1:1:Nspp;];
    plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

    function f(f, x, Wi = spp_data.Wi, a = results[it, :a])
        sum1 = Array{BigFloat}(undef, length(x)) ## initialize empty array for summed terms
        sum2 = Array{BigFloat}(undef, length(x))
        #x[x .<0] .= 0
        for i in [1:1:length(x);]
            sum1[i] = x[i] *  spp_data.G[i]^l
        end
        for j in [1:1:length(x);]
            if j == 1
                sum2[j] = (W0 - Wi[j]) / (sum(sum1[[j:1:length(x);]]))
            else
                sum2[j] = (Wi[j-1] - Wi[j]) / (sum(sum1[[j:1:length(x);]]))
            end
        end
        for k in [1:1:length(x);]
            f[k] = (((l+1) * F * (spp_data.G[k]^b) / E) * (sum(sum2[[1:1:k;]])) * (x[k] / (1 + a * x[k])))
        end
        f
    end

    x = fixedpoint(f, repeat([1.0], inner = nrow(spp_data)),
                   iterations = 500000, ftol = 1e-10)

    if converged(x)
        results[it, :prop] = count(x.zero .> 1e-10)
    end
    println("finished iteration: ")
    println(it)

end


plot(results.prop, seriestype = :histogram, bins = 20)








CSV.write("data/cndd_coex_annuals.csv", results)

res = CSV.File("data/cndd_coex_annuals.csv") |> DataFrame

p_data = combine(groupby(res, ["a", "Nspp"]), :prop => mean => :mean)

plot(p_data.a, p_data.mean, seriestype = :scatter, group = p_data.Nspp,
     linewidth = 3, grid = false, framestyle = :box, legend = :bottomright, ylim = [0, 200])









## simulator


## start with 2 species
Nspp = 50
Nyr = 2000
Ninit = 1
W0 = 0.6
spp_data = DataFrame(C1 = rand(Uniform(0, 0.01), Nspp),
                     C2 = rand(Uniform(0, 0.0005), Nspp));
spp_data.tau = broadcast(calc_tau, spp_data.C1, spp_data.C2, 1, mu, 40, b)
mu = 0.1

#spp_data.Wi = (0.03 .* spp_data.tau .- sqrt(0.55)).^2 .+ 0.05
spp_data.Wi = 0.6 .* exp.(-0.1 * spp_data.tau) .+ rand(Normal(0, 0.015), Nspp)
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
    t = calc_t(sum_L(biomass_data, n_data)[:,2], spp_data.Wi)
    ## calculate average growth rate for each species
    g = Vector{Float64}(undef, Nspp);
    for s in spp_data.spp
        g[s] = calc_g(t[s], spp_data[s, :C1], spp_data[s, :C2])
    end

    ## grow species
    for c in live_cohorts
        for s in spp_data.spp
            new_b = biomass_data[c,s+1]^(1/b) + (g[s]*T)
            if new_b < 0
                biomass_data[c,s+1] = 0
            else
                biomass_data[c,s+1] = new_b^b
            end
        end
    end

    ## mortality
    if yr < Nyr
        n_data[[1:1:yr+1;], [2:1:Nspp+1;]] = n_data[[1:1:yr+1;], [2:1:Nspp+1;]] .* (1-mu)
    else
        n_data[[1:1:yr;], [2:1:Nspp+1;]] = n_data[[1:1:yr;], [2:1:Nspp+1;]] .* (1-mu)
    end

    ## reproduction
    total_biomass = sum_B(biomass_data, n_data)
    f = total_biomass.sum_B .* F
    f = f ./ (1 .+ (0.05 .* f))
    n_data[yr+1,[2:1:Nspp+1;]] = f

    ## record population and biomass in dynamic data
    biomass_dynamics[yr, [2:1:Nspp+1;]] = total_biomass.sum_B
    nstack = stack(n_data, [2:1:Nspp+1;])
    n_dynamics[yr, [2:1:Nspp+1;]] = combine(groupby(nstack, :variable), :value => sum => :sum_n).sum_n

end

t
spp_data.tau

calc_eqN(spp_data)
n_data[nrow(n_data),:]

plot(vcat(0,spp_data.tau), vcat(W0, spp_data.Wi), seriestypes = :scatter)

plot(stack(biomass_dynamics).rowkey, stack(biomass_dynamics).value.+1, group = stack(biomass_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, xlims = (0,Nyr), yscale = :log10)

plot(stack(n_dynamics).rowkey, stack(n_dynamics).value, group = stack(n_dynamics).variable,
     linewidth = 3, grid = false, framestyle = :box, yscale = :log10, xlims = (0, Nyr))

check_eq(biomass_dynamics, 1e-12)



ex = function(z, lambda = 5)
    -lambda^2 * z * exp(-lambda*z)
end

x = [0:0.001:20;]
y = (1 ./ (x .+ 1))
plot(x,y, seriestype = :scatter)
ex(0.001)
