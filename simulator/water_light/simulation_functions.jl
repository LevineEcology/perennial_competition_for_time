
E = 0.05
## define parameters
## calculate growth time within season
function calc_t!(t::Vector{Float64}, wcL::Vector{Float64}, Wᵢ::Vector{Float64},
                 W₀::Float64 = 0.6, T::Float64 = 40.0)

    ## only perform calculations for spp with Wᵢ < W₀
    g0 = findall(W₀ .> Wᵢ)
    t[1:g0[1]-1] .= 0.0
    t[g0] .= T
    if length(g0) != 0
        for s in g0[1]:length(Wᵢ)
            if s == 1
                t[s] = (W₀ - Wᵢ[s]) / (E * sum(wcL))
            else
                t[s] = ((Wᵢ[s-1] - Wᵢ[s]) / (E * sum(wcL[[s:1:length(wcL);]]))) + t[s-1]
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
function calc_g(t::Float64, C₁::Float64, C₂::Float64, k::Float64, T::Float64 = 40.0)
    k * (t * (C₁ + C₂) - T * C₂)/T
end;

function calc_τ(C₁::Float64, C₂::Float64, F::Float64, T::Float64, expon::Float64)
  (1.0 / (C₁ + C₂)) * ((1.0 / (F^(1/b) * 30.0)) + T * C₂)
end

## grow plant biomass
function grow(biomass::Float64, b::Float64, g::Float64, T::Float64)
    replace(x -> isless(x, 0.0) ? 0.0 : x, [(biomass ^ (1/b) + (g .* T))])[1] ^ b
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

    spp_data = DataFrame(C₁ = rand(Uniform(0, C₁max), Int64(Nspp/4)),
                         C₂ = rand(Uniform(0, C₂max), Int64(Nspp/4)),
                         h = rand(Uniform(0.4, 0.6), Int64(Nspp/4)));
    spp_data.τ = broadcast(calc_τ, spp_data.C₁, spp_data.C₂, F, T, b);
    spp_data.Wᵢ = Wmax .* exp.(-tradeoff_exp .* spp_data.τ) .+
        rand(Normal(0, tradeoff_sd), Int64(Nspp/4))
    spp_data = vcat(spp_data, spp_data, spp_data, spp_data)
    spp_data.k = repeat([0.025:0.025:0.1;], inner = Int64(Nspp/4))
    spp_data = sort(spp_data, :Wᵢ, rev = true)
    spp_data.Wᵢ = replace(x -> isless(x, 0) ? 0 : x, spp_data.Wᵢ)
    spp_data.Wᵢ = replace(x -> isless(Wmax, x) ? Wmax : x, spp_data.Wᵢ)
    spp_data.spp = Int[1:1:Nspp;];
    return spp_data
end;


F
Nspp = 8
Nyr = 100
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

function ΣLwc(L::DataFrame, wcL::Vector{Float64}, out)
    for i in 1:length(unique(L.Wᵢ))
        wcL[i] = sum(L[L.Wᵢ .== unique(L.Wᵢ)[i] .&& n_data .> 0.0 .&& L.spp .∉ Ref(out), :L])
    end
    return wcL
end

t₀::Int64 = 20 * 40 ## disturbance interval
F = 100.0

theme(:dark)
my_cgrad = cgrad(:acton)

spp_data = generate_spp_data(8, 0.6, 40.0, 100.0, 0.1, 2.5, 0.4, 0.0)
out = sim_water_light(spp_data, 100, nrow(spp_data), 1.0, 0.8, 40.0, 25.0*40.0)
plot_simulation_dynamics(out)
spp_data

## simulate water only
function sim_water_light(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                         Ninit::Float64 = 1.0, W₀ = 0.6,
                         T::Float64 = 40.0, t₀ = 20.0*40.0)

    ## setup for simulations
    n_data = repeat([Ninit], Nspp)
    biomass_data = repeat([0.0], Nspp)
    L = spp_data[:,[:k, :Wᵢ, :spp]]
    L.L = repeat([0.0], inner = Nspp)
    L.canopy_status = repeat(["NA"], inner = Nspp)
    hcL = zeros(length(unique(spp_data.k)))
    wcL = zeros(length(unique(spp_data.Wᵢ)))

    g = Vector{Float64}(undef, Nspp);
    t = Vector{Float64}(undef, length(unique(spp_data.Wᵢ)));
    tfull = Vector{Float64}(undef, Nspp);
    Tvec = repeat([T], inner = Nyr);
    W₀vec = repeat([W₀], inner = Nyr);

    n_dynamics = Array{Float64}(undef, Nyr, Nspp)

    iterate_water_light_sim(Nyr, spp_data, t₀,
                           biomass_data,
                           n_data, n_dynamics,
                           g, t, L, hcL, wcL, tfull, Nspp,
                           repeat([W₀], inner = Nyr),
                           repeat([T], inner = Nyr))

end;

## perform iterations and return output
function iterate_water_light_sim(Nyr::Int64, spp_data::DataFrame, t₀::Float64,
                                 biomass_data::Vector{Float64},
                                 n_data::Vector{Float64}, n_dynamics::Matrix{Float64},
                                 g::Vector{Float64}, t::Vector{Float64}, L::DataFrame,
                                 hcL::Vector{Float64}, wcL::Vector{Float64}, tfull::Vector{Float64},
                                 Nspp::Int64,
                                 W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                                 Tvec::Vector{Float64} = repeat([40.0], inner = Nyr))

    shortest = 1

    for yr in 1:Nyr

        fc = zeros(Nspp)
        biomass_data = zeros(Nspp)
        out = Vector{Int64}(undef,0)
        hcL = Vector{Float64}(undef, length(unique(spp_data[:, :k])))

        for r in 1:Int64(t₀/T)

            #### WATER

            ## sort L by water class
            L = sort(L, :spp)
            L = ΣL(biomass_data, n_data, L)
            wcL = zeros(length(unique(spp_data.Wᵢ)))
            wcL = ΣLwc(L, wcL, out) ## sum leaf area by water strategy

            ## calculate season length for each species
            if any(spp_data.Wᵢ .< W₀vec[yr])
                calc_t!(t, wcL, unique(L.Wᵢ), W₀vec[yr], Tvec[yr])
                for i in 1:length(t)
                    tfull[L.Wᵢ .== unique(L.Wᵢ)[i]] .= t[i]
                end
            else
                tfull = repeat([0.0], inner = Nspp)
            end

            ## calculate average growth rate for each species
            g = calc_g.(tfull, spp_data.C₂, spp_data.C₂, spp_data.k, Tvec[yr])

            for s in [1:1:Nspp;][[1:1:Nspp;] .∉ Ref(out)]
                biomass_data[s] = grow(biomass_data[s], b, g[s], Tvec[yr])
            end

            ## determine total leaf area for each height class
            L = sort(L, :spp)
            L = ΣL(biomass_data, n_data, L)
            L = sort(L, :k, rev = true)
            hcL = ΣLhc(L, hcL)

            rsum = copy(hcL)
            for i in  reverse(2:length(hcL))
                rsum[i-1] = hcL[i-1] + rsum[i]
            end
            rsum

            for i in [length(hcL):-1:1;]
                if rsum[i] < 1.0 && rsum[i] > 0.0
                    L[L.k .== unique(L.k)[i], :canopy_status] .= "c"
                elseif rsum[i] > 1.0 && i == length(hcL)
                    L[L.k .== unique(L.k)[i], :canopy_status] .= "c"
                elseif rsum[i] > 1.0 && all(L[L.k .∈ Ref(unique(L.k)[i+1:length(hcL)]), :canopy_status] .== "u")
                    L[L.k .== unique(L.k)[i], :canopy_status] .= "c"
                elseif rsum[i] > 1.0 && rsum[i+1] < 1.0 && rsum[i+1] > 0.0
                    L[L.k .== unique(L.k)[i], :canopy_status] .= "p"
                else
                    L[L.k .== unique(L.k)[i], :canopy_status] .= "u"
                end
            end

            for i in 1:nrow(L)
                if n_data[L[i, :spp]] > 0.0
                    if L[i, :canopy_status] == "c"
                        fc[L[i, :spp]] = fc[L[i, :spp]] +
                            L[i, :L] / sum(L[L.canopy_status .== "c", :L]) * F
                    elseif L[i, :canopy_status] == "p"
                        fc[L[i, :spp]] = fc[L[i, :spp]] +
                            L[i, :L] / sum(L[L.canopy_status .== "p", :L]) *
                            (1 - sum(L[L.canopy_status .== "c", :L])) .^ (b/l) .* F
                    else
                        biomass_data[L[i, :spp]] = 0.0
                        out = unique(append!(out, L[i, :spp]))
                    end
                end
            end

        end

        n_dynamics[yr,:] = fc
        n_data = fc

    end

    return n_dynamics

end;




function plot_simulation_dynamics(results, save::Bool = false)

    yd = Vector{Float64}(undef, size(results, 1) * size(results,2))
    xd = repeat(collect(1:size(results,1)), outer = size(results,2))
    gd = repeat(collect(1:size(results,2)), inner = size(results,1))

    for i in 1:size(results,2)
        yd[(i-1)*size(results,1)+1:i*size(results,1)] .= results[:,i]
    end

    p = plot(xd, yd, group = gd, line_z = gd,
             ylim = [1, round(maximum(yd) + 100)], xlim = [0, maximum(xd)],
             seriescolor = my_cgrad, seriestype = :line, lw = 2,
             legend = :topleft, frame = :box, grid = false, linewidth = 1.5,
             colorbar = false)

    if save
        savefig(p, filename)
    end

    return p
end
