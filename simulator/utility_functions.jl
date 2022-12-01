using Base: _grow_filter!
E::Float64 = 0.02
G::Float64 = 0.03
l::Float64 = 1.5
b::Float64 = 2.5
mu::Float64 = 0.1
F::Float64 = 100
T::Float64 = 40.0

function calc_t_hr(L, ## total leaf area of each spp
                Wi, W0 = 0.6, T = 40)
    df = DataFrame(L = L, Wi = Wi, t = repeat(Float64[T], inner = length(L)), rowkey = [1:1:length(L);])
    df = sort(df, :Wi, rev = true) ## sort in increasing order of drought tolerance
    for s in [1,2]
        if s == 1
            ti = (W0 - df[s, :Wi]) / (E * sum(df.L))
        else
            ti = ((df[s-1, :Wi] - df[s, :Wi]) / (E * sum(df.L[[s:1:nrow(df);]]))) + df.t[s-1]
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

function sum_L_hr(biomass_data, n_data, yr)
    bstack = stack(biomass_data, [2:1:ncol(biomass_data);])
    bstack.L = bstack.value.^(l/b)
    bstack[:, :total_L] = bstack[:,:L] .* Array(n_data[yr,[2:1:ncol(n_data);]])
    combine(groupby(bstack, :variable), :total_L => sum => :sum_L)
end

function calc_height(biomass_data)
    height_data = Array{Int64}(biomass_data[:,1])
    for i in [1:1:Nspp;]
        height_data = hcat(height_data, biomass_data[:,i+1] .^ spp_data[i,:h])
    end
    DataFrame(height_data, names(biomass_data))
end

function canopy_mem(biomass_data, c_data, n_data)
    height_data = calc_height(biomass_data)
    h_order = stack(height_data, [2:1:ncol(height_data);])
    h_order = sort(h_order, :value, rev = true)
    h_order[:,:sum_canopy] = repeat([0.0], inner = nrow(h_order))
    n_stack = stack(n_data, [2:1:ncol(biomass_data);])
    for i in [2:1:nrow(h_order);]
        h = spp_data[parse(Int64, h_order[i,2]),:h]
        h_order[i,:sum_canopy] = h_order[i-1,:sum_canopy] +
            (h_order[i,:value]^(l/h) * n_stack[i,3])
        if h_order[i,:sum_canopy] < 1
            c_data[c_data.rowkey .== h_order[i,:rowkey], parse(Int64, h_order[i,2])+1] .= 1
        else
            c_data[c_data.rowkey .== h_order[i,:rowkey], parse(Int64, h_order[i,2])+1] .= 0
        end
    end
    k = findfirst(>=(1), h_order.sum_canopy)
    if isnothing(k)
        zstar = 0
    else
        zstar = h_order[findfirst(>=(1), h_order.sum_canopy), :value]
    end
    [c_data, zstar]
end


function B2L(B)
    B^(l/b)
end

function lrs(biomass_data)
    bstack = stack(biomass_data, [2:1:ncol(biomass_data);])
    x = [(nrow(bstack))-1:-1:0;]
    y = (1-mu).^x
    bstackcor = bstack[bstack.variable .== "1",:value] .* y
    sum(bstackcor)
end

function grow(df, g, cohort = 1)
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
## for CNDD
## solves for equilibrium densities for annual system with CNDD
function f(f, x, Wi = spp_data.Wi, tau = spp_data.tau)
    G = (1 ./ (tau .* F)).^(1/b)
    sum1 = Array{Float64}(undef, length(x)) ## initialize empty array for summed terms
    sum2 = Array{Float64}(undef, length(x))
    for i in [1:1:length(x);]
        sum1[i] = (x[i] / (1 + (a * x[i]))) * G[i]^l
    end
    for j in [1:1:length(x);]
        if j == 1
            sum2[j] = (W0 - Wi[j]) / (E * sum(sum1[[j:1:length(x);]]))
        else
            sum2[j] = (Wi[j-1] - Wi[j]) / (E * sum(sum1[[j:1:length(x);]]))
        end
    end
    for k in [1:1:length(x);]
        f[k] = F * G[k]^b * (sum(sum2[[1:1:k;]])) * (1 / (1 + a * x[k])) - 1
    end
    f
end

function calc_eqN_annuals(sd, F = 1, E = E, W0 = 0.6)
    eqN = Vector{Float64}(undef, nrow(sd))
    sd = sort(sd, :Wi, rev = true)

    fe = feas(sd)
    eqN[(!in).(sd.spp,Ref(fe))] .= 0

    for s in fe
        if s == fe[1]
            eqN[s] = (((l+1)) / (E * sd.G[s]^l)) * ((W0 - sd[s,:Wi])/(sd[s,:tau]) - ((sd[s,:Wi] - sd[s+1,:Wi]) / (sd[s+1,:tau] - sd[s,:tau])))
        elseif s == fe[length(fe)]
            eqN[s] = (((l+1)) / (E * sd.G[s]^l)) * ((sd[s-1,:Wi] - sd[s,:Wi])/(sd[s,:tau] - sd[s-1, :tau]))
        else
            eqN[s] = (((l+1)) / (E * sd.G[s]^l)) * (((sd[s-1,:Wi] - sd[s,:Wi])/(sd[s,:tau]-sd[s-1,:tau])) - ((sd[s,:Wi] - sd[s+1,:Wi])/(sd[s+1,:tau]-sd[s,:tau])))
        end
    end

    sd.eqN = eqN
    sort(sd, :spp)
end
