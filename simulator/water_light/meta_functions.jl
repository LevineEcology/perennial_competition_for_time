##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- meta_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## November 2022
##
## This script defines functions used to perform analyses across environmental
## gradients, typically performing calculations over a range of environmental values
## for many communities.
##---------------------------------------------------------------

"""
    multi_eq(Nspp::Int64 = 10, Niter::Int64 = 10, θ_fc = 0.4,
             mintotal::Float64 = 0.1, maxtotal::Float64 = 3.0,
             lengthtotal::Int64 = 10,
             minP::Int64 = 8, maxP::Int64 = 50, lengthP::Int64 = 10,
             F::Float64 = 10.0, μ::Float64 = 0.03, n_hts::Int64 = 3, uf::Float64 = 0.1,
             constant_ss::Bool = false, ss::Float64 = 0.02)

Calculates equilibrium population density, equilibrium biomass, and total ecosystem evapotranspiration
over a gradient of precipitation regimes supplied by the user. `Nspp` sets the size of the communities, while
`Niter` dictates how many communities will be generated.

The precipitation regimes are defined by two parameters, `total`, which describes the total monthly
precipitation in the volumetric water content equivalent, and `P` whichdescribes the number of storms that occur
per month. Both of these parameters are assumed to be invariable. Therange of the parameters over which
calculations are performed are set by `mintototal`/`maxtotal` and `minP`/`maxP`, while the number of values over
which to perform the calculations is specified by `lengthtotal` and `lengthP`. The `constant_ss` flag lets the
user override the parameterization and dictate a constant storm size, which is provided by `ss`.
"""
function multi_eq(Nspp::Int64 = 10, Niter::Int64 = 10, θ_fc = 0.4,
                  mintotal::Float64 = 0.1, maxtotal::Float64 = 3.0,
                  lengthtotal::Int64 = 10,
                  minP::Int64 = 8, maxP::Int64 = 50, lengthP::Int64 = 10,
                  F::Float64 = 10.0, μ::Float64 = 0.03, n_hts::Int64 = 3, uf::Float64 = 0.1,
                  constant_ss::Bool = false, ss::Float64 = 0.02)

    total_list = collect(range(mintotal, stop = maxtotal, length = lengthtotal));
    P_list = collect(range(minP, stop = maxP, length = lengthP));
    params = reshape(collect(Base.product(total_list, P_list)), (length(total_list) * length(P_list), 1));

    W₀_list = Vector{Float64}(undef, length(params))
    for i in 1:length(params)
        W₀_list[i] = params[i][1] / params[i][2]
    end

    params = hcat(W₀_list, repeat(P_list, inner = length(total_list)),
                  repeat(total_list, outer = length(P_list)));

    params[params[:,1] .> θ_fc, 1] .= θ_fc;
    if constant_ss
        params[:,1] .= params[:,2] .* ss
    end

    sd = []
    for i in 1:Niter
        sd = push!(sd, generate_spp_data(Nspp, 0.9, n_hts, 1.0 / ((15 - 1) / 2), F, μ,
                                         3.0, 0.4, 0.0, 0.0001, 0.00005, 11.0, 0.3, 0.6))
    end

    full_results = Array{Float64}(undef, Nspp*Niter, size(params)[1])
    biomass_results = Array{Float64}(undef, Nspp*Niter, size(params)[1])
    transpir_results = Array{Float64}(undef, Niter, size(params)[1])
    Threads.@threads for i in 1:size(params)[1]
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        biomass_sub = Vector{Float64}(undef, Nspp*Niter)

        for j in 1:Niter
            spp_data = sd[j]
            ## need to recalculate τ, because T changes with params
            spp_data.τ = calc_τ.(spp_data.C₁, spp_data.C₂, F, μ, 1.0 / params[i,2], b)
            eqN = Vector(calc_eqN(spp_data, 1.0 / params[i,2], params[i,1], F, E, θ_fc, μ, false, false, 0.1)[:,:eqN])
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] = eqN
            eqB = calc_eq_biomass(spp_data,
                                  params[i, 1], E, 1.0 / params[i,2], F, μ, uf)
            biomass_sub[[((j-1)*Nspp)+1:1:j*Nspp;]] = eqB
            if sum(eqN .> 0.0) == 0
                transpir_results[j,i] = 0
            else
                transpir_results[j,i] = params[i,1] * params[i, 2]
            end
        end

        println("completed iteration: " * string(i) * " of " * string(size(params)[1]))
        full_results[:,i] = sub_results
        biomass_results[:,i] = biomass_sub
    end

    return [full_results, Niter, params, Nspp, biomass_results, transpir_results]

end


"""
    summarize_multi_eq(multi_eq_output::Vector{Any})

Summarizes the output of `multi_eq`, creating a dataframe of mean equilibrium population densities
across iterations for each precipitation regime. The standard deviation in density is also included
in this dataframe.
"""
function summarize_multi_eq(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3], outer = 4),
                        P = repeat(multi_eq_output[3][:,2], outer = 4),
                        var = repeat(["n", "min", "max", "avg"], inner = Npar),
                        mean = Vector{Float64}(undef, Npar*4),
                        sd = Vector{Float64}(undef, Npar*4))

    summary.total = round.(summary.total .* 100)

    nfeas_temp = Vector{Float64}(undef, Niter); maxfeas_temp = Vector{Float64}(undef, Niter);
    minfeas_temp = Vector{Float64}(undef, Niter); avgfeas_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[1])[2]
        for j in 1:Niter
            tmp = findall(multi_eq_output[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 0.01)
            maxfeas_temp[j] = maximum(tmp, init = 0)
            minfeas_temp[j] = minimum(tmp, init = 0)
            avgfeas_temp[j] = mean(tmp)
            nfeas_temp[j] = length(tmp)
        end
        summary[i, :mean] = mean(nfeas_temp); summary[i, :sd] = std(nfeas_temp)
        summary[Npar+i, :mean] = mean(minfeas_temp); summary[Npar+i, :sd] = std(minfeas_temp)
        summary[Npar*2+i, :mean] = mean(maxfeas_temp); summary[Npar*2+i, :sd] = std(maxfeas_temp)
        summary[Npar*3+i, :mean] = mean(avgfeas_temp); summary[Npar*3+i, :sd] = std(avgfeas_temp)
    end

    return summary

end;


"""
    summarize_multi_eq_biomass(multi_eq_output::Vector{Any})

Summarizes the output of `multi_eq`, creating a dataframe of mean equilibrium biomass
across iterations for each precipitation regime. The standard deviation in biomass is also included
in this dataframe.
"""
function summarize_multi_eq_biomass(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3]),
                        P = repeat(multi_eq_output[3][:,2]),
                        mean = Vector{Float64}(undef, Npar),
                        sd = Vector{Float64}(undef, Npar))

    summary.total = round.(summary.total .* 100)

    la_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[5])[2]
        for j in 1:Niter
            la_temp[j] = sum(multi_eq_output[5][[(Nspp*(j-1))+1:1:Nspp*j;],i])
        end
        summary[i, :mean] = mean(la_temp); summary[i, :sd] = std(la_temp)
    end

    return summary

end;


"""
    summarize_multi_eq_transpiration(multi_eq_output::Vector{Any})

Summarizes the output of `multi_eq`, creating a dataframe of mean evapotranspiration rate
across iterations for each precipitation regime. The standard deviation in evapotranspiration is also included
in this dataframe.
"""
function summarize_multi_eq_transpiration(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3]),
                        P = repeat(multi_eq_output[3][:,2]),
                        mean = Vector{Float64}(undef, Npar),
                        sd = Vector{Float64}(undef, Npar))

    summary.total = round.(summary.total .* 100)

    for i in 1:size(multi_eq_output[6])[2]
        summary[i, :mean] = mean(multi_eq_output[6][:,i]); summary[i, :sd] = std(multi_eq_output[6][:,i])
    end

    return summary

end;

##---------------------------------------------------------------
## Stochastic
##---------------------------------------------------------------

function multi_eq_variable_map(Nspp::Int64 = 10, Niter::Int64 = 10, Nyr::Int64 = 400, θ_fc::Float64 = 0.4,
                               min_map_mean::Float64 = 0.1 * 15, max_map_mean::Float64 = 0.6 * 15,
                               length_map_mean::Int64 = 10,
                               min_map_sd::Float64 = 0.0, max_map_sd::Float64 = 1.5, length_map_sd::Int64 = 5,
                               Pmean::Float64 = 10.0, Pdisp::Float64 = 10.0, n_ht::Int64 = 1,
                               F::Float64 = 10.0, μ::Float64 = 0.03, cluster::Bool = false, b::Float64 = 3.0)

    ## create combinatorial grid of parameter values
    map_mean_list = collect(range(min_map_mean, stop = max_map_mean, length = length_map_mean));
    map_sd_list = collect(range(min_map_sd, stop = max_map_sd, length = length_map_sd));
    pars = hcat(repeat(map_mean_list, inner = length(map_sd_list)),
                repeat(map_sd_list, outer = length(map_mean_list)));

    ## generate communities for simulations (# communities = Niter)
    sd = []
    for i in 1:Niter
        sd = push!(sd, generate_spp_data(Nspp, 0.7, n_ht, 1.0 / Pmean, F, μ, b, 0.4, 0.0, 0.0001, 0.00005, 11.0))
    end

    ## create empty arrays for results
    full_results = Array{Float64}(undef, Nspp*Niter, size(pars)[1])
    leaf_area_results = Array{Float64}(undef, Niter, size(pars)[1])
    transpir_results = Array{Float64}(undef, Niter, size(pars)[1])

    ## loop through parameters and communities in parallel
    Threads.@threads for i in 1:size(pars)[1]

        println("starting iteration: " * string(i) * " of ", string(size(pars)[1]))

        ## create new rainfall regime
        rr = generate_rainfall_regime(Nyr, Pmean, Pdisp, pars[i, 1], pars[i, 2])

        ## create mortality_table in advance to save time
        mort = mortality_table(length(rr[2]), μ, rr[2])

        ## empty vector for parameter set i's results
        sub_results = Vector{Float64}(undef, Nspp*Niter)

        ## loop over communities
        for j in 1:Niter

            println("for community: " * string(j))

            ## conduct simulation
            spp_data = sd[j]
            result = sim_water_ppa_stochastic(spp_data, length(rr[1]), Nspp, 1.0, rr, F, θ_fc, μ,
                                              mort, θ_fc, b, 0.1, false)

            ## equilibrium population density
            eqN = Matrix(result[2])[length(rr[2]), 2:Nspp+1]
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] .= eqN

            ## average leaf area
            la = Matrix(result[1])
            leaf_area_results[j,i] = mean(sum(la[size(la)[1]-
                Int(round(0.1*size(la)[1])):size(la)[1], 2:size(la)[2]], dims = 2))

            ## average annual transpiration
            transpiration = (result[6] - result[5]) ./ rr[2]
            l = length(transpiration)
            transpir_results[j,i] = mean(transpiration[length(transpiration)-
                Int(round(0.1*length(transpiration))):length(transpiration)])

        end

        full_results[:,i] = sub_results
    end

    return [full_results, Niter, pars, Nspp, leaf_area_results, transpir_results]

end




function multi_eq_variable_P(Nspp::Int64 = 10, Niter::Int64 = 10, Nyr::Int64 = 400, θ_fc::Float64 = 0.4,
                             minPmean::Int64 = 0.1 * 15, maxPmean::Int64 = 0.6 * 15,
                             lengthPmean::Int64 = 10,
                             minPdisp::Float64 = 0.0, maxPdisp::Float64 = 1.5, lengthPdisp::Int64 = 5,
                             totalmean::Float64 = 4.0, totalsd::Float64 = 1.5,
                             F::Float64 = 10.0, μ::Float64 = 0.03, cluster::Bool = false)

    Pmean_list = collect(range(minPmean, stop = maxPmean, length = lengthPmean));
    Pdisp_list = collect(range(minPdisp, stop = maxPdisp, length = lengthPdisp));
    pars = hcat(repeat(Pmean_list, inner = length(Pdisp_list)),
                repeat(Pdisp_list, outer = length(Pmean_list)));

    sd = []
    for i in 1:Niter
        sd = push!(sd, generate_spp_data(Nspp, 0.6, 1.0 / mean(Pmean_list), F, μ, 2.5, 0.05, 0.0))
    end

    full_results = Array{Float64}(undef, Nspp*Niter, size(pars)[1])
    leaf_area_results = Array{Float64}(undef, Nspp*Niter, size(pars)[1])
    transpir_results = Array{Float64}(undef, Niter, size(pars)[1])
    Threads.@threads for i in ProgressBar(1:size(pars)[1])
        rr = generate_rainfall_regime(Nyr, Int(round(pars[i, 1])), pars[i, 2], totalmean, totalsd)
        mort = mortality_table(length(rr[2]), μ, rr[2])
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        leaf_area_sub = Vector{Float64}(undef, Nspp*Niter)
        for j in 1:Niter
            spp_data = sd[j]
            eqN = Matrix(sim_water_only_stochastic(spp_data, Nspp, 1.0, rr, θ_fc,
                                                   μ, F, mort, false)[2])[length(rr[2]), 2:Nspp+1]
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] .= eqN
        end
        full_results[:,i] = sub_results
    end

    return [full_results, Niter, pars, Nspp]

end



function summarize_multi_eq_variable_map(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(mapmean = repeat(multi_eq_output[3][:,1], outer = 4),
                        mapsd = repeat(multi_eq_output[3][:,2], outer = 4),
                        var = repeat(["n", "min", "max", "avg"], inner = Npar),
                        mean = Vector{Float64}(undef, Npar*4),
                        sd = Vector{Float64}(undef, Npar*4))

    summary.mapmean = round.(summary.mapmean .* 100)
    summary.mapsd = round.(summary.mapsd .* 100)

    nfeas_temp = Vector{Float64}(undef, Niter); maxfeas_temp = Vector{Float64}(undef, Niter);
    minfeas_temp = Vector{Float64}(undef, Niter); avgfeas_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[1])[2]
        for j in 1:Niter
            tmp = findall(multi_eq_output[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 1e-3)
            maxfeas_temp[j] = maximum(tmp, init = 0)
            minfeas_temp[j] = minimum(tmp, init = 31)
            avgfeas_temp[j] = mean(tmp)
            nfeas_temp[j] = length(tmp)
        end
        summary[i, :mean] = mean(nfeas_temp); summary[i, :sd] = std(nfeas_temp)
        summary[Npar+i, :mean] = mean(minfeas_temp); summary[Npar+i, :sd] = std(minfeas_temp)
        summary[Npar*2+i, :mean] = mean(maxfeas_temp); summary[Npar*2+i, :sd] = std(maxfeas_temp)
        summary[Npar*3+i, :mean] = mean(avgfeas_temp); summary[Npar*3+i, :sd] = std(avgfeas_temp)
    end

    return summary

end;




function plot_multi_eq_variable(data::DataFrame, yvar::String = "n", xvar::Symbol = :W₀,
                                Nspp::Int = 10, save::Bool = true, filename = "")

    if xvar == :W₀
        groupvar = :T
        xl = "Mean W₀"
        yl = "σ W₀"
    else
        groupvar = :W₀
        xl = "Mean T₀ (days)"
        yl = "σ T₀ (days)"
    end

    subdata = data[data.var .== yvar, :]
    subdata = subdata[subdata[:, Symbol(string(groupvar) * "sd")] .== 0.0, :]
    subdata = subdata[.!(isnan.(subdata[:, :mean])), :]

    x1 = Symbol(string(xvar) * "mean")
    x2 = Symbol(string(xvar) * "sd")
    groupvar = Symbol(string(groupvar), "mean")
    subdata[:,groupvar]

    p = plot(subdata[:,x1], subdata[:, x2], subdata.mean, group = subdata[:,groupvar],
             zcolor = wrap(subdata[:,groupvar]),
             st = :surface,
             #surfacecolor = subdata[:,groupvar],
             seriescolor = my_cgrad,
             zlim = [0, Nspp], xlim = [minimum(subdata[:,x1]), maximum(subdata[:,x1])],
             ylim = [minimum(subdata[:,x2]), maximum(subdata[:,x2])],
             xflip = true,
             legend = :topleft, frame = :box,  linewidth = 3, fillalpha = 0.7, colorbar = false,
             xlab = xl, ylab = yl, zlab = "# species persisting")

    subdata
    if save
        savefig(p, filename)
    end

    return p

end;



