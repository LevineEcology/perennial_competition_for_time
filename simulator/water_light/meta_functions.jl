##---------------------------------------------------------------
## WATER_AND_LIGHT_SIMULATOR -- meta_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## November 2022
##---------------------------------------------------------------

function multi_eq(Nspp::Int64 = 10, Niter::Int64 = 10, Nyr::Int64 = 100, θ_fc = 0.4,
                  mintotal::Float64 = 0.1 * 15, maxtotal::Float64 = 0.6 * 15,
                  lengthtotal::Int64 = 10,
                  minP::Int64 = 8, maxP::Int64 = 50, lengthP::Int64 = 10,
                  F::Float64 = 10.0, μ::Float64 = 0.03)

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
    params

    sd = []
    for i in 1:Niter
        sd = push!(sd, generate_spp_data(Nspp, 0.6, 1.0 / ((maxP - minP) / 2), F, μ,
                                         2.5, 0.5, 0.0, 0.0001, 0.00005))
    end

    full_results = Array{Float64}(undef, Nspp*Niter, size(params)[1])
    leaf_area_results = Array{Float64}(undef, Nspp*Niter, size(params)[1])
    transpir_results = Array{Float64}(undef, Niter, size(params)[1])
    Threads.@threads for i in 1:size(params)[1]
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        leaf_area_sub = Vector{Float64}(undef, Nspp*Niter)
        for j in 1:Niter
            spp_data = sd[j]
            eqN = Matrix(sim_water_ppa(spp_data, Nyr, nrow(spp_data), 1.0, μ, F,
                                    Int(round(params[i,2])), params[i,1], θ_fc, zeros(1,1),
                                       false)[2])[Int(Nyr*round(params[i,2])),2:nrow(spp_data)+1]
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] = eqN
            leaf_area_sub[[((j-1)*Nspp)+1:1:j*Nspp;]] = calc_eq_leaf_area(eqN, F, μ)
            if sum(eqN .> 0.0) == 0
                transpir_results[j,i] = 0
            else
                transpir_results[j,i] = (params[i,1] - minimum(spp_data[eqN .> 0.0, :Wᵢ])) * params[i, 2]
            end
        end
        println("completed iteration: " * string(i) * " of ", * string(size(params)[1]))
        full_results[:,i] = sub_results
        leaf_area_results[:,i] = leaf_area_sub
    end

    return [full_results, Niter, params, Nspp, leaf_area_results, transpir_results]

end


function summarize_multi_eq(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3], outer = 4),
                        T = repeat(multi_eq_output[3][:,2], outer = 4),
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


function summarize_multi_eq_leafarea(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3]),
                        T = repeat(multi_eq_output[3][:,2]),
                        mean = Vector{Float64}(undef, Npar),
                        sd = Vector{Float64}(undef, Npar))

    la_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[5])[2]
        for j in 1:Niter
            la_temp[j] = sum(multi_eq_output[5][[(Nspp*(j-1))+1:1:Nspp*j;],i])
        end
        summary[i, :mean] = mean(la_temp); summary[i, :sd] = std(la_temp)
    end

    return summary

end;

function summarize_multi_eq_transpiration(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(total = repeat(multi_eq_output[3][:,3]),
                        T = repeat(multi_eq_output[3][:,2]),
                        mean = Vector{Float64}(undef, Npar),
                        sd = Vector{Float64}(undef, Npar))

    for i in 1:size(multi_eq_output[6])[2]
        summary[i, :mean] = mean(multi_eq_output[6][:,i]); summary[i, :sd] = std(multi_eq_output[6][:,i])
    end

    return summary

end;



## NEED TO EDIT
function plot_multi_eq(data::DataFrame, yvar::String = "n", xvar::Symbol = :W₀, Nspp::Int = 10, save::Bool = true, filename = "")

    if xvar == :W₀
        groupvar = :T
        xl = "W₀"
        ll = "T₀"
    else
        groupvar = :W₀
        xl = "T₀"
        ll = "W₀"
    end

    if yvar == "n"
        yl = "# species coexisting"
    elseif yvar == "min"
        yl = "Min. ID of coexisting species"
    else
        yl = "Avg. ID of coexisting species"
    end

    subdata = data[data.var .== yvar, :]
    p = plot(subdata[:,xvar], subdata.mean, group = subdata[:,groupvar], line_z = subdata[:,groupvar],
             fill_z = subdata[:,groupvar], ribbon = subdata.sd, seriescolor = my_cgrad,
             seriestype = :line, ylim = [0, Nspp], xlim = [minimum(subdata[:,xvar]), maximum(subdata[:,xvar])],
             legend = :topleft, frame = :box, grid = false, linewidth = 3, fillalpha = 0.3,
             xlab = xl, legendtitle = ll, ylab = yl, colorbar = false)

    if save
        savefig(p, filename)
    end

    return p

end;


##---------------------------------------------------------------
## Stochastic
##---------------------------------------------------------------

function multi_eq_variable_total(Nspp::Int64 = 10, Niter::Int64 = 10, Nyr::Int64 = 400, θ_fc::Float64 = 0.4,
                                 mintotalmean::Float64 = 0.1 * 15, maxtotalmean::Float64 = 0.6 * 15,
                                 lengthtotalmean::Int64 = 10,
                                 mintotalsd::Float64 = 0.0, maxtotalsd::Float64 = 1.5, lengthtotalsd::Int64 = 5,
                                 Pmean::Int64 = 10, Pdisp::Float64 = 10.0,
                                 F::Float64 = 10.0, μ::Float64 = 0.03, cluster::Bool = false)

    totalmean_list = collect(range(mintotalmean, stop = maxtotalmean, length = lengthtotalmean));
    totalsd_list = collect(range(mintotalsd, stop = maxtotalsd, length = lengthtotalsd));
    pars = hcat(repeat(totalmean_list, inner = length(totalsd_list)),
                repeat(totalsd_list, outer = length(totalmean_list)));

    sd = []
    for i in 1:Niter
        sd = push!(sd, generate_spp_data(Nspp, 0.6, 1.0 / Pmean, F, μ, 2.5, 0.05, 0.0))
    end

    full_results = Array{Float64}(undef, Nspp*Niter, size(pars)[1])
    leaf_area_results = Array{Float64}(undef, Niter, size(pars)[1])
    transpir_results = Array{Float64}(undef, Niter, size(pars)[1])
    Threads.@threads for i in ProgressBar(1:size(pars)[1])
        rr = generate_rainfall_regime(Nyr, Pmean, Pdisp, pars[i, 1], pars[i, 2])
        mort = mortality_table(length(rr[2]), μ, rr[2])
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        for j in 1:Niter
            spp_data = sd[j]
            result = sim_water_only_stochastic(spp_data, Nspp, 1.0, rr, θ_fc, μ,
                                                   F, mort, false)
            ## equilibrium population density
            eqN = Matrix(result[2])[length(rr[2]), 2:Nspp+1]
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] .= eqN

            ## average leaf area
            la = Matrix(result[1])
            leaf_area_results[j,i] = mean(sum(la[size(la)[1]-
                Int(round(0.1*size(la)[1])):size(la)[1],2:size(la)[2]], dims = 2))

            ## average annual transpiration
            transpiration = (result[6] - result[5]) ./ rr[2]
            l = length(transpiration)
            transpiration_results = mean(transpiration[length(transpiration)-
                Int(round(0.1*length(transpiration))):length(transpiration)])

        end
        full_results[:,i] = sub_results
    end

    return [full_results, Niter, pars, Nspp, leaf_area_results, transpiration_results]

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



function summarize_multi_eq_variable_total(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(totalmean = repeat(multi_eq_output[3][:,1], outer = 4),
                        totalsd = repeat(multi_eq_output[3][:,2], outer = 4),
                        var = repeat(["n", "min", "max", "avg"], inner = Npar),
                        mean = Vector{Float64}(undef, Npar*4),
                        sd = Vector{Float64}(undef, Npar*4))

    summary.totalmean = round.(summary.totalmean .* 100)
    summary.totalsd = round.(summary.totalsd .* 100)

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




"""
    dens_est_constant_water(Nspp::Int64 = 10, Niter::Int64 = 10,
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)

TBW
"""
function multi_eq_geography(sm_grid, Nspp::Int64 = 40)

    spp_data = generate_spp_data(Nspp, 0.8, 40.0, 100.0, 0.1, 2.5, 0.4, 0.0) ## initialize spp_data
    diversity_results = Matrix{Int64}(undef, size(sm_grid[:,:,1], 1), size(sm_grid[:,:,1], 2))
    biomass_results = Matrix{Float64}(undef, size(sm_grid[:,:,1], 1), size(sm_grid[:,:,1], 2))
    phenology_results = Matrix{Float64}(undef, size(sm_grid[:,:,1], 1), size(sm_grid[:,:,1], 2))
    @Threads.threads for i in 1:size(sm_grid[:,:,1], 2)
        for j in 1:size(sm_grid[:,:,1], 1)
            result = calc_eqN(spp_data, F, E, Float64(sm_grid[:,:,1][j,i]))
            winners = result[result.eqN .> 1e-20, :spp]
            diversity_results[j,i] = count(result.eqN .> 1e-20)
            biomass_results[j,i] = sum(eq_biomass(result[result.eqN .> 1e-20,:], result[result.eqN .> 1e-20,:eqN]))
            phenology_results[j,i] = mean(winners)
        end
    end

    return [diversity_results, biomass_results, phenology_results]

end

function eq_biomass(spp_data::DataFrame, eqN::Vector{Float64}, T = 40.0, b = 2.5, μ = 0.1)

    n_data = DataFrame(rowkey = repeat([1:1:1e3;], inner = length(eqN)),
                       spp = repeat(string.([1:1:length(eqN);]), outer = Int(1e3)),
                       N = repeat(repeat(Float64[0], inner = length(eqN)), inner = Int(1e3)));
    n_data = unstack(n_data, :spp, :N);
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]]);

    biomass_data = DataFrame(rowkey = repeat([1:1:1e3;], inner = length(eqN)),
                       spp = repeat(string.([1:1:length(eqN);]), outer = Int(1e3)),
                             B = repeat(repeat(Float64[0], inner = length(eqN)), inner = Int(1e3)));
    biomass_data = unstack(biomass_data, :spp, :B);
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]]);
    g = calc_g.(spp_data.τ, spp_data[:,:C₁], spp_data[:,:C₂], 40.0)

    for i in 1:length(eqN)
        for j in 1:Int(1e3)
            n_data[j,i+1] = (1-μ)^j * eqN[i]
            biomass_data[j,i+1] = g[i]^b * (j*T)^b
        end
    end

    ΣB(biomass_data, n_data, Vector{Float64}(undef, length(eqN)))

end
