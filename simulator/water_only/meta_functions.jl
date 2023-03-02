##---------------------------------------------------------------
## WATER_ONLY_SIMULATOR -- meta_functions.jl
##
## By: Jacob Levine -- jacoblevine@princeton.edu
## November 2022
##---------------------------------------------------------------

"""
    dens_est_constant_water(Nspp::Int64 = 10, Niter::Int64 = 10,
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)

TBW
"""
function multi_eq_constant_water(Nspp::Int64 = 10, Niter::Int64 = 10, θ_fc = 0.4,
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)

    W₀_list = collect(range(minW₀, stop = maxW₀, length = lengthW₀));
    T_list = collect(range(minT, stop = maxT, length = lengthT));
    params = collect(Base.product(W₀_list, T_list));

    spp_data = generate_spp_data(Nspp, 0.6) ## initialize spp_data
    full_results = Array{Float64}(undef, Nspp*Niter, length(params))
    Threads.@threads for i in 1:length(params)
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        for j in 1:Niter
            spp_data = generate_spp_data(Nspp, 0.8, params[i][2], 100.0, 0.1, 2.5, 0.4, 0.0)
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] = calc_eqN(spp_data, F, E, params[i][1], θ_fc)[:,:eqN]
        end
        full_results[:,i] = sub_results
    end

    return [full_results, Niter, params, Nspp]

end

function summarize_multi_eq(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(W₀ = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:2][:,1], outer = 4),
                        T = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:2][:,2], outer = 4),
                        var = repeat(["n", "min", "max", "avg"], inner = Npar),
                        mean = Vector{Float64}(undef, Npar*4),
                        sd = Vector{Float64}(undef, Npar*4))

    nfeas_temp = Vector{Float64}(undef, Niter); maxfeas_temp = Vector{Float64}(undef, Niter);
    minfeas_temp = Vector{Float64}(undef, Niter); avgfeas_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[1])[2]
        for j in 1:Niter
            tmp = findall(multi_eq_output[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 0.0)
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

"""
    multi_eq_variable_water(Nspp::Int64 = 10, Niter::Int64 = 10,
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)
                                 minW₀::Float64 = 0.1, maxW₀::Float64 = 0.6, lengthW₀::Int64 = 10,
                                 minT::Float64 = 1.0, maxT::Float64 = 100.0, lengthT::Int64 = 10)

TBW
"""
function multi_eq_variable_water(Nspp::Int64 = 10, Niter::Int64 = 10,
                                 Nyr_sim = 4000,
                                 minW₀mean::Float64 = 0.1, maxW₀mean::Float64 = 0.6,
                                 lengthW₀mean::Int64 = 10,
                                 minW₀sd::Float64 = 0.0, maxW₀sd::Float64 = 0.2,
                                 lengthW₀sd::Int64 = 10,
                                 minTmean::Float64 = 1.0, maxTmean::Float64 = 100.0,
                                 lengthTmean::Int64 = 10,
                                 minTsd::Float64 = 0.0, maxTsd::Float64 = 30.0,
                                 lengthTsd::Int64 = 10,
                                 θₘ = 0.4)

    W₀mean_list = collect(range(minW₀mean, stop = maxW₀mean, length = lengthW₀mean));
    W₀sd_list = collect(range(minW₀sd, stop = maxW₀sd, length = lengthW₀sd));
    Tmean_list = collect(range(minTmean, stop = maxTmean, length = lengthTmean));
    Tsd_list = collect(range(minTsd, stop = maxTsd, length = lengthTsd));
    params = vcat(collect(Base.product(W₀mean_list, W₀sd_list, Tmean_list, Tsd_list)))

    spp_data = generate_spp_data(Nspp) ## initialize spp_data
    full_results = Array{Float64}(undef, Nspp*Niter, length(params))
    Threads.@threads for i in 1:length(params)
        sub_results = Vector{Float64}(undef, Nspp*Niter)
        sim = Vector{Float64}(undef, Nspp)
        for j in 1:Niter
            spp_data = generate_spp_data(Nspp, 0.8, params[i][3], 100.0, 0.1, 2.5, 0.4, 0.0)
            sim .= Vector(sim_water_only_stochastic(spp_data, Nyr_sim, nrow(spp_data), 1.0, true, true,
                                                    params[i][1], params[i][2], params[i][3], params[i][4])[2][Nyr_sim, 2:Nspp+1], θₘ)
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] .= replace(x -> isless(x, 1e-15) ? 0 : x, sim)
        end
        full_results[:,i] = sub_results
    end

    summarize_multi_eq_variable([full_results, Niter, params, Nspp])

end

function summarize_multi_eq_variable(multi_eq_output::Vector{Any})

    ## initialize data
    Niter = multi_eq_output[2]; Nspp = multi_eq_output[4]; tmp = Vector{Int64}; Npar = size(multi_eq_output[1])[2];
    summary = DataFrame(W₀mean = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:4][:,1], outer = 4),
                        W₀sd = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:4][:,2], outer = 4),
                        Tmean = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:4][:,3], outer = 4),
                        Tsd = repeat([multi_eq_output[3][i][j] for i = 1:Npar, j = 1:4][:,4], outer = 4),
                        var = repeat(["n", "min", "max", "avg"], inner = Npar),
                        mean = Vector{Float64}(undef, Npar*4),
                        sd = Vector{Float64}(undef, Npar*4))

    nfeas_temp = Vector{Float64}(undef, Niter); maxfeas_temp = Vector{Float64}(undef, Niter);
    minfeas_temp = Vector{Float64}(undef, Niter); avgfeas_temp = Vector{Float64}(undef, Niter);

    for i in 1:size(multi_eq_output[1])[2]
        for j in 1:Niter
            tmp = findall(multi_eq_output[1][[(Nspp*(j-1))+1:1:Nspp*j;],i] .> 0.0)
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
