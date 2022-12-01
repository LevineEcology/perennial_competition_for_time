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
function multi_eq_constant_water(Nspp::Int64 = 10, Niter::Int64 = 10,
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
            sub_results[[((j-1)*Nspp)+1:1:j*Nspp;]] = calc_eqN(spp_data, F, E, params[i][1])[:,:eqN]
        end
        full_results[:,i] = sub_results
    end

    return [full_results, Niter, params, Nspp]

end

function summarize_multi_eq(multo_eq_output)

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
