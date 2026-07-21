using NetCDF, HDF5, Printf

function _ensure_dims_order(G::Int, M::Int, T::Int)
    # placeholder
    return nothing
end

function save_simulation_netCDF(epi_params::Epidemic_Params, population::Population_Params, output_fname::String; G_coords=String[], M_coords=String[], T_coords=String[])
    
    
    

    try
        M = population.M
        T = epi_params.T


        if length(M_coords) != M
            M_coords = collect(1:M)
        end
        if length(T_coords) != T
            T_coords = collect(1:T) 
        end

        m_dim = NcDim("M", M, atts=Dict("description" => "Region", "Unit" => "unitless"), values=M_coords, unlimited=false)
        t_dim = NcDim("T", T, atts=Dict("description" => "Time", "Unit" => "unitless"), values=T_coords, unlimited=false)
        dimlist = [m_dim, t_dim]

        S  = NcVar("S" , dimlist; atts=Dict("description" => "Suceptibles"), t=Float64, compress=-1)
        I  = NcVar("I" , dimlist; atts=Dict("description" => "Infected"), t=Float64, compress=-1)
        R  = NcVar("R" , dimlist; atts=Dict("description" => "Recovered"), t=Float64, compress=-1)
        
        varlist = [S, I, R,]

        data_dict = Dict()
        #data_dict["S"]  = epi_params.ρˢ 
        #data_dict["I"]  = epi_params.ρᴵ 
        #data_dict["R"]  = epi_params.ρᴿ 
        total_population = sum(population.nᵢ)
        data_dict["S"]  = epi_params.ρˢ .* population.nᵢ / total_population
        data_dict["I"]  = epi_params.ρᴵ .* population.nᵢ / total_population
        data_dict["R"]  = epi_params.ρᴿ .* population.nᵢ / total_population

        isfile(output_fname) && rm(output_fname)

        NetCDF.create(filename, varlist, mode=NC_NETCDF4)
        for (var_label, data) in data_dict
            ncwrite(data, filename, var_label)
        end
    catch e
        @error "Error saving simulation output" exception=(e, catch_backtrace())
        rethrow(e)
    end
    @debug "- Done saving"
end

function save_simulation_hdf5(epi_params::Epidemic_Params, population::Population_Params, output_fname::String; export_time_t=-1)
    h5open(output_fname, "w") do f
        write(f, "S", epi_params.ρˢᵍ)
        write(f, "I", epi_params.ρᴵᵍ)
        write(f, "R", epi_params.ρᴿᵍ)
        write(f, "arrival_times", epi_params.arrival_times)
    end
    return nothing
end

function save_observables_netCDF(epi_params::Epidemic_Params, population::Population_Params, output_fname::String; G_coords=String[], M_coords=String[], T_coords=String[])
    # Simple observables: S_total, I_total, R_total, arrival_times
    T = epi_params.T
    S_tot = zeros(Float64, T)
    I_tot = zeros(Float64, T)
    R_tot = zeros(Float64, T)
    for t in 1:T
        S_tot[t] = sum(epi_params.ρˢᵍ[:, :, t] .* population.nᵢᵍ) / sum(population.nᵢ)
        I_tot[t] = sum(epi_params.ρᴵᵍ[:, :, t] .* population.nᵢᵍ) / sum(population.nᵢ)
        R_tot[t] = sum(epi_params.ρᴿᵍ[:, :, t] .* population.nᵢᵍ) / sum(population.nᵢ)
    end

    nc = NetCDF.create(output_fname)
    try
        defDim(nc, "T", T)
        defVar(nc, "S_total", S_tot)
        defVar(nc, "I_total", I_tot)
        defVar(nc, "R_total", R_tot)
        defVar(nc, "arrival_times", epi_params.arrival_times)
    finally
        NetCDF.close(nc)
    end
    return nothing
end
