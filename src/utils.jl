using DataFrames, CSV



struct Epidemic_Params
    T::Int
    dt::Float64
    dt_save::Float64
    μ::Float64
    R₀::Float64
    β::Float64
    ρˢ::Array{Float64,2}   # M x T
    ρᴵ::Array{Float64,2}
    ρᴿ::Array{Float64,2}
    arrival_times::Vector{Float64}
end


"""
Initialize epidemic parameters struct. T is number of saved timesteps.
"""
function init_epi_parameters_struct(G::Int,M::Int, T::Int, epi_params_dict::AbstractDict)
    μ = Float64.(epi_params_dict["μ"])
    R₀ = Float64.(epi_params_dict["R₀"])
    β = μ .* R₀ 
    dt = Float64(epi_params_dict["dt"])
    dt_save = Float64(epi_params_dict["dt_save"])

    ρˢ = zeros(Float64, M, T)
    ρᴵ = zeros(Float64, M, T)
    ρᴿ = zeros(Float64, M, T)
    arrival_times = ones(Float64, M) .* (T)

    epi = Epidemic_Params(T, dt, dt_save, μ, R₀, β, ρˢ, ρᴵ, ρᴿ, arrival_times)
    return epi
end

struct Population_Params
    M::Int
    nᵢ::Vector{Float64}        # total per patch (length M)
    nᵢᵍ::Array{Float64,2}      # per strata (G x M)
    mobility::Array{Float64,2} # M x M
    mobility_baseline::Union{Nothing, Array{Float64,2}}
    commuter_population::Array{Int,2} # M x M
    mobility_normalized::Array{Float64,2}
end

"""
Initialize Population_Params from provided dicts and DataFrames.
Expects metapop_df to contain at least a column with population counts
and network_df either mobility matrix or edgelist information.
"""
function init_pop_param_struct(G::Int, M::Int, G_coords::Array{String, 1}, pop_params_dict::AbstractDict, 
                               metapop_df::DataFrame, network_df::DataFrame)
    # Extract patch labels and populations
    #This has to be done in EpiSim.jl
   
    # Total column is always required
    if !(:"total" in names(metapop_df))
        error("metapopulation DataFrame must contain a 'total' column")
    end
    nᵢ = Float64.(metapop_df[:, "total"])

    # Subpopulation by strata
    nᵢᵍ = copy(transpose(Array{Float64,2}(metapop_df[:, G_coords])))
  
    
    # Network: try to read a mobility matrix from network_df
    # network_df may be a square matrix DataFrame or an edgelist with columns source,target,weight
    mobility = zeros(Float64, M, M)

    #Checking if network_df is a square matric MxM
    if size(network_df, 1) == M && size(network_df, 2) == M
        # assume this is a mobility matrix
        for i in 1:M
            for j in 1:M
                mobility[i, j] = parse(Float64, string(network_df[i, j]))
            end
        end
    else
        # assume edgelist with columns source,target,weight
        src_col = hasproperty(network_df, :source_idx) ? :source_idx : (hasproperty(network_df, :i) ? :i : names(network_df)[1])
        tgt_col = hasproperty(network_df, :target_idx) ? :target_idx : (hasproperty(network_df, :j) ? :j : names(network_df)[2])
        w_col   = hasproperty(network_df, :weight) ? :weight : (hasproperty(network_df, :w) ? :w : names(network_df)[3])
        for row in eachrow(network_df)
            i = Int(row[src_col])
            j = Int(row[tgt_col])
            mobility[i, j] = Float64(row[w_col])
        end
    end

    commuter_population, mobility_normalized = compute_commuter_population(mobility, nᵢ, M)

    pop = Population_Params(M, nᵢ, nᵢᵍ, mobility, nothing, commuter_population, mobility_normalized)
    return pop
end

"""
Set initial compartments into epi_params arrays. initial_compartments_dict should have keys "S","I","R" with arrays of shape (G,M).
If scale_by_population true, treat values as absolute counts and convert to fractions.
"""
function set_compartments!(epi_params::Epidemic_Params, population::Population_Params, npi_params, initial_compartments_dict::Dict; scale_by_population::Bool=true)
    M = population.M
    # Expect arrays or scalars
    S0 = haskey(initial_compartments_dict, "S") ? initial_compartments_dict["S"] : population.nᵢ
    I0 = haskey(initial_compartments_dict, "I") ? initial_compartments_dict["I"] : zeros(Float64, M)
    R0 = haskey(initial_compartments_dict, "R") ? initial_compartments_dict["R"] : zeros(Float64, M)

    for i in 1:M
        if scale_by_population
            denom = population.nᵢ[i] > 0 ? population.nᵢ[i] : 1.0
            epi_params.ρˢ[i, 1] = S0[i] / denom
            epi_params.ρᴵ[i, 1] = I0[i] / denom
            epi_params.ρᴿ[i, 1] = R0[i] / denom
        else
            epi_params.ρˢ[i, 1] = S0[i]
            epi_params.ρᴵ[i, 1] = I0[i]
            epi_params.ρᴿ[i, 1] = R0[i]
        end
    end

    # Build internal commuter counts per strata: commuter_population is G x M x M
    # Distribute susceptibles/infecteds across commuter compartments proportionally to commuter population
    # store as integer arrays for stochastic updates
    # We will store these in local refs inside the run function when called
    return nothing
end

"""
Create an initial_compartments_dict from seeds/conditions. A simple helper that
builds zero arrays with a seed I0 in a random patch or in patches_idxs if provided.
"""
function create_initial_compartments_dict(M_coords::Array{String,1}, G_coords::Array{String,1}, nᵢᵍ::Array{Float64,2},
                                          conditions₀, patches_idxs; scale_seeds = nothing)
    M = length(M_coords)
    G = length(G_coords)
    NDIMS = length((G, M))
    
    comp_coords = ["S", "I", "R"]
    
    if scale_seeds !== nothing
        @info "- Scaling initial seeds by factor $(scale_seeds)"
        conditions₀ .= conditions₀ .* scale_seeds
    end

    @debug "- Creating initial compartment dict using arrays of size ($(G), $(M))"
    init_compartments_dict = Dict{String, Array{Float64, NDIMS}}(label => zeros(G, M) for label in comp_coords)
    init_compartments_dict["I"][:, patches_idxs] .= conditions₀
    init_compartments_dict["S"][:,:] .= nᵢᵍ - init_compartments_dict["I"][:, :]
    @info "- Setting remaining population $(sum(init_compartments_dict["S"])) in compartment S" 
    
    return init_compartments_dict
end

"""
Given mobility matrix and per-stratum population nᵢᵍ (G x M), compute commuter
population per strata: commuter_population[g, i, j] = mobility_normalized[i,j] * nᵢᵍ[g, i]
Returns Array{Int,3} of shape (G,M,M)
"""
function compute_commuter_population(mobility::Array{Float64,2}, nᵢ::Array{Float64,1}, M::Int)
    mobility_subpops = mobility * ones(M)
    mobility_normalized = copy(mobility)
    for i in 1:M
        s = mobility_subpops[i]
        if s > 0
            mobility_normalized[i, :] .= mobility[i, :] ./ s
        else
            mobility_normalized[i, :] .= 0.0
        end
    end
    commuter = zeros(Int, M, M)
    
    for i in 1:M
        for j in 1:M
            val = round(Int, mobility_normalized[i,j] * nᵢ[i])
            commuter[i, j] = max(val, 0)
        end
    end
    
    return commuter, mobility_normalized
end
