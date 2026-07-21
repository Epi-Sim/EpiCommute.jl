using Distributions, Random, Printf
using DataFrames




function run_epidemic_spreading!(epi_params::Epidemic_Params, population::Population_Params, npi_params::NPI_Params; verbose::Bool=false, rng = MersenneTwister())
    M = population.M
    T = epi_params.T
    dt = epi_params.dt
    dt_save = epi_params.dt_save
    μ = epi_params.μ
    β = epi_params.β

    # Initialize commuter counts from epi_params initial fractions
    # commuter_population: M x M
    commuter = copy(population.commuter_population)

    # Initialize S, I, R as integer arrays G x M x M
    S = zeros(Int, M, M)
    I = zeros(Int, M, M)
    R = zeros(Int, M, M)

    if npi_params.quarantine_mode == "isolation"
        S = round(Int, commuter .* kappa)
        R = round(Int, commuter .* (1 .- kappa))
    else
        S = commuter
    end
    
    # Initialize the infection seed

    for i in 1:M
        seed_i = round(Int, epi_params.ρᴵ[i, 1] * population.nᵢ[i])
        for n in 1:seed_i
            j = rand(1:M)
            while I[i, j] == commuter[i, j]
                j = rand(1:M)
            end
            if I[i, j] < commuter[i, j]
                I[i, j] += 1
                S[i, j] -= 1
            end
        end
    end

    for t in 1:(T-1)

        t_count = 0.0
        while t_count < dt_save        
            
            # Home work of infection
            I_ij_sumj = vec(sum(I, dims=2))
            N_ij_sumj = vec(sum(S, dims=2) .+ sum(I, dims=2) .+ sum(R, dims=2))
            lambda_home = 0.5 * β .* (I_ij_sumj ./ N_ij_sumj)
            # Work force of infection
            I_ji_sumj = vec(sum(I, dims=1))
            N_ji_sumj = vec(sum(S, dims=1) .+ sum(I, dims=1) .+ sum(R, dims=1))
            lambda_work = 0.5 * β .* (I_ji_sumj ./ N_ji_sumj)

            for i in 1:M
                # Normal infection rate
                if npi_params == nothing || npi_params.quarantine_mode == "isolation" || npi_params.quarantine_mode == nothing
                    lambda_home_eff = lambda_home[i]
                    lambda_work_eff = copy(lambda_work)
                elseif npi_params.quarantine_mode == "distancing"
                    lambda_home_eff = kappa[i, :] .* lambda_home[i] 
                    lambda_work_eff = kappa[i, :] .* lambda_work_eff
                end   
                # Calculate infections
                # Home force of infection
                p_inf = 1.0 .- exp.(-(lambda_home_eff .+ lambda_work_eff) .* dt )
                dSI_i = rand.(Binomial.(S[i,:], p_inf ))
                # Calculate recoveries
                p_rec = 1.0 .- exp.(-μ * dt)
                dIR_i = rand.(Binomial.(I[i,:], p_rec))

                # Update system
                S[i, :] = S[i, :] .- dSI_i
                I[i, :] = I[i, :] .+ dSI_i .- dIR_i
                R[i, :] = R[i, :] .+ dIR_i
                
            end            

            t_count += dt
        end

        epi_params.ρˢ[:, t + 1] = sum(S, dims=2)[:, 1] ./ population.nᵢ
        epi_params.ρᴵ[:, t + 1] = sum(I, dims=2)[:, 1] ./ population.nᵢ
        epi_params.ρᴿ[:, t + 1] = sum(R, dims=2)[:, 1] ./ population.nᵢ
    end

    # # Initialize commuter counts from epi_params initial fractions
    # # commuter_population: M x M
    # commuter = copy(population.commuter_population)

    # # Initialize S_counts, I_counts, R_counts as integer arrays G x M x M
    # S_counts = zeros(Int, M, M)
    # I_counts = zeros(Int, M, M)
    # R_counts = zeros(Int, M, M)

    
    # for i in 1:M
    #     # fraction susceptible in patch i at t=1
    #     s_frac = epi_params.ρˢ[i, 1]
    #     i_frac = epi_params.ρᴵ[i, 1]
    #     r_frac = epi_params.ρᴿ[i, 1]
    #     total = sum(commuter[ i, :])
    #     if total == 0
    #         continue
    #     end
    #     for j in 1:M
    #         # distribute counts proportionally to commuter population for origin i
    #         proportion = commuter[i, j] / max(total, 1)
    #         S_counts[i, j] = round(Int, proportion * s_frac * population.nᵢ[i])
    #         I_counts[i, j] = round(Int, proportion * i_frac * population.nᵢ[i])
    #         R_counts[i, j] = round(Int, proportion * r_frac * population.nᵢ[i])
    #     end
    # end
    

    # # Prepare arrival time threshold
    # I_threshold = 0.001

    # # kappa factor (if distancing)
    # kappa = npi_params == nothing ? ones(M, M) : ones(M, M)
    # if npi_params != nothing && npi_params.quarantine_mode == "distancing"
    #     # Build kappa matrix from mobility and baseline if available: use mobility_normalized ratio if baseline present
    #     if population.mobility_baseline != nothing
    #         base = population.mobility_baseline
    #         kappa_mat = copy(population.mobility_normalized)
    #         for i in 1:M
    #             for j in 1:M
    #                 denom = sum(base[i, :])
    #                 kappa_mat[i, j] = denom > 0 ? population.mobility[i,j] / base[i,j] : 1.0
    #             end
    #         end
    #         kappa = kappa_mat
    #     else
    #         kappa = ones(M, M)
    #     end
    # end


    # # Time loop
    # for t in 1:(T-1)
    #     # Aggregate infected and populations by origin/destination
    #     # sum over destination j for origin i
    #     sum_I_origin = zeros(Float64, M)
    #     sum_S_origin = zeros(Float64, M)
    #     sum_R_origin = zeros(Float64, M)
    #     for i in 1:M
    #         sum_I_origin[i] += sum(I_counts[i, :])
    #         sum_S_origin[i] += sum(S_counts[i, :])
    #         sum_R_origin[i] += sum(R_counts[i, :])
    #     end
    #     # sum over origin i for destination j
    #     sum_I_dest = zeros(Float64, M)
    #     sum_S_dest = zeros(Float64, M)
    #     sum_R_dest = zeros(Float64, M)
    #     for j in 1:M
    #         sum_I_dest[j] += sum(I_counts[:, j])
    #         sum_S_dest[j] += sum(S_counts[:, j])
    #         sum_R_dest[j] += sum(R_counts[:, j])
    #     end

    #     N_origin = sum_S_origin .+ sum_I_origin .+ sum_R_origin
    #     N_dest = sum_S_dest .+ sum_I_dest .+ sum_R_dest
    #     # Avoid zeros
    #     for i in 1:M
    #         if N_origin[i] == 0
    #             N_origin[i] = 1e-7
    #         end
    #         if N_dest[i] == 0
    #             N_dest[i] = 1e-7
    #         end
    #     end

    #     lambda_home = 0.5 * β .* (sum_I_origin ./ N_origin)
    #     lambda_work = 0.5 * β .* (sum_I_dest ./ N_dest)
     

    #     # For each origin i, update each commuter compartment j
    #     for i in 1:M
    #         # Compute lambda_home_eff and lambda_work_eff arrays of length M
    #         if npi_params == nothing || npi_params.quarantine_mode == "isolation" || npi_params.quarantine_mode == nothing
    #             lambda_home_eff_scalar = lambda_home[i]
    #             lambda_work_eff_vec = copy(lambda_work)
    #             for j in 1:M
    #                 p_inf = 1.0 - exp(- (lambda_home_eff_scalar + lambda_work_eff_vec[j]) * dt)
    #                 p_inf = clamp(p_inf, 0.0, 1.0)
    #                 nS = S_counts[i, j]
    #                 if nS > 0
    #                     dSI = rand(rng, Binomial(nS, p_inf))
    #                 else
    #                     dSI = 0
    #                 end
    #                 nI = I_counts[i, j]
    #                 p_rec = 1.0 - exp(- μ * dt)
    #                 dIR = nI > 0 ? rand(rng, Binomial(nI, p_rec)) : 0

    #                 S_counts[i, j] -= dSI
    #                 I_counts[i, j] += dSI - dIR
    #                 R_counts[i, j] += dIR
    #             end
    #         elseif npi_params.quarantine_mode == "distancing"
    #             # lambda_home_eff is vector kappa[i, :] * lambda_home[i]
    #             lambda_home_eff_vec = lambda_home[i] .* kappa[i, :]
    #             # lambda_work_eff is vector kappa[:, i] .* lambda_work
    #             lambda_work_eff_vec = lambda_work .* view(kappa, :, i)
    #             for j in 1:M
    #                 p_inf = 1.0 - exp(- (lambda_home_eff_vec[j] + lambda_work_eff_vec[j]) * dt)
    #                 p_inf = clamp(p_inf, 0.0, 1.0)
    #                 nS = S_counts[i, j]
    #                 dSI = nS > 0 ? rand(rng, Binomial(nS, p_inf)) : 0
    #                 nI = I_counts[i, j]
    #                 p_rec = 1.0 - exp(- μ * dt)
    #                 dIR = nI > 0 ? rand(rng, Binomial(nI, p_rec)) : 0

    #                 S_counts[i, j] -= dSI
    #                 I_counts[i, j] += dSI - dIR
    #                 R_counts[i, j] += dIR
    #             end
    #         end
    #     end

    #     # After updates, store fractions into epi_params for t+1
    #     for i in 1:M
    #         tot_comm = sum(commuter[i, :])
    #         denom = population.nᵢ[i] > 0 ? population.nᵢ[i] : max(tot_comm, 1)
    #         s_sum = sum(S_counts[i, :])
    #         i_sum = sum(I_counts[i, :])
    #         r_sum = sum(R_counts[i, :])
    #         epi_params.ρˢ[i, t+1] = s_sum / denom
    #         epi_params.ρᴵ[i, t+1] = i_sum / denom
    #         epi_params.ρᴿ[i, t+1] = r_sum / denom

    #         # Update arrival times
    #         if epi_params.arrival_times[i] > T && (i_sum / max(denom,1) > 0.001)
    #             epi_params.arrival_times[i] = t+1
    #         end
    #     end

    #     if verbose
    #         @printf("Time step %d completed\n", t)
    #     end
    # end

    return nothing
end
