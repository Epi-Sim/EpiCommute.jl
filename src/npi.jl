using Dates, CSV, DataFrames

struct NPI_Params
    quarantine_mode::Union{Nothing, String}
end

function init_NPI_parameters_struct(data_path::String, npi_params_dict::AbstractDict, kappa0_filename::Union{String,Nothing}, first_day::Date)
    # Basic parsing: expect either vectors in dict or single values
    if haskey(npi_params_dict, "are_there_npi") && npi_params_dict["are_there_npi"] == false
        return NPI_Params(nothing)
    end

    # κ₀s = get(npi_params_dict, "κ₀s", [0.0])
    # tcs = get(npi_params_dict, "tᶜs", [typemax(Int)])
    qm = get(npi_params_dict, "quarantine_mode", nothing)

    # # If kappa0_filename provided, try to read it (not required)
    # if kappa0_filename !== nothing
    #     kappa0_df = CSV.read(kappa0_filename, DataFrame)
    #     kappa0s = kappa0_df[:, "kappa0"]
    #     tcs = kappa0_df[:, "tᶜ"]
    # else
    #     kappa0s = get(npi_params_dict, "κ₀s", [0.0])
    #     tcs = get(npi_params_dict, "tᶜs", [typemax(Int)])
    # end

    return NPI_Params(qm)
end
