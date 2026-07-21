module EpiCommute

using DataFrames, CSV, Distributions, Random, Dates, Printf
using NetCDF, HDF5

include("utils.jl")
include("npi.jl")
include("model.jl")
include("io.jl")

# Exported types and functions (minimum API expected by EpiSim)
export Epidemic_Params, Population_Params, NPI_Params,
       init_pop_param_struct, init_epi_parameters_struct, init_NPI_parameters_struct,
       create_config_template,
       set_compartments!, create_initial_compartments_dict,
       run_epidemic_spreading!,
       save_simulation_netCDF, save_simulation_hdf5, save_observables_netCDF,
       correct_self_loops

end # module
