include("Include.jl")
compound="ethanol"
# load the data dictionary -
# data_dictionary = maximize_acetate_data_dictionary(0,0,0)
# data_dictionary = maximize_atp_data_dictionary(0,0,0)
# data_dictionary = maximize_cellmass_data_dictionary(0,0,0)
# data_dictionary = maximize_formate_data_dictionary(0,0,0)
data_dictionary = maximize_ethanol_data_dictionary(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

#calculates the residual
path_to_atom_file="./Atom.txt"
atom_matrix=generate_atom_matrix(path_to_atom_file,data_dictionary);
# path_to_atom_array = "./atom_array.dat"
# writedlm(path_to_atom_array, atom_matrix)
residual=transpose(atom_matrix)*uptake_array[:,1]
file_path_residual = "./figs/residual_"*compound*".dat"
writedlm(file_path_residual,residual)
@show residual
generate_net_reaction_string(uptake_array,0.01,data_dictionary)


# find_missing_metabolites_in_atom_dataset(path_to_atom_file,data_dictionary)
