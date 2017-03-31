# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[24] = -1.0

	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]
	default_flux_bounds_array[21,2] = 0.0

	# ATP maintenance -
	# default_flux_bounds_array[20,1:2] = 7.6

	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	20.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	100.0	;	# 76 M_co2_b
		0.0	0.0	;	# 77 M_etoh_b
		0.0	0.0	;	# 78 M_for_b
		0.0	0.0	;	# 79 M_fru_b
		0.0	1.0	;	# 80 M_fum_b
		-10.0	0.0	;	# 81 M_glc_D_b
		0.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-100.0	100.0	;	# 84 M_h2o_b
		-100.0	100.0	;	# 85 M_h_b
		0.0	0.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-10.0	10.0	;	# 88 M_nh4_b
		-15.0	0.0	;	# 89 M_o2_b
		-10.0	10.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]

	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end

	return data_dictionary
end


function maximize_atp_data_dictionary(time_start,time_stop,time_step)

	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[20] = -1.0

	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]
	default_flux_bounds_array[21,2] = 0.0

	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	0.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	100.0	;	# 76 M_co2_b
		0.0	0.0	;	# 77 M_etoh_b
		0.0	1.0	;	# 78 M_for_b
		0.0	1.0	;	# 79 M_fru_b
		0.0	1.0	;	# 80 M_fum_b
		-1.0	0.0	;	# 81 M_glc_D_b
		0.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-1.0	1.0	;	# 84 M_h2o_b
		-1.0	1.0	;	# 85 M_h_b
		0.0	0.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-1.0	1.0	;	# 88 M_nh4_b
		-100.0	0.0	;	# 89 M_o2_b
		-1.0	1.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]

	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end

	return data_dictionary
end

function maximize_acetate_data_dictionary(time_start,time_stop,time_step)
compound ="acetate"
	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[35] = -1.0

	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	10.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	1.0	;	# 76 M_co2_b
		0.0	0.0	;	# 77 M_etoh_b
		0.0	1.0	;	# 78 M_for_b
		0.0	1.0	;	# 79 M_fru_b
		0.0	1.0	;	# 80 M_fum_b
		-1.0	0.0	;	# 81 M_glc_D_b
		-1.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-1.0	1.0	;	# 84 M_h2o_b
		-1.0	1.0	;	# 85 M_h_b
		0.0	0.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-1.0	1.0	;	# 88 M_nh4_b
		-1.0	1.0	;	# 89 M_o2_b
		-1.0	1.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]

	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end

	return data_dictionary
end

function maximize_formate_data_dictionary(time_start,time_stop,time_step)

	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[41] = -1.0

	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	10.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	10.0	;	# 76 M_co2_b
		0.0	0.0	;	# 77 M_etoh_b
		0.0	10.0	;	# 78 M_for_b
		0.0	1.0	;	# 79 M_fru_b
		0.0	1.0	;	# 80 M_fum_b
		-1.0	0.0	;	# 81 M_glc_D_b
		-1.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-10.0	10.0	;	# 84 M_h2o_b
		-10.0	10.0	;	# 85 M_h_b
		0.0	10.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-10.0	10.0	;	# 88 M_nh4_b
		-10.0	10.0	;	# 89 M_o2_b
		-10.0	10.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]

	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end

	return data_dictionary
end
function maximize_ethanol_data_dictionary(time_start,time_stop,time_step)

	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[40] = -1.0

	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	1.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	10.0	;	# 76 M_co2_b
		0.0	10.0	;	# 77 M_etoh_b
		0.0	1.0	;	# 78 M_for_b
		0.0	1.0	;	# 79 M_fru_b
		0.0	1.0	;	# 80 M_fum_b
		-1.0	0.0	;	# 81 M_glc_D_b
		-1.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-10.0	10.0	;	# 84 M_h2o_b
		-10.0	10.0	;	# 85 M_h_b
		0.0	0.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-10.0	10.0	;	# 88 M_nh4_b
		-10.0	10.0	;	# 89 M_o2_b
		-10.0	10.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]

	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end

	return data_dictionary
end

# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-03-30T21:03:35.549
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# Setup default flux bounds array -
	default_bounds_array = [
		0	100.0	;	# 1 M_acald_c+M_coa_c+M_nad_c --> M_accoa_c+M_h_c+M_nadh_c
		0	100.0	;	# 2 M_accoa_c+M_h_c+M_nadh_c --> M_acald_c+M_coa_c+M_nad_c
		0	100.0	;	# 3 M_acald_e --> M_acald_c
		0	100.0	;	# 4 M_acald_c --> M_acald_e
		0	100.0	;	# 5 M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0	100.0	;	# 6 M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0	100.0	;	# 7 M_cit_c --> M_acon_C_c+M_h2o_c
		0	100.0	;	# 8 M_acon_C_c+M_h2o_c --> M_cit_c
		0	100.0	;	# 9 M_acon_C_c+M_h2o_c --> M_icit_c
		0	100.0	;	# 10 M_icit_c --> M_acon_C_c+M_h2o_c
		0	100.0	;	# 11 M_ac_e+M_h_e --> M_ac_c+M_h_c
		0	100.0	;	# 12 M_ac_c+M_h_c --> M_ac_e+M_h_e
		0	100.0	;	# 13 M_amp_c+M_atp_c --> 2.0*M_adp_c
		0	100.0	;	# 14 2.0*M_adp_c --> M_amp_c+M_atp_c
		0	100.0	;	# 15 M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0	100.0	;	# 16 M_akg_e+M_h_e --> M_akg_c+M_h_c
		0	100.0	;	# 17 M_akg_c+M_h_c --> M_akg_e+M_h_e
		0	100.0	;	# 18 M_etoh_c+M_nad_c --> M_acald_c+M_h_c+M_nadh_c
		0	100.0	;	# 19 M_acald_c+M_h_c+M_nadh_c --> M_etoh_c+M_nad_c
		0	100.0	;	# 20 M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0	100.0	;	# 21 M_adp_c+M_pi_c --> M_atp_c+M_h2o_c
		0	100.0	;	# 22 M_adp_c+4.0*M_h_e+M_pi_c --> M_atp_c+M_h2o_c+4.0*M_h_c
		0	100.0	;	# 23 M_atp_c+M_h2o_c+4.0*M_h_c --> M_adp_c+4.0*M_h_e+M_pi_c
		0	100.0	;	# 24 1.496*M_3pg_c+3.7478*M_accoa_c+59.81*M_atp_c+0.361*M_e4p_c+0.0709*M_f6p_c+0.129*M_g3p_c+0.205*M_g6p_c+0.2557*M_gln_L_c+4.9414*M_glu_L_c+59.81*M_h2o_c+3.547*M_nad_c+13.0279*M_nadph_c+1.7867*M_oaa_c+0.5191*M_pep_c+2.8328*M_pyr_c+0.8977*M_r5p_c --> 59.81*M_adp_c+4.1182*M_akg_c+3.7478*M_coa_c+59.81*M_h_c+3.547*M_nadh_c+13.0279*M_nadp_c+59.81*M_pi_c
		0	100.0	;	# 25 M_co2_e --> M_co2_c
		0	100.0	;	# 26 M_co2_c --> M_co2_e
		0	100.0	;	# 27 M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0	100.0	;	# 28 2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+2.0*M_h_e+M_q8_c
		0	100.0	;	# 29 M_h_e+M_lac_D_e --> M_h_c+M_lac_D_c
		0	100.0	;	# 30 M_h_c+M_lac_D_c --> M_h_e+M_lac_D_e
		0	100.0	;	# 31 M_2pg_c --> M_h2o_c+M_pep_c
		0	100.0	;	# 32 M_h2o_c+M_pep_c --> M_2pg_c
		0	100.0	;	# 33 M_etoh_e+M_h_e --> M_etoh_c+M_h_c
		0	100.0	;	# 34 M_etoh_c+M_h_c --> M_etoh_e+M_h_e
		0	100.0	;	# 35 M_ac_e --> M_ac_b
		0	100.0	;	# 36 M_acald_e --> M_acald_b
		0	100.0	;	# 37 M_akg_e --> M_akg_b
		0	100.0	;	# 38 M_co2_e --> M_co2_b
		0	100.0	;	# 39 M_co2_b --> M_co2_e
		0	100.0	;	# 40 M_etoh_e --> M_etoh_b
		0	100.0	;	# 41 M_for_e --> M_for_b
		0	100.0	;	# 42 M_fru_e --> M_fru_b
		0	100.0	;	# 43 M_fum_e --> M_fum_b
		0	100.0	;	# 44 M_glc_D_e --> M_glc_D_b
		0	100.0	;	# 45 M_glc_D_b --> M_glc_D_e
		0	100.0	;	# 46 M_gln_L_e --> M_gln_L_b
		0	100.0	;	# 47 M_glu_L_e --> M_glu_L_b
		0	100.0	;	# 48 M_h_e --> M_h_b
		0	100.0	;	# 49 M_h_b --> M_h_e
		0	100.0	;	# 50 M_h2o_e --> M_h2o_b
		0	100.0	;	# 51 M_h2o_b --> M_h2o_e
		0	100.0	;	# 52 M_lac_D_e --> M_lac_D_b
		0	100.0	;	# 53 M_mal_L_e --> M_mal_L_b
		0	100.0	;	# 54 M_nh4_e --> M_nh4_b
		0	100.0	;	# 55 M_nh4_b --> M_nh4_e
		0	100.0	;	# 56 M_o2_e --> M_o2_b
		0	100.0	;	# 57 M_o2_b --> M_o2_e
		0	100.0	;	# 58 M_pi_e --> M_pi_b
		0	100.0	;	# 59 M_pi_b --> M_pi_e
		0	100.0	;	# 60 M_pyr_e --> M_pyr_b
		0	100.0	;	# 61 M_succ_e --> M_succ_b
		0	100.0	;	# 62 M_fdp_c --> M_dhap_c+M_g3p_c
		0	100.0	;	# 63 M_dhap_c+M_g3p_c --> M_fdp_c
		0	100.0	;	# 64 M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0	100.0	;	# 65 M_for_e+M_h_e --> M_for_c+M_h_c
		0	100.0	;	# 66 M_for_c --> M_for_e
		0	100.0	;	# 67 M_fum_c+M_q8h2_c --> M_q8_c+M_succ_c
		0	100.0	;	# 68 M_fru_e+M_pep_c --> M_f6p_c+M_pyr_c
		0	100.0	;	# 69 M_fum_c+M_h2o_c --> M_mal_L_c
		0	100.0	;	# 70 M_mal_L_c --> M_fum_c+M_h2o_c
		0	100.0	;	# 71 M_fum_e+2.0*M_h_e --> M_fum_c+2.0*M_h_c
		0	100.0	;	# 72 M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0	100.0	;	# 73 M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0	100.0	;	# 74 M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0	100.0	;	# 75 M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0	100.0	;	# 76 M_glc_D_e+M_pep_c --> M_g6p_c+M_pyr_c
		0	100.0	;	# 77 M_atp_c+M_glu_L_c+M_nh4_c --> M_adp_c+M_gln_L_c+M_h_c+M_pi_c
		0	100.0	;	# 78 M_atp_c+M_gln_L_e+M_h2o_c --> M_adp_c+M_gln_L_c+M_pi_c
		0	100.0	;	# 79 M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh4_c
		0	100.0	;	# 80 M_akg_c+M_nadph_c+M_nh4_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0	100.0	;	# 81 M_gln_L_c+M_h2o_c+M_h_c --> M_glu_L_c+M_nh4_c
		0	100.0	;	# 82 M_akg_c+M_gln_L_c+M_h_c+M_nadph_c --> 2.0*M_glu_L_c+M_nadp_c
		0	100.0	;	# 83 M_glu_L_e+M_h_e --> M_glu_L_c+M_h_c
		0	100.0	;	# 84 M_glu_L_c+M_h_c --> M_glu_L_e+M_h_e
		0	100.0	;	# 85 M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0	100.0	;	# 86 M_h2o_e --> M_h2o_c
		0	100.0	;	# 87 M_h2o_c --> M_h2o_e
		0	100.0	;	# 88 M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0	100.0	;	# 89 M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0	100.0	;	# 90 M_icit_c --> M_glx_c+M_succ_c
		0	100.0	;	# 91 M_lac_D_c+M_nad_c --> M_h_c+M_nadh_c+M_pyr_c
		0	100.0	;	# 92 M_h_c+M_nadh_c+M_pyr_c --> M_lac_D_c+M_nad_c
		0	100.0	;	# 93 M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0	100.0	;	# 94 2.0*M_h_e+M_mal_L_e --> 2.0*M_h_c+M_mal_L_c
		0	100.0	;	# 95 M_mal_L_c+M_nad_c --> M_h_c+M_nadh_c+M_oaa_c
		0	100.0	;	# 96 M_h_c+M_nadh_c+M_oaa_c --> M_mal_L_c+M_nad_c
		0	100.0	;	# 97 M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0	100.0	;	# 98 M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0	100.0	;	# 99 4.0*M_h_c+M_nadh_c+M_q8_c --> 3.0*M_h_e+M_nad_c+M_q8h2_c
		0	100.0	;	# 100 M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0	100.0	;	# 101 M_nh4_e --> M_nh4_c
		0	100.0	;	# 102 M_nh4_c --> M_nh4_e
		0	100.0	;	# 103 M_o2_e --> M_o2_c
		0	100.0	;	# 104 M_o2_c --> M_o2_e
		0	100.0	;	# 105 M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0	100.0	;	# 106 M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0	100.0	;	# 107 M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0	100.0	;	# 108 M_g6p_c --> M_f6p_c
		0	100.0	;	# 109 M_f6p_c --> M_g6p_c
		0	100.0	;	# 110 M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0	100.0	;	# 111 M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0	100.0	;	# 112 M_6pgl_c+M_h2o_c --> M_6pgc_c
		0	100.0	;	# 113 M_2pg_c --> M_3pg_c
		0	100.0	;	# 114 M_3pg_c --> M_2pg_c
		0	100.0	;	# 115 M_h_e+M_pi_e --> M_h_c+M_pi_c
		0	100.0	;	# 116 M_h_c+M_pi_c --> M_h_e+M_pi_e
		0	100.0	;	# 117 M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0	100.0	;	# 118 M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0	100.0	;	# 119 M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0	100.0	;	# 120 M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0	100.0	;	# 121 M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0	100.0	;	# 122 M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0	100.0	;	# 123 M_h_e+M_pyr_e --> M_h_c+M_pyr_c
		0	100.0	;	# 124 M_h_c+M_pyr_c --> M_h_e+M_pyr_e
		0	100.0	;	# 125 M_ru5p_D_c --> M_xu5p_D_c
		0	100.0	;	# 126 M_xu5p_D_c --> M_ru5p_D_c
		0	100.0	;	# 127 M_r5p_c --> M_ru5p_D_c
		0	100.0	;	# 128 M_ru5p_D_c --> M_r5p_c
		0	100.0	;	# 129 2.0*M_h_e+M_succ_e --> 2.0*M_h_c+M_succ_c
		0	100.0	;	# 130 M_h_e+M_succ_c --> M_h_c+M_succ_e
		0	100.0	;	# 131 M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0	100.0	;	# 132 M_atp_c+M_coa_c+M_succ_c --> M_adp_c+M_pi_c+M_succoa_c
		0	100.0	;	# 133 M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0	100.0	;	# 134 M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0	100.0	;	# 135 M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 136 2.0*M_h_e+M_nadh_c+M_nadp_c --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0	100.0	;	# 137 M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 138 M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0	100.0	;	# 139 M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0	100.0	;	# 140 M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0	100.0	;	# 141 M_dhap_c --> M_g3p_c
		0	100.0	;	# 142 M_g3p_c --> M_dhap_c
	];

	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1 M_13dpg_c
		0.0	0.0	;	# 2 M_2pg_c
		0.0	0.0	;	# 3 M_3pg_c
		0.0	0.0	;	# 4 M_6pgc_c
		0.0	0.0	;	# 5 M_6pgl_c
		0.0	0.0	;	# 6 M_ac_c
		0.0	0.0	;	# 7 M_ac_e
		0.0	0.0	;	# 8 M_acald_c
		0.0	0.0	;	# 9 M_acald_e
		0.0	0.0	;	# 10 M_accoa_c
		0.0	0.0	;	# 11 M_acon_C_c
		0.0	0.0	;	# 12 M_actp_c
		0.0	0.0	;	# 13 M_adp_c
		0.0	0.0	;	# 14 M_akg_c
		0.0	0.0	;	# 15 M_akg_e
		0.0	0.0	;	# 16 M_amp_c
		0.0	0.0	;	# 17 M_atp_c
		0.0	0.0	;	# 18 M_cit_c
		0.0	0.0	;	# 19 M_co2_c
		0.0	0.0	;	# 20 M_co2_e
		0.0	0.0	;	# 21 M_coa_c
		0.0	0.0	;	# 22 M_dhap_c
		0.0	0.0	;	# 23 M_e4p_c
		0.0	0.0	;	# 24 M_etoh_c
		0.0	0.0	;	# 25 M_etoh_e
		0.0	0.0	;	# 26 M_f6p_c
		0.0	0.0	;	# 27 M_fdp_c
		0.0	0.0	;	# 28 M_for_c
		0.0	0.0	;	# 29 M_for_e
		0.0	0.0	;	# 30 M_fru_e
		0.0	0.0	;	# 31 M_fum_c
		0.0	0.0	;	# 32 M_fum_e
		0.0	0.0	;	# 33 M_g3p_c
		0.0	0.0	;	# 34 M_g6p_c
		0.0	0.0	;	# 35 M_glc_D_e
		0.0	0.0	;	# 36 M_gln_L_c
		0.0	0.0	;	# 37 M_gln_L_e
		0.0	0.0	;	# 38 M_glu_L_c
		0.0	0.0	;	# 39 M_glu_L_e
		0.0	0.0	;	# 40 M_glx_c
		0.0	0.0	;	# 41 M_h2o_c
		0.0	0.0	;	# 42 M_h2o_e
		0.0	0.0	;	# 43 M_h_c
		0.0	0.0	;	# 44 M_h_e
		0.0	0.0	;	# 45 M_icit_c
		0.0	0.0	;	# 46 M_lac_D_c
		0.0	0.0	;	# 47 M_lac_D_e
		0.0	0.0	;	# 48 M_mal_L_c
		0.0	0.0	;	# 49 M_mal_L_e
		0.0	0.0	;	# 50 M_nad_c
		0.0	0.0	;	# 51 M_nadh_c
		0.0	0.0	;	# 52 M_nadp_c
		0.0	0.0	;	# 53 M_nadph_c
		0.0	0.0	;	# 54 M_nh4_c
		0.0	0.0	;	# 55 M_nh4_e
		0.0	0.0	;	# 56 M_o2_c
		0.0	0.0	;	# 57 M_o2_e
		0.0	0.0	;	# 58 M_oaa_c
		0.0	0.0	;	# 59 M_pep_c
		0.0	0.0	;	# 60 M_pi_c
		0.0	0.0	;	# 61 M_pi_e
		0.0	0.0	;	# 62 M_pyr_c
		0.0	0.0	;	# 63 M_pyr_e
		0.0	0.0	;	# 64 M_q8_c
		0.0	0.0	;	# 65 M_q8h2_c
		0.0	0.0	;	# 66 M_r5p_c
		0.0	0.0	;	# 67 M_ru5p_D_c
		0.0	0.0	;	# 68 M_s7p_c
		0.0	0.0	;	# 69 M_succ_c
		0.0	0.0	;	# 70 M_succ_e
		0.0	0.0	;	# 71 M_succoa_c
		0.0	0.0	;	# 72 M_xu5p_D_c
		-1.0	1.0	;	# 73 M_ac_b
		-1.0	1.0	;	# 74 M_acald_b
		-1.0	1.0	;	# 75 M_akg_b
		-1.0	1.0	;	# 76 M_co2_b
		-1.0	1.0	;	# 77 M_etoh_b
		-1.0	1.0	;	# 78 M_for_b
		-1.0	1.0	;	# 79 M_fru_b
		-1.0	1.0	;	# 80 M_fum_b
		-1.0	1.0	;	# 81 M_glc_D_b
		-1.0	1.0	;	# 82 M_gln_L_b
		-1.0	1.0	;	# 83 M_glu_L_b
		-1.0	1.0	;	# 84 M_h2o_b
		-1.0	1.0	;	# 85 M_h_b
		-1.0	1.0	;	# 86 M_lac_D_b
		-1.0	1.0	;	# 87 M_mal_L_b
		-1.0	1.0	;	# 88 M_nh4_b
		-1.0	1.0	;	# 89 M_o2_b
		-1.0	1.0	;	# 90 M_pi_b
		-1.0	1.0	;	# 91 M_pyr_b
		-1.0	1.0	;	# 92 M_succ_b
	];

	# Setup the objective coefficient array -
	objective_coefficient_array = [
		0.0	;	# 1 R_ACALD::M_acald_c+M_coa_c+M_nad_c --> M_accoa_c+M_h_c+M_nadh_c
		0.0	;	# 2 R_ACALD_reverse::M_accoa_c+M_h_c+M_nadh_c --> M_acald_c+M_coa_c+M_nad_c
		0.0	;	# 3 R_ACALDt::M_acald_e --> M_acald_c
		0.0	;	# 4 R_ACALDt_reverse::M_acald_c --> M_acald_e
		0.0	;	# 5 R_ACKr::M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0.0	;	# 6 R_ACKr_reverse::M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0.0	;	# 7 R_ACONTa::M_cit_c --> M_acon_C_c+M_h2o_c
		0.0	;	# 8 R_ACONTa_reverse::M_acon_C_c+M_h2o_c --> M_cit_c
		0.0	;	# 9 R_ACONTb::M_acon_C_c+M_h2o_c --> M_icit_c
		0.0	;	# 10 R_ACONTb_reverse::M_icit_c --> M_acon_C_c+M_h2o_c
		0.0	;	# 11 R_ACt2r::M_ac_e+M_h_e --> M_ac_c+M_h_c
		0.0	;	# 12 R_ACt2r_reverse::M_ac_c+M_h_c --> M_ac_e+M_h_e
		0.0	;	# 13 R_ADK1::M_amp_c+M_atp_c --> 2.0*M_adp_c
		0.0	;	# 14 R_ADK1_reverse::2.0*M_adp_c --> M_amp_c+M_atp_c
		0.0	;	# 15 R_AKGDH::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0.0	;	# 16 R_AKGt2r::M_akg_e+M_h_e --> M_akg_c+M_h_c
		0.0	;	# 17 R_AKGt2r_reverse::M_akg_c+M_h_c --> M_akg_e+M_h_e
		0.0	;	# 18 R_ALCD2x::M_etoh_c+M_nad_c --> M_acald_c+M_h_c+M_nadh_c
		0.0	;	# 19 R_ALCD2x_reverse::M_acald_c+M_h_c+M_nadh_c --> M_etoh_c+M_nad_c
		0.0	;	# 20 R_ATPM::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0.0	;	# 21 R_ATPM_reverse::M_adp_c+M_pi_c --> M_atp_c+M_h2o_c
		0.0	;	# 22 R_ATPS4r::M_adp_c+4.0*M_h_e+M_pi_c --> M_atp_c+M_h2o_c+4.0*M_h_c
		0.0	;	# 23 R_ATPS4r_reverse::M_atp_c+M_h2o_c+4.0*M_h_c --> M_adp_c+4.0*M_h_e+M_pi_c
		0.0	;	# 24 R_Biomass_Ecoli_core_w_GAM::1.496*M_3pg_c+3.7478*M_accoa_c+59.81*M_atp_c+0.361*M_e4p_c+0.0709*M_f6p_c+0.129*M_g3p_c+0.205*M_g6p_c+0.2557*M_gln_L_c+4.9414*M_glu_L_c+59.81*M_h2o_c+3.547*M_nad_c+13.0279*M_nadph_c+1.7867*M_oaa_c+0.5191*M_pep_c+2.8328*M_pyr_c+0.8977*M_r5p_c --> 59.81*M_adp_c+4.1182*M_akg_c+3.7478*M_coa_c+59.81*M_h_c+3.547*M_nadh_c+13.0279*M_nadp_c+59.81*M_pi_c
		0.0	;	# 25 R_CO2t::M_co2_e --> M_co2_c
		0.0	;	# 26 R_CO2t_reverse::M_co2_c --> M_co2_e
		0.0	;	# 27 R_CS::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0.0	;	# 28 R_CYTBD::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+2.0*M_h_e+M_q8_c
		0.0	;	# 29 R_D_LACt2::M_h_e+M_lac_D_e --> M_h_c+M_lac_D_c
		0.0	;	# 30 R_D_LACt2_reverse::M_h_c+M_lac_D_c --> M_h_e+M_lac_D_e
		0.0	;	# 31 R_ENO::M_2pg_c --> M_h2o_c+M_pep_c
		0.0	;	# 32 R_ENO_reverse::M_h2o_c+M_pep_c --> M_2pg_c
		0.0	;	# 33 R_ETOHt2r::M_etoh_e+M_h_e --> M_etoh_c+M_h_c
		0.0	;	# 34 R_ETOHt2r_reverse::M_etoh_c+M_h_c --> M_etoh_e+M_h_e
		0.0	;	# 35 R_EX_ac_e::M_ac_e --> M_ac_b
		0.0	;	# 36 R_EX_acald_e::M_acald_e --> M_acald_b
		0.0	;	# 37 R_EX_akg_e::M_akg_e --> M_akg_b
		0.0	;	# 38 R_EX_co2_e::M_co2_e --> M_co2_b
		0.0	;	# 39 R_EX_co2_e_reverse::M_co2_b --> M_co2_e
		0.0	;	# 40 R_EX_etoh_e::M_etoh_e --> M_etoh_b
		0.0	;	# 41 R_EX_for_e::M_for_e --> M_for_b
		0.0	;	# 42 R_EX_fru_e::M_fru_e --> M_fru_b
		0.0	;	# 43 R_EX_fum_e::M_fum_e --> M_fum_b
		0.0	;	# 44 R_EX_glc_e::M_glc_D_e --> M_glc_D_b
		0.0	;	# 45 R_EX_glc_e_reverse::M_glc_D_b --> M_glc_D_e
		0.0	;	# 46 R_EX_gln_L_e::M_gln_L_e --> M_gln_L_b
		0.0	;	# 47 R_EX_glu_L_e::M_glu_L_e --> M_glu_L_b
		0.0	;	# 48 R_EX_h_e::M_h_e --> M_h_b
		0.0	;	# 49 R_EX_h_e_reverse::M_h_b --> M_h_e
		0.0	;	# 50 R_EX_h2o_e::M_h2o_e --> M_h2o_b
		0.0	;	# 51 R_EX_h2o_e_reverse::M_h2o_b --> M_h2o_e
		0.0	;	# 52 R_EX_lac_D_e::M_lac_D_e --> M_lac_D_b
		0.0	;	# 53 R_EX_mal_L_e::M_mal_L_e --> M_mal_L_b
		0.0	;	# 54 R_EX_nh4_e::M_nh4_e --> M_nh4_b
		0.0	;	# 55 R_EX_nh4_e_reverse::M_nh4_b --> M_nh4_e
		0.0	;	# 56 R_EX_o2_e::M_o2_e --> M_o2_b
		0.0	;	# 57 R_EX_o2_e_reverse::M_o2_b --> M_o2_e
		0.0	;	# 58 R_EX_pi_e::M_pi_e --> M_pi_b
		0.0	;	# 59 R_EX_pi_e_reverse::M_pi_b --> M_pi_e
		0.0	;	# 60 R_EX_pyr_e::M_pyr_e --> M_pyr_b
		0.0	;	# 61 R_EX_succ_e::M_succ_e --> M_succ_b
		0.0	;	# 62 R_FBA::M_fdp_c --> M_dhap_c+M_g3p_c
		0.0	;	# 63 R_FBA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c
		0.0	;	# 64 R_FBP::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0.0	;	# 65 R_FORt2::M_for_e+M_h_e --> M_for_c+M_h_c
		0.0	;	# 66 R_FORti::M_for_c --> M_for_e
		0.0	;	# 67 R_FRD7::M_fum_c+M_q8h2_c --> M_q8_c+M_succ_c
		0.0	;	# 68 R_FRUpts2::M_fru_e+M_pep_c --> M_f6p_c+M_pyr_c
		0.0	;	# 69 R_FUM::M_fum_c+M_h2o_c --> M_mal_L_c
		0.0	;	# 70 R_FUM_reverse::M_mal_L_c --> M_fum_c+M_h2o_c
		0.0	;	# 71 R_FUMt2_2::M_fum_e+2.0*M_h_e --> M_fum_c+2.0*M_h_c
		0.0	;	# 72 R_G6PDH2r::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0.0	;	# 73 R_G6PDH2r_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0.0	;	# 74 R_GAPD::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0.0	;	# 75 R_GAPD_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0.0	;	# 76 R_GLCpts::M_glc_D_e+M_pep_c --> M_g6p_c+M_pyr_c
		0.0	;	# 77 R_GLNS::M_atp_c+M_glu_L_c+M_nh4_c --> M_adp_c+M_gln_L_c+M_h_c+M_pi_c
		0.0	;	# 78 R_GLNabc::M_atp_c+M_gln_L_e+M_h2o_c --> M_adp_c+M_gln_L_c+M_pi_c
		0.0	;	# 79 R_GLUDy::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh4_c
		0.0	;	# 80 R_GLUDy_reverse::M_akg_c+M_nadph_c+M_nh4_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0.0	;	# 81 R_GLUN::M_gln_L_c+M_h2o_c+M_h_c --> M_glu_L_c+M_nh4_c
		0.0	;	# 82 R_GLUSy::M_akg_c+M_gln_L_c+M_h_c+M_nadph_c --> 2.0*M_glu_L_c+M_nadp_c
		0.0	;	# 83 R_GLUt2r::M_glu_L_e+M_h_e --> M_glu_L_c+M_h_c
		0.0	;	# 84 R_GLUt2r_reverse::M_glu_L_c+M_h_c --> M_glu_L_e+M_h_e
		0.0	;	# 85 R_GND::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0.0	;	# 86 R_H2Ot::M_h2o_e --> M_h2o_c
		0.0	;	# 87 R_H2Ot_reverse::M_h2o_c --> M_h2o_e
		0.0	;	# 88 R_ICDHyr::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0.0	;	# 89 R_ICDHyr_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0.0	;	# 90 R_ICL::M_icit_c --> M_glx_c+M_succ_c
		0.0	;	# 91 R_LDH_D::M_lac_D_c+M_nad_c --> M_h_c+M_nadh_c+M_pyr_c
		0.0	;	# 92 R_LDH_D_reverse::M_h_c+M_nadh_c+M_pyr_c --> M_lac_D_c+M_nad_c
		0.0	;	# 93 R_MALS::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0.0	;	# 94 R_MALt2_2::2.0*M_h_e+M_mal_L_e --> 2.0*M_h_c+M_mal_L_c
		0.0	;	# 95 R_MDH::M_mal_L_c+M_nad_c --> M_h_c+M_nadh_c+M_oaa_c
		0.0	;	# 96 R_MDH_reverse::M_h_c+M_nadh_c+M_oaa_c --> M_mal_L_c+M_nad_c
		0.0	;	# 97 R_ME1::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0.0	;	# 98 R_ME2::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0.0	;	# 99 R_NADH16::4.0*M_h_c+M_nadh_c+M_q8_c --> 3.0*M_h_e+M_nad_c+M_q8h2_c
		0.0	;	# 100 R_NADTRHD::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0.0	;	# 101 R_NH4t::M_nh4_e --> M_nh4_c
		0.0	;	# 102 R_NH4t_reverse::M_nh4_c --> M_nh4_e
		0.0	;	# 103 R_O2t::M_o2_e --> M_o2_c
		0.0	;	# 104 R_O2t_reverse::M_o2_c --> M_o2_e
		0.0	;	# 105 R_PDH::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0.0	;	# 106 R_PFK::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0.0	;	# 107 R_PFL::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0.0	;	# 108 R_PGI::M_g6p_c --> M_f6p_c
		0.0	;	# 109 R_PGI_reverse::M_f6p_c --> M_g6p_c
		0.0	;	# 110 R_PGK::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0.0	;	# 111 R_PGK_reverse::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0.0	;	# 112 R_PGL::M_6pgl_c+M_h2o_c --> M_6pgc_c
		0.0	;	# 113 R_PGM::M_2pg_c --> M_3pg_c
		0.0	;	# 114 R_PGM_reverse::M_3pg_c --> M_2pg_c
		0.0	;	# 115 R_PIt2r::M_h_e+M_pi_e --> M_h_c+M_pi_c
		0.0	;	# 116 R_PIt2r_reverse::M_h_c+M_pi_c --> M_h_e+M_pi_e
		0.0	;	# 117 R_PPC::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0.0	;	# 118 R_PPCK::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0.0	;	# 119 R_PPS::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0.0	;	# 120 R_PTAr::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0.0	;	# 121 R_PTAr_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0.0	;	# 122 R_PYK::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0.0	;	# 123 R_PYRt2r::M_h_e+M_pyr_e --> M_h_c+M_pyr_c
		0.0	;	# 124 R_PYRt2r_reverse::M_h_c+M_pyr_c --> M_h_e+M_pyr_e
		0.0	;	# 125 R_RPE::M_ru5p_D_c --> M_xu5p_D_c
		0.0	;	# 126 R_RPE_reverse::M_xu5p_D_c --> M_ru5p_D_c
		0.0	;	# 127 R_RPI::M_r5p_c --> M_ru5p_D_c
		0.0	;	# 128 R_RPI_reverse::M_ru5p_D_c --> M_r5p_c
		0.0	;	# 129 R_SUCCt2_2::2.0*M_h_e+M_succ_e --> 2.0*M_h_c+M_succ_c
		0.0	;	# 130 R_SUCCt3::M_h_e+M_succ_c --> M_h_c+M_succ_e
		0.0	;	# 131 R_SUCDi::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0.0	;	# 132 R_SUCOAS::M_atp_c+M_coa_c+M_succ_c --> M_adp_c+M_pi_c+M_succoa_c
		0.0	;	# 133 R_SUCOAS_reverse::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0.0	;	# 134 R_TALA::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0.0	;	# 135 R_TALA_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0.0	;	# 136 R_THD2::2.0*M_h_e+M_nadh_c+M_nadp_c --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0.0	;	# 137 R_TKT1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0.0	;	# 138 R_TKT1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0.0	;	# 139 R_TKT2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0.0	;	# 140 R_TKT2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0.0	;	# 141 R_TPI::M_dhap_c --> M_g3p_c
		0.0	;	# 142 R_TPI_reverse::M_g3p_c --> M_dhap_c
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"R_ACALD::M_acald_c+M_coa_c+M_nad_c --> M_accoa_c+M_h_c+M_nadh_c"
		"R_ACALD_reverse::M_accoa_c+M_h_c+M_nadh_c --> M_acald_c+M_coa_c+M_nad_c"
		"R_ACALDt::M_acald_e --> M_acald_c"
		"R_ACALDt_reverse::M_acald_c --> M_acald_e"
		"R_ACKr::M_ac_c+M_atp_c --> M_actp_c+M_adp_c"
		"R_ACKr_reverse::M_actp_c+M_adp_c --> M_ac_c+M_atp_c"
		"R_ACONTa::M_cit_c --> M_acon_C_c+M_h2o_c"
		"R_ACONTa_reverse::M_acon_C_c+M_h2o_c --> M_cit_c"
		"R_ACONTb::M_acon_C_c+M_h2o_c --> M_icit_c"
		"R_ACONTb_reverse::M_icit_c --> M_acon_C_c+M_h2o_c"
		"R_ACt2r::M_ac_e+M_h_e --> M_ac_c+M_h_c"
		"R_ACt2r_reverse::M_ac_c+M_h_c --> M_ac_e+M_h_e"
		"R_ADK1::M_amp_c+M_atp_c --> 2.0*M_adp_c"
		"R_ADK1_reverse::2.0*M_adp_c --> M_amp_c+M_atp_c"
		"R_AKGDH::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c"
		"R_AKGt2r::M_akg_e+M_h_e --> M_akg_c+M_h_c"
		"R_AKGt2r_reverse::M_akg_c+M_h_c --> M_akg_e+M_h_e"
		"R_ALCD2x::M_etoh_c+M_nad_c --> M_acald_c+M_h_c+M_nadh_c"
		"R_ALCD2x_reverse::M_acald_c+M_h_c+M_nadh_c --> M_etoh_c+M_nad_c"
		"R_ATPM::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c"
		"R_ATPM_reverse::M_adp_c+M_pi_c --> M_atp_c+M_h2o_c"
		"R_ATPS4r::M_adp_c+4.0*M_h_e+M_pi_c --> M_atp_c+M_h2o_c+4.0*M_h_c"
		"R_ATPS4r_reverse::M_atp_c+M_h2o_c+4.0*M_h_c --> M_adp_c+4.0*M_h_e+M_pi_c"
		"R_Biomass_Ecoli_core_w_GAM::1.496*M_3pg_c+3.7478*M_accoa_c+59.81*M_atp_c+0.361*M_e4p_c+0.0709*M_f6p_c+0.129*M_g3p_c+0.205*M_g6p_c+0.2557*M_gln_L_c+4.9414*M_glu_L_c+59.81*M_h2o_c+3.547*M_nad_c+13.0279*M_nadph_c+1.7867*M_oaa_c+0.5191*M_pep_c+2.8328*M_pyr_c+0.8977*M_r5p_c --> 59.81*M_adp_c+4.1182*M_akg_c+3.7478*M_coa_c+59.81*M_h_c+3.547*M_nadh_c+13.0279*M_nadp_c+59.81*M_pi_c"
		"R_CO2t::M_co2_e --> M_co2_c"
		"R_CO2t_reverse::M_co2_c --> M_co2_e"
		"R_CS::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c"
		"R_CYTBD::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+2.0*M_h_e+M_q8_c"
		"R_D_LACt2::M_h_e+M_lac_D_e --> M_h_c+M_lac_D_c"
		"R_D_LACt2_reverse::M_h_c+M_lac_D_c --> M_h_e+M_lac_D_e"
		"R_ENO::M_2pg_c --> M_h2o_c+M_pep_c"
		"R_ENO_reverse::M_h2o_c+M_pep_c --> M_2pg_c"
		"R_ETOHt2r::M_etoh_e+M_h_e --> M_etoh_c+M_h_c"
		"R_ETOHt2r_reverse::M_etoh_c+M_h_c --> M_etoh_e+M_h_e"
		"R_EX_ac_e::M_ac_e --> M_ac_b"
		"R_EX_acald_e::M_acald_e --> M_acald_b"
		"R_EX_akg_e::M_akg_e --> M_akg_b"
		"R_EX_co2_e::M_co2_e --> M_co2_b"
		"R_EX_co2_e_reverse::M_co2_b --> M_co2_e"
		"R_EX_etoh_e::M_etoh_e --> M_etoh_b"
		"R_EX_for_e::M_for_e --> M_for_b"
		"R_EX_fru_e::M_fru_e --> M_fru_b"
		"R_EX_fum_e::M_fum_e --> M_fum_b"
		"R_EX_glc_e::M_glc_D_e --> M_glc_D_b"
		"R_EX_glc_e_reverse::M_glc_D_b --> M_glc_D_e"
		"R_EX_gln_L_e::M_gln_L_e --> M_gln_L_b"
		"R_EX_glu_L_e::M_glu_L_e --> M_glu_L_b"
		"R_EX_h_e::M_h_e --> M_h_b"
		"R_EX_h_e_reverse::M_h_b --> M_h_e"
		"R_EX_h2o_e::M_h2o_e --> M_h2o_b"
		"R_EX_h2o_e_reverse::M_h2o_b --> M_h2o_e"
		"R_EX_lac_D_e::M_lac_D_e --> M_lac_D_b"
		"R_EX_mal_L_e::M_mal_L_e --> M_mal_L_b"
		"R_EX_nh4_e::M_nh4_e --> M_nh4_b"
		"R_EX_nh4_e_reverse::M_nh4_b --> M_nh4_e"
		"R_EX_o2_e::M_o2_e --> M_o2_b"
		"R_EX_o2_e_reverse::M_o2_b --> M_o2_e"
		"R_EX_pi_e::M_pi_e --> M_pi_b"
		"R_EX_pi_e_reverse::M_pi_b --> M_pi_e"
		"R_EX_pyr_e::M_pyr_e --> M_pyr_b"
		"R_EX_succ_e::M_succ_e --> M_succ_b"
		"R_FBA::M_fdp_c --> M_dhap_c+M_g3p_c"
		"R_FBA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c"
		"R_FBP::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c"
		"R_FORt2::M_for_e+M_h_e --> M_for_c+M_h_c"
		"R_FORti::M_for_c --> M_for_e"
		"R_FRD7::M_fum_c+M_q8h2_c --> M_q8_c+M_succ_c"
		"R_FRUpts2::M_fru_e+M_pep_c --> M_f6p_c+M_pyr_c"
		"R_FUM::M_fum_c+M_h2o_c --> M_mal_L_c"
		"R_FUM_reverse::M_mal_L_c --> M_fum_c+M_h2o_c"
		"R_FUMt2_2::M_fum_e+2.0*M_h_e --> M_fum_c+2.0*M_h_c"
		"R_G6PDH2r::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c"
		"R_G6PDH2r_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c"
		"R_GAPD::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c"
		"R_GAPD_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c"
		"R_GLCpts::M_glc_D_e+M_pep_c --> M_g6p_c+M_pyr_c"
		"R_GLNS::M_atp_c+M_glu_L_c+M_nh4_c --> M_adp_c+M_gln_L_c+M_h_c+M_pi_c"
		"R_GLNabc::M_atp_c+M_gln_L_e+M_h2o_c --> M_adp_c+M_gln_L_c+M_pi_c"
		"R_GLUDy::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh4_c"
		"R_GLUDy_reverse::M_akg_c+M_nadph_c+M_nh4_c --> M_glu_L_c+M_h2o_c+M_nadp_c"
		"R_GLUN::M_gln_L_c+M_h2o_c+M_h_c --> M_glu_L_c+M_nh4_c"
		"R_GLUSy::M_akg_c+M_gln_L_c+M_h_c+M_nadph_c --> 2.0*M_glu_L_c+M_nadp_c"
		"R_GLUt2r::M_glu_L_e+M_h_e --> M_glu_L_c+M_h_c"
		"R_GLUt2r_reverse::M_glu_L_c+M_h_c --> M_glu_L_e+M_h_e"
		"R_GND::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c"
		"R_H2Ot::M_h2o_e --> M_h2o_c"
		"R_H2Ot_reverse::M_h2o_c --> M_h2o_e"
		"R_ICDHyr::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c"
		"R_ICDHyr_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c"
		"R_ICL::M_icit_c --> M_glx_c+M_succ_c"
		"R_LDH_D::M_lac_D_c+M_nad_c --> M_h_c+M_nadh_c+M_pyr_c"
		"R_LDH_D_reverse::M_h_c+M_nadh_c+M_pyr_c --> M_lac_D_c+M_nad_c"
		"R_MALS::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c"
		"R_MALt2_2::2.0*M_h_e+M_mal_L_e --> 2.0*M_h_c+M_mal_L_c"
		"R_MDH::M_mal_L_c+M_nad_c --> M_h_c+M_nadh_c+M_oaa_c"
		"R_MDH_reverse::M_h_c+M_nadh_c+M_oaa_c --> M_mal_L_c+M_nad_c"
		"R_ME1::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c"
		"R_ME2::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c"
		"R_NADH16::4.0*M_h_c+M_nadh_c+M_q8_c --> 3.0*M_h_e+M_nad_c+M_q8h2_c"
		"R_NADTRHD::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c"
		"R_NH4t::M_nh4_e --> M_nh4_c"
		"R_NH4t_reverse::M_nh4_c --> M_nh4_e"
		"R_O2t::M_o2_e --> M_o2_c"
		"R_O2t_reverse::M_o2_c --> M_o2_e"
		"R_PDH::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c"
		"R_PFK::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c"
		"R_PFL::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c"
		"R_PGI::M_g6p_c --> M_f6p_c"
		"R_PGI_reverse::M_f6p_c --> M_g6p_c"
		"R_PGK::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c"
		"R_PGK_reverse::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c"
		"R_PGL::M_6pgl_c+M_h2o_c --> M_6pgc_c"
		"R_PGM::M_2pg_c --> M_3pg_c"
		"R_PGM_reverse::M_3pg_c --> M_2pg_c"
		"R_PIt2r::M_h_e+M_pi_e --> M_h_c+M_pi_c"
		"R_PIt2r_reverse::M_h_c+M_pi_c --> M_h_e+M_pi_e"
		"R_PPC::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c"
		"R_PPCK::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c"
		"R_PPS::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c"
		"R_PTAr::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c"
		"R_PTAr_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c"
		"R_PYK::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c"
		"R_PYRt2r::M_h_e+M_pyr_e --> M_h_c+M_pyr_c"
		"R_PYRt2r_reverse::M_h_c+M_pyr_c --> M_h_e+M_pyr_e"
		"R_RPE::M_ru5p_D_c --> M_xu5p_D_c"
		"R_RPE_reverse::M_xu5p_D_c --> M_ru5p_D_c"
		"R_RPI::M_r5p_c --> M_ru5p_D_c"
		"R_RPI_reverse::M_ru5p_D_c --> M_r5p_c"
		"R_SUCCt2_2::2.0*M_h_e+M_succ_e --> 2.0*M_h_c+M_succ_c"
		"R_SUCCt3::M_h_e+M_succ_c --> M_h_c+M_succ_e"
		"R_SUCDi::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c"
		"R_SUCOAS::M_atp_c+M_coa_c+M_succ_c --> M_adp_c+M_pi_c+M_succoa_c"
		"R_SUCOAS_reverse::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c"
		"R_TALA::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c"
		"R_TALA_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c"
		"R_THD2::2.0*M_h_e+M_nadh_c+M_nadp_c --> 2.0*M_h_c+M_nad_c+M_nadph_c"
		"R_TKT1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c"
		"R_TKT1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c"
		"R_TKT2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c"
		"R_TKT2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c"
		"R_TPI::M_dhap_c --> M_g3p_c"
		"R_TPI_reverse::M_g3p_c --> M_dhap_c"
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"M_13dpg_c"
		"M_2pg_c"
		"M_3pg_c"
		"M_6pgc_c"
		"M_6pgl_c"
		"M_ac_c"
		"M_ac_e"
		"M_acald_c"
		"M_acald_e"
		"M_accoa_c"
		"M_acon_C_c"
		"M_actp_c"
		"M_adp_c"
		"M_akg_c"
		"M_akg_e"
		"M_amp_c"
		"M_atp_c"
		"M_cit_c"
		"M_co2_c"
		"M_co2_e"
		"M_coa_c"
		"M_dhap_c"
		"M_e4p_c"
		"M_etoh_c"
		"M_etoh_e"
		"M_f6p_c"
		"M_fdp_c"
		"M_for_c"
		"M_for_e"
		"M_fru_e"
		"M_fum_c"
		"M_fum_e"
		"M_g3p_c"
		"M_g6p_c"
		"M_glc_D_e"
		"M_gln_L_c"
		"M_gln_L_e"
		"M_glu_L_c"
		"M_glu_L_e"
		"M_glx_c"
		"M_h2o_c"
		"M_h2o_e"
		"M_h_c"
		"M_h_e"
		"M_icit_c"
		"M_lac_D_c"
		"M_lac_D_e"
		"M_mal_L_c"
		"M_mal_L_e"
		"M_nad_c"
		"M_nadh_c"
		"M_nadp_c"
		"M_nadph_c"
		"M_nh4_c"
		"M_nh4_e"
		"M_o2_c"
		"M_o2_e"
		"M_oaa_c"
		"M_pep_c"
		"M_pi_c"
		"M_pi_e"
		"M_pyr_c"
		"M_pyr_e"
		"M_q8_c"
		"M_q8h2_c"
		"M_r5p_c"
		"M_ru5p_D_c"
		"M_s7p_c"
		"M_succ_c"
		"M_succ_e"
		"M_succoa_c"
		"M_xu5p_D_c"
		"M_ac_b"
		"M_acald_b"
		"M_akg_b"
		"M_co2_b"
		"M_etoh_b"
		"M_for_b"
		"M_fru_b"
		"M_fum_b"
		"M_glc_D_b"
		"M_gln_L_b"
		"M_glu_L_b"
		"M_h2o_b"
		"M_h_b"
		"M_lac_D_b"
		"M_mal_L_b"
		"M_nh4_b"
		"M_o2_b"
		"M_pi_b"
		"M_pyr_b"
		"M_succ_b"
	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
