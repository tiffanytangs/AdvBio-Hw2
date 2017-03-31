# Script to estimate the acetate and celmass in W3110 under AEROBIC conditions
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("include.jl")

# setup the time-scale -
time_start = 0.0
time_stop = 10.0
time_step = 0.1
time_array = collect(time_start:time_step:time_stop)
number_of_timesteps = length(time_array)

# Fire up the max cellmass -
data_dictionary = maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 10.5 #mmol/g/h
K_glucose_uptake = 1.0 #mmol/L
vmax_acetate_uptake=10.5 #mmol/g/h
K_acetate_uptake=1.0 #mmol/L

# initialize the problem -
number_of_external_states =3
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.5     # 2 acetate
state_array[1,3] = 0.001   # 3 cellmass

# flux storage
length_flux=142 #number of reactions
length_time=length(time_array)
mega_flux_array=zeros(length_flux,length_time)

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1
  #check glucose levels
  if state_array[time_step_index,1] < 1
    copy_data_dictionary = deepcopy(data_dictionary)

    # grab the state -
    glucose = state_array[time_step_index,1]
    acetate = state_array[time_step_index,2]
    cellmass = state_array[time_step_index,3]

    # calculate acetate uptake -
    qAct = vmax_acetate_uptake*(acetate)/(K_acetate_uptake+acetate) #mmol/g/h

    # setup the species bounds to use acetate for cell growth
    species_bounds_array = copy_data_dictionary["species_bounds_array"]
    species_bounds_array[81,1] = 0.0
    species_bounds_array[81,2] = 0.0
    species_bounds_array[73,1] = -qAct
    species_bounds_array[73,2] = -0.99*qAct

    # calculate the fluxes using the LP -
    (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)
    #stores the flux_array into a giant flux_array
    mega_flux_array[:,time_step_index]=flux_array

    # grab the growth rate from the flux_array -
    mu = 1.23*flux_array[24]

    # update the external state -
    state_array[time_step_index+1,1] = 0
    state_array[time_step_index+1,2] = acetate - qAct*cellmass*time_step #mMol
    state_array[time_step_index+1,3] = cellmass + mu*cellmass*time_step #g/L

  else
    #Consuming glucose

    # make a deepcopy of the data_dictionary -
    copy_data_dictionary = deepcopy(data_dictionary)

    # grab the state -
    glucose = state_array[time_step_index,1]
    acetate = state_array[time_step_index,2]
    cellmass = state_array[time_step_index,3]

    # calculate glucose uptake -
    qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

    # setup the species bounds -
    species_bounds_array = copy_data_dictionary["species_bounds_array"]
    species_bounds_array[81,1] = -qGlc
    species_bounds_array[81,2] = -0.99*qGlc

    # calculate the fluxes using the LP -
    (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)
    #stores flux_array in a giant matrix
    mega_flux_array[:,time_step_index]=flux_array

    # grab the respective fluxes from the flux_array -
    mu = 1.23*flux_array[24]
    qAcetate_production = flux_array[35]

    # update the external state -
    state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass*time_step
    state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass*time_step
    state_array[time_step_index+1,3] = cellmass + mu*cellmass*time_step
  end

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)
end

#Save the mega flux array
file_path = "./mega_flux_array_aerobic.dat"
writedlm(file_path, mega_flux_array)

#Plotting Goodness
using PyPlot
#biomass experimental values - Figure 7
x_b=[1.78,3.43,5.18, 6.32, 7.22, 8.78, 9.709];
y_b=[.006,.03, .15, .31, .5967, .7213, .739];

#acetate experimental values - Figure 7
x_a=[1.733,3.06,4.17,5.23,6.27,7.39,8.266,9.324,10.31];
y_a=[0.339,.528,.8592,1.38,2.27,3.46,3.25,1.11,.537];

#glucose experimental values - Figure 7
x_g=[1.709,3.019,4.298,5.375,6.77,7.829,9.285,10.194];
y_g=[10.626,10.009,9.54,8.19,4.32,.54,0,0];


#Plots experimental data using a scatterplot
subplot(2,2,1)
scatter(x_g,y_g,label="Glucose Expt", color="lightgreen")
plot(time_array, state_array[:,1], label="Glucose FBA", color="green")
xlabel("Time [hr]")
ylabel("Glucose Concentration [mMol]")
legend();

subplot(2,2,2)
scatter(x_a,y_a,label="Acetate Expt", color="skyblue")
plot(time_array, state_array[:,2], label="Acetate FBA", color="skyblue")
xlabel("Time [hr]")
ylabel("Acetate Concentration [mMol]")
legend();

subplot(2,2,3)
scatter(x_b,y_b,label="Biomass Expt", color="navajowhite")
plot(time_array, state_array[:,3], label="Biomass FBA", color="orange")
xlabel("Time [hr]")
ylabel("Cell Density [g/L]")
legend();

savefig("./figs/Aerobic_profiles.pdf")
