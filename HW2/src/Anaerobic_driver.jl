# Script to estimate the acetate and celmass in W3110 in ANAEROBIC conditions
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
vmax_glucose_uptake = 18.5 #mmol/g/h
K_glucose_uptake = 1.0 #mmol/L
vmax_acetate_uptake=10.5 #mmol/g/h
K_acetate_uptake=1.0 #mmol/L

#Flux array storage
length_flux=142 #number of reactions
length_time=length(time_array)
mega_flux_array=zeros(length_flux,length_time)

# initialize the problem -
number_of_external_states = 5
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.001   # 2 acetate
state_array[1,3] = 0.001   # 3 cellmass
state_array[1,4] = 0.001    # 4 formate
state_array[1,5] = 0.001     # 5 ethanol
# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1

  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  glucose = state_array[time_step_index,1]
  acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,3]
  formate = state_array[time_step_index,4]
  ethanol = state_array[time_step_index,5]

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

  # setup the species bounds -
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  species_bounds_array[81,1] = -qGlc #glucose
  species_bounds_array[81,2] = -0.99*qGlc
  species_bounds_array[73,2] = 20.0 #acetate
  species_bounds_array[77,2] = 15.0 #ethnaol
  species_bounds_array[78,2] = 100.0 #formate`
  species_bounds_array[86,2] = 100.0 #lactate
  species_bounds_array[89,:] = [0.0 0.0] #oxygen

  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)
#stores the flux
  mega_flux_array[:,time_step_index]=flux_array

  # grab the fluxes from the flux_array -
  mu = 1.23*flux_array[24]
  qAcetate_production = flux_array[35]
  qFormate_production = flux_array[41]
  qEthanol_production = flux_array[40]

  # update the external state -
  state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass*time_step
  state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass*time_step
  state_array[time_step_index+1,3] = cellmass + mu*cellmass*time_step
  state_array[time_step_index+1,4] = acetate + qFormate_production*cellmass*time_step
  state_array[time_step_index+1,5] = acetate + qEthanol_production*cellmass*time_step

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)
end

#Saving the flux array
file_path = "./mega_flux_array_anaerobic.dat"
writedlm(file_path, mega_flux_array)

#Experimental points
#glucose
x_g=[1.046153846
2.994871795
4.574358974
6.112820513
7.569230769
8.594871795
9.641025641
10.60512821
];
y_g=[10.96549188
10.16500153
10.05994484
8.903340484
7.149616917
4.776953724
2.380263561
0.103585657
];
#acetate points
x_a=[3.984999657
5.001495917
6.40270363
7.559267138
8.512317299
9.403691759
10.53921636
];
y_a=[0.647773279
0.870445344
1.943319838
3.299595142
4.71659919
6.457489879
8.157894737
];
#biomass points
x_b=[0.0433,2.011,4.537,5.95, 7.49, 9.018,9.899,10.6695];
y_b=[0.008178,0.01143,0.026,0.05,0.0991,0.177,0.2357,0.241];
#formate points
x_f=[4.036166365
5.077757685
6.054249548
7.05244123
8.159132007
8.983725136
10.04701627
10.5244123
];
y_f=[0.536362922
1.061875934
2.282805252
4.025552323
5.594779464
10.85863669
13.34075006
15.12508845
];
#ethanol points
x_e=[4.009784201
4.981027341
6.012741889
7.022420959
7.949340828
8.976429608
9.98064186
10.53691153
];
y_e=[0.758710112
0.98503634
1.520223719
2.315603558
3.111095536
4.540598126
6.392895941
6.847426741
];

using PyPlot
subplot(3,2,1)
scatter(x_g,y_g,label="Glucose Expt", color="lightgreen")
plot(time_array, state_array[:,1], label="Glucose FBA", color="green")
xlabel("Time [hr]")
ylabel("Glucose Concentration [mMol]")
legend();

subplot(3,2,2)
scatter(x_a,y_a,label="Acetate Expt", color="skyblue")
plot(time_array, state_array[:,2], label="Acetate FBA", color="blue")
xlabel("Time [hr]")
ylabel("Acetate Concentration [mMol]")
legend();

subplot(3,2,3)
scatter(x_b,y_b,label="Cell Density Expt", color="navajowhite")
plot(time_array, state_array[:,3], label="Cell Density FBA", color="orange")
xlabel("Time [hr]")
ylabel("Cell Density [g/L]")
legend();

subplot(3,2,4)
scatter(x_f,y_f,label="Formate Expt", color="lightcoral")
plot(time_array, state_array[:,4], label="Formate FBA", color="crimson")
xlabel("Time [hr]")
ylabel("Formate Concentration [mMol]")
legend();

subplot(3,2,5)
scatter(x_e,y_e,label="Ethanol Expt", color="thistle")
plot(time_array, state_array[:,5], label="Ethanol FBA", color="orchid")
xlabel("Time [hr]")
ylabel("Ethanol Concentration [mMol]")
legend();

# savefig("./figs/Anaerobic_profiles.pdf")
