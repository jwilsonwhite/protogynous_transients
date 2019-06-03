# protogynous_transients
Model code for the transient population dynamics of sex-changing fish. Accompanies Easter et al. (in review)

The code can simulate two different sets of population dynamics: demographically 'open' with constant larval recruitment, and demographically 'closed' in which all recruitment is spawned within the model population. There is parallel code to produce the two sets of results.

All code was originally created in and run using Matlab R2018. No additional toolboxes are necessary.

File descriptions:
Transient_params.m - writes parameters values to several .mat files that are called by other functions.

Transient_Struct.m - a wrapper function that loops over several model scenarios and calls Transient_Model to perform simulations. Results are stored in a structure variable that is written to a .mat file.

Transient_Model.m - runs the closed population dynamics model, returns population trajectory and several statistics. Designed to be called from Transient_Struct.

Transient_Struct_Open.m - a wrapper function that loops over several model scenarios and calls Transient_Model_Open to perform simulations. Results are stored in a structure variable that is written to a .mat file.

Transient_Model_Open.m - runs the open population dynamics model, returns population trajectory and several statistics. Designed to be called from Transient_Struct_Open.

beta_fert_curve.m - Creates Figure 1

Fig2_open_pop_trajectories.m - reads the file created by Transient_Struct_Open and plots various results

Fig3_closed_pop_trajectories.m - reads the file created by Transient_Struct and plots various results


