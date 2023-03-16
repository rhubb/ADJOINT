This repository contains supporting code for "Identifying sources of disparities in surveillance mammography performance and personalized  recommendations for supplemental breast imaging: A simulation study" by Hubbard, Pujol, Alhajjar, Edoh, and Martin

Code in this repository generates simulation inputs based on Breast Cancer Surveillance Consortium (BCSC) data, reads in simulation inputs, and executes simulations as described in the paper. Code should be run in the following order:

0_bcsc_simulation_prepare_data.R - this file requires BCSC data, outputs simulation parameters to inputs directory
1_bcsc_simulation_inputs.R - reads in inputs from inputs directory
2_bcsc_simulation_funcn.R - contains functions that execute the simulations
3_bcsc_simulation_run.R - runs simulations as described in the manuscript

Simulations can be run based on inputs in the inputs directory without requiring raw BCSC data by omitting 0_bcsc_simulation_prepare_data.R. 

BCSC data is available by request to the authors and BCSC Steering Committee, with appropriate regulatory approvals.
