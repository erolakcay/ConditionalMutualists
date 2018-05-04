# ConditionalMutualists
Code for Brown and Akcay, 2018, Evolution of transmission mode in conditional mutualisms with spatial variation

### To run simulations:
Run `HostControlMain_N200.jl`, `HostcontrolMain_N1000.jl`, or  `SymbiontControlMain.jl` in Julia. These simulate transmission evolving in hosts with population sizes 200 and 1000 and in symbionts with population size 200, respectively. Simulations rely on the functions in `sim_host_control.jl` and `sim_symb_control.jl`.

### To summarize simulation results:
Run `summarize_simulations_transmission.R` and `summarize_simulations_infection.R` in R.
These return the average transmission rate in M- and P-patches at each time step and the fraction of infected hosts, respectively. These scripts should be run from the command line, with folders of simulations as input. Run either script from the command line with no input to get a more detailed description of how to enter input.

### To run the analytical model and plot simulations results:
Run `run_analytical_model_and_all_plots.wls` in Mathematica. This script runs in parallel -- to change the number of cores used, change the number inside `LaunchKernels[#]` on line 9. This script can also be run in sections from within Mathematica. This script requires `adaptiveDynamics.wl`, `ecologicalEquilibrium.wl`, `hostcontrolAD.wl`, `symbiontControlAD.wl`, and `plottingFunctions.wl` (just need to be in the same directory as `run_analytical_model_and_all_plots.wls`).
