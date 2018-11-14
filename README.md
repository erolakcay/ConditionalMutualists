# ConditionalMutualists
Code for Brown and Akcay, 2018, Evolution of transmission mode in conditional mutualisms with spatial variation

### To run simulations:
Run `HostControlMain_N200.jl`, `HostControlMain_N1000.jl`, `HostControlMain_3Patches.jl`, or  `SymbiontControlMain.jl` in Julia. These simulate transmission evolution in the following populations:
1. `HostControlMain_N200.jl`: transmission evolution in hosts; 2 environmental patches with 100 hosts each
2. `HostControlMain_N1000.jl`: transmisssion evolution in hosts; 2 environmental patches with 500 hosts each
3. `HostControlMain_3Patches.jl`: transmision evolution in hosts; 3 environmental patches with 100 hosts each (there are unequal numbers of patches where the symbiont is beneficial and patches where it is harmful)
4. `SymbiontControlMain.jl`: transmission evolution in hosts; 2 environmental patches with 100 hosts each

To work, these scripts must be in the same directory as `sim_host_control.jl` (simulations 1-3 in the list above) or `sim_symb_control.jl` (simulation 4).

### To summarize simulation results:
Run `summarize_simulations.py` in Python v2.7.

### To run the analytical model and plot simulations results:
Run `run_analytical_model_and_all_plots.wls` in Mathematica. This script runs in parallel. To change the number of cores used, change the number inside `LaunchKernels[#]` on line 9. This script can also be run in sections from within Mathematica, once the first section (lines  1-12) has been run to load the packages. This script requires `adaptiveDynamics.wl`, `ecologicalEquilibrium.wl`, `hostcontrolAD.wl`, `symbiontControlAD.wl`, and `plottingFunctions.wl` to be in the same directory. In order to plot simulation results, `summarize_simulations.py` must have already been run on the simulation results.

### Software versions used:
Julia v0.5
Python v2.7.15
Mathematica v11
