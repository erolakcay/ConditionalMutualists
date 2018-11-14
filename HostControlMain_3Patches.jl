# This file is `HostControlMain_3Patches.jl`.
# Purpose of this is to simulation transmission evolution in a finite population
# of hosts where one type of environment is more common than the other.
# This code is meant to run many simulations in parallel -- run with
#    "julia -p # HostControlMain_3Patches.jl", where # is the number of cores to use
# To run a single simulation, use `myfun` (below) without `pmap`, or
#    use the function `simHostCntrl` from sim_host_control.jl

# Parameters
# - Host transmission evolution
# - Population initialized at 100% infection
# - Population size = 300 (100 in each of 3 patches)
# - Time to equilibrate before evolution = 6000 time steps
# - Time to evolve = 1.5*10^7 time steps
# - Chance of spontaneous infection = 0.0033
# - Mutation rate = 0.02
# - Mutation standard deviation = 0.05

# Add the functions for simulating host evolution to every core
@everywhere include("sim_host_control.jl")

# Create lists of the parameter values for each simulation. The idea is to have
# all possible combinations of parameters for the following:
# - flist = fecundity, slist = establishment probability, mlist = mortality;
#   unlike all the other parameters, the values for the symbiont effects only
#   show up in certain combinations: f = 0.5, s = 1.0, m = 1.0 (symbiont affects fecundity)
#   or f = 1.0, s = 0.5, m = 1.0 (symbiont affects newborn establishment = lifespan)
# - dlist = dispersal rates (note: the maximum dispersal rate is 1/(# patches))
# - nlist = population size
# - pMlist = whether there are 2 (true) or 1 (false) patches where the symbiont is mutualistic
# - hlist = initial horizontal transmission probabilities
# - vlist = initial vertical transmssion probabilities
# - replist = an identifier for each replicate of the simulation
# - filelist = names of the files to save the simulation results in

@everywhere flist = [f for r=1:2 for f=(0.5, 1.0) for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere slist = [s for r=1:2 for s=(1.0, 0.5) for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere mlist = [m for r=1:2 for m=(1.0, 1.0) for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere dlist = Float64[d for r=1:2 for fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere nlist = Int[n for r=1:2 for fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere pMlist = Bool[pM for r=1:2 for fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere hlist = Float64[h for r=1:2 for fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere vlist = Float64[v for r=1:2 for fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere replist = [r for r=("a", "b") for
	fsm=1:2 for d=((1/150), (1/15), (2/3)) for n=(300) for pM=(true, false) for h=0:0.1:1 for v=0:0.1:1]
@everywhere filelist = ["Simulation results/Host_Control_3_patches_rep_$(replist[i])/f = $(flist[i]), s = $(slist[i]), m = $(mlist[i])/d = $(round(dlist[i], 4))/N = $(nlist[i]), $(["P-patches 1 & 2, M-patch 3" "P-patch 1, M-patches 2 & 3"][pMlist[i] + 1])/hinit = $(hlist[i]), vinit = $(vlist[i]).csv" for i=1:length(dlist)]

# This function initializes and runs a single simulation
# Input:
# - n: population size
# - twoMpatches: boolean for whether patch 2 is the same as patch 1 or patch 3.
#      If true, the symbiont will be parasitic in patch 1 and mutualistic in
#      patches 2 & 3. If false, the symbiont will be parasitic in patches 1 & 2
#      and mutualistic in patch 3. Assumes f, m <= 1 (s must be <= 1 for sim to
#      run). If f, m > 1, symbiont will decrease fecundity & increase
#      mortality in "M" patches
# - startH: starting horizontal transmission probability
# - startV: starting vertical transmission probability
# - d: dispersal rate
# - f: fecundity cost of infection (P-patches) or lack of infection (M-patches)
# - s: establishment cost of infection (P-patches) or lack of infection (M-patches)
# - m: mortality decrease due to infection (M-patches) or lack of infection (P-patches)
# - myfile: file to save results to
# Output:
# - simulation results saved to myfile

@everywhere function myfun(n::Int, twoMpatches::Bool, startH::Float64, startV::Float64,
	d::Float64, f::Float64, s::Float64, m::Float64, myfile::String)
	# if the directories that contain myfile do not exist, make them
	if !isdir(dirname(myfile))
		mkpath(dirname(myfile))
	end

	# initialize the population so all hosts are infected
	initPop = makePop(ones(Int, 3) * Int[0 n/3], startH, startV);

	# create matrices of the effects of infection; rows = patches,
	# columns = infection statuses (1st = uninfected, 2nd = infected)
	# fec = fecundity, est = establishment, mort = mortality
	if twoMpatches # make the symbiont parasitic in patch 1 and mutualistic
		# in patches 2 & 3 (assumes f, s, m <= 1)
		fec  = [1 f; f 1; f 1];
		est  = [1 s; s 1; s 1];
		mort = [m 1; 1 m; 1 m];
	else # make the symbiont parasitic in patches 1 & 2 and mutualistic in
		# patch 3 (assumes f, s, m <= 1)
		fec  = [1 f; 1 f; f 1];
		est  = [1 s; 1 f; s 1];
		mort = [m 1; m 1; 1 m];
	end

	# allow the population to reach an ecological equilibrium (no evolution)
	ecoEquil = simHostCntrl(initPop, d, fec, est, mort,
		0.0, 				# mutation rate = 0
		0.0, 0.0, 	# standard deviation of horiz. and vert. transmission probability mutations = 0
		0.0033, 			# spontaneous infection probability = 0.005
		1, 					# number of potentially infectious contacts = 1
		6000, 			# equilibration time = 4000 time steps
		1, 5999); 	# record last time step (record interval = 1, burnin = 3999)

	# simulate transmission evolution, initializing with population (hopefully) at ecological equilibrium

	return writeData(simHostCntrl(vec(ecoEquil), d, fec, est, mort,
		0.02, 			# mutation rate = 0.02
		0.05, 0.05, # standard deviation of horiz. and vert. transmission probability mutations = 0.05
		0.0033, 			# spontaneous infection probability = 0.005
		1, 					# number of potentially infectious contacts = 1
		15000000, 	# simulation run time = 1.5*10^7 time steps
		30000, 0), 	# record every 30,000th time step (record interval = 30000, burnin = 0)
		myfile, 30000, 0)
end

# Map the function to every parameter value (pmap does this in parallel)
pmap((n, twoMpatches, startH, startV, d, f, s, m, myfile) -> myfun(n, twoMpatches, startH, startV, d, f, s, m, myfile),
	nlist, pMlist, hlist, vlist, dlist, flist, slist, mlist, filelist)
