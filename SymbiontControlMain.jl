# This file is `SymbiontControlMain.jl`.
# Purpose of this is to simulation transmission evolution in a finite population of symbionts
# This code is meant to run many simulations in parallel -- run with
#    "julia -p # SymbiontControlMain.jl", where # is the number of cores to use
# To run a single simulation, use `myfun` (below) without `pmap`, or
#    use the function `simSymbCntrl` from sim_symb_control.jl

# Parameters
# - Host transmission evolution
# - Population initialized at 100% infection
# - Population size = 200 (100 in each of 2 patches)
# - Time to equilibrate before evolution = 4000 time steps
# - Time to evolve = 10^7 time steps
# - Mutation rate = 0.02
# - Mutation standard deviation = 0.05
# Note: symbiont evolution simulations do not allow spontaneous infection of
# newborn hosts not infected vertically or horizontally

# Add the functions for simulating host evolution to every core
@everywhere include("sim_symb_control.jl")

# Create lists of the parameter values for each simulation. The idea is to have
# all possible combinations of parameters for the following:
# - flist = fecundity, slist = establishment probability, mlist = mortality;
#   unlike all the other parameters, the values for the symbiont effects only
#   show up in certain combinations: f = 0.5, s = 1.0, m = 1.0 (symbiont affects fecundity)
#   or f = 1.0, s = 0.5, m = 1.0 (symbiont affects newborn establishment = lifespan)
# - dlist = dispersal rates
# - nlist = population size
# - patchlist = number of patches (note: the maximum dispersal rate is 1/(# patches))
# - hlist = initial horizontal transmission probabilities
# - vlist = initial vertical transmssion probabilities
# - filelist = names of the files to save the simulation results in

@everywhere flist = [f for f=(0.5, 1.0) for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere slist = [s for s=(1.0, 0.5) for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere mlist = [m for m=(1.0, 1.0) for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere dlist = Float64[d for fs=1:2 for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere nlist = Int[n for fs=1:2 for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere patchlist = Int[p for fs=1:2 for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere hlist = Float64[h for fs=1:2 for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere vlist = Float64[v for fs=1:2 for d=(0.005, 0.05, 0.5) for n=(200) for p=(2) for h=0:0.1:1 for v=0:0.1:1]
@everywhere filelist = ["Simulation results/Symbiont_Control/f = $(flist[i]), s = $(slist[i]), m = $(mlist[i])/d = $(dlist[i])/N = $(nlist[i]), patches = $(patchlist[i])/hinit = $(hlist[i]), vinit = $(vlist[i]).csv" for i=1:length(dlist)]


# This function initializes and runs a single simulation
# Input:
# - n: population size
# - patches: number of patches
# - startH: starting horizontal transmission probability
# - startV: starting vertical transmission probability
# - d: dispersal rate
# - f: fecundity cost of infection (P-patches) or lack of infection (M-patches)
# - s: establishment cost of infection (P-patches) or lack of infection (M-patches)
# - m: mortality decrease due to infection (M-patches) or lack of infection (P-patches)
# - myfile: file to save results to
# Output:
# - simulation results saved to myfile

@everywhere function myfun(n::Int, patches::Int, startH::Float64, startV::Float64, d::Float64, f::Float64, s::Float64, m::Float64, myfile::String)
	# if the directories that contain myfile do not exist, make them
	if !isdir(dirname(myfile))
		mkpath(dirname(myfile))
	end

	# initialize the population so all hosts are infected
	initPop = makePop(ones(Int, patches) * Int[0 n/patches], startH, startV);

	# create matrices of the effects of infection; rows = patches, columns = infection statuses
	# fec = fecundity, est = establishment, mort = mortality
	fec = [ones(Float64, Int(patches/2)) * Float64[1 f]; ones(Float64, Int(patches/2)) * Float64[f 1]];
	est = [ones(Float64, Int(patches/2)) * Float64[1 s]; ones(Float64, Int(patches/2)) * Float64[s 1]];
	mort = [ones(Float64, Int(patches/2)) * Float64[m 1]; ones(Float64, Int(patches/2)) * Float64[1 m]];

	# allow the population to reach an ecological equilibrium (no evolution)
	ecoEquil = simSymbCntrl(initPop, d, fec, est, mort,
		0.0, 				# mutation rate = 0
		0.0, 0.0, 	# standard deviation of horiz. and vert. transmission probability mutations = 0
		1, 					# number of potentially infectious contacts = 1
		4000, 			# equilibration time = 4000 time steps
		1, 3999); 	# record last time step (record interval = 1, burnin = 3999)

	# simulate transmission evolution, initializing with population (hopefully) at ecological equilibrium

	return writeData(simSymbCntrl(vec(ecoEquil), d, fec, est, mort,
		0.02, 			# mutation rate = 0.02
		0.05, 0.05, # standard deviation of horiz. and vert. transmission probability mutations = 0.05
		1, 					# number of potentially infectious contacts = 1
		10000000, 	# simulation run time = 10^7 time steps
		20000, 0), 	# record every 20,000th time step (record interval = 20000, burnin = 0)
		myfile, 20000, 0)
end

# Map the function to every parameter value (pmap does this in parallel)
pmap((n, patches, startH, startV, d, f, s, m, myfile) -> myfun(n, patches, startH, startV, d, f, s, m, myfile),
	nlist, patchlist, hlist, vlist, dlist, flist, slist, mlist, filelist)
