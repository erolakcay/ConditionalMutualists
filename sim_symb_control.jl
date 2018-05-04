
using Distributions

# Define the immutable type host to have the follow attributes:
#  - patch: integer indicating the patch the host lives in
#  - symbiont: integer indicating the host's infection status (1 = uninfected, 2 = infected)
#  - ht: *symbiont's* transmission horizontal transmission probability
#  - vt: *symbiont's* vertical transmission probability

immutable host
    patch :: Int
    symb :: Int
    ht :: Float64
    vt :: Float64
end

Base.show(io::IO, z::host) = print(io, z.patch, ", ", z.symb, ", ", z.ht, ", " , z.vt)

# Produces a function for host birth.
#  - Input: fec, the fecundity of hosts in each patch (rows) for uninfected hosts and hosts infected with each symbiont (columns)
#  - Output: the function `birth`
#      - `birth` input: pop,  current population
#      - `birth` output: a randomly chosen parent (based on fecundity)

function makeBirth(fec::Array{Float64})
    function birth(pop::Vector{host})
        fecundities = [fec[pop[i].patch, pop[i].symb] for i=1:length(pop)];
    return pop[rand(Categorical(fecundities / sum(fecundities)))]
    end
    return birth
end

# Produces a function for symbiont mutation.
#  - Input:
#      - Mu: mutation rate
#      - sdH: standard deviation of mutations in horizontal transmission probability
#      - sdV: standard deviation of mutations in vertical transmission probability
#  - Output: the function `hostMut`
#      - `hostMut` input: newborn, a newborn host
#      - `hostMut` output: if newborn infected, the input newborn with its symbiont's transmission probabilities mutated; if newborn is uninfected, it is returned unchanged

function makeSymbMut(mu::Float64, sdH::Float64, sdV::Float64)
    function symbMut(newborn::host)
        h = newborn.ht;
        v = newborn.vt;
        if !isnan(h) & (rand() < mu)
            h = max(min(h + sdH * randn(), 1), 0);
        end
        if !isnan(v) & (rand() < mu)
            v = max(min(v + sdV * randn(), 1), 0);
        end
        return host(newborn.patch, newborn.symb, h, v)
    end
    return symbMut
end

# Produces a function for host dispersal.
#  - Input: d, the dispersal rate, and patches, the total number of patches
#  - Output: the function `dispersal`
#      - `dispersal` input: parent, a host parent, and newborn, the parent's offspring
#      - `dispersal` output: the input newborn with its patch equal to its parent's patch with probability 1- d, and in a random other patch with probability d

function makeDispersal(d::Float64, patches::Int)
    function dispersal(parent::host, newborn::host)
        dispProbs = [d / (patches - 1) for i=1:patches];
        dispProbs[parent.patch] = 1 - d;
        return host(rand(Categorical(dispProbs)), newborn.symb, newborn.ht, newborn.vt)
    end
    return dispersal
end

# Function for newborn infection.
# - Input:
#      - c, the number of potentially infectious contacts a newborn has
# - Output: the function 'transmission'
#      - `transmission` input: pop, the current population of hosts (vector); parent, a host parent; and newborn, a host which is the offspring of parent
#      - `transmission` output: input newborn with infection status 1 or 2 to reflect uninfected or infected (respectively). Infection probability is determined by the parent/contact's transmission probabilities.
#           If infected, input newborn is returned with the transmission probabilities of the symbiont that infected it.

function makeTransmission(c::Int)
    function transmission(pop::Vector{host}, parent::host, newborn::host)
        if (parent.symb == 2) & (rand() < parent.vt)
            (symb, ht, vt) = (2, parent.ht, parent.vt);
        else
            patchmates = [pop[i].patch == newborn.patch for i=1:length(pop)]; # 1 if in newborn's patch, 0 otherwise
            contacts = [pop[rand(Categorical(patchmates/sum(patchmates)))] for i=1:c]; # newborn's potentially infectious contacts

            (symb, ht, vt) = (1, NaN, NaN);
            for i=1:c
              if rand() < contacts[i].ht
                (symb, ht, vt) = (2, contacts[i].ht, contacts[i].vt);
                break
              end
            end
        end
        return host(newborn.patch, symb, ht, vt)
    end
    return transmission
end



# Produces a function for newborn establishment/host death.
#  - Input:
#      - est, the establishment probability of hosts in each patch (rows) for uninfected hosts and hosts infected with each symbiont (columns)
#      - mort, the chance of dying (relative to hosts in other circumstances) of hosts in each patch (rows) for uninfected hosts and hosts infected with each symbiont (columns)
#  - Output: the function `establishment`
#      - `establishment` input: pop, the current population of hosts, and newborn, a newborn host
#      - `establishment` output: the population with the newborn established in place of an adult in its patch (with probability given by est) or the population unchanged (with probability 1 - est[newborn's patch, newborn's infection status])

function makeEstablishment(est::Array{Float64}, mort::Array{Float64})
    function establishment(pop::Vector{host}, newborn::host)
        if rand() < est[newborn.patch, newborn.symb]
            deathProbs = [(pop[i].patch == newborn.patch) * mort[pop[i].patch, pop[i].symb] for i=1:length(pop)];
            pop[rand(Categorical(deathProbs / sum(deathProbs)))] = newborn;
        end
        return pop
    end
    return establishment
end



# Function for recording data.
#  - Input:
#      - simHistory: an array of hosts present in the population at various time points (columns indicate different time points)
#      - myfile: a string giving the file data should be saved in
#      - record: the interval between recordings (1 = no interval between recordings)
#      - burnin: the interval before recording was started (0 = recording started immediately)
#  - Output: a file with the hosts in simHistory recorded with time points assigned to them based on burnin and record

function writeData(simHistory::Array{host}, myfile::String, record::Int=1, burnin::Int=0)
open(myfile, "w") do f
    write(f, "time, patch, symb, h, v \r\n");
        write(f, ["$(j*record + burnin), $(simHistory[i,j]) \r\n" for i=1:size(simHistory,1), j=1:size(simHistory,2)]);
    end
end

# Function for producing a population of hosts infected with a monomorphic population of symbionts.
#  - Input:
#      - popData: array of the number of uninfected and infected hosts in each patch. Rows correspond to patches and columns to infection status (first column = uninfected, second column = infected)
#      - ht: horizontal transmission probability of symbionts
#      - vt: vertical transmission probability of symbionts
#  - Output: a vector of hosts with the numbers of infected and uninfected hosts in each patch determined by popData

function makePop(popData::Array{Int}, ht::Float64, vt::Float64)
    return [host(i, j, ht, vt) for i=1:size(popData, 1) for j=1:size(popData, 2)  for n=1:popData[i, j]]
end

# Function to run simulation.
#  - Input:
#      - pop = starting population
#      - d = dispersal rate
#      - fec = fecundity of hosts by patch (rows) and infection status (columns)
#      - est = fecundity of hosts by patch (rows) and infection status (columns)
#      - mort = relative mortality of hosts by patch (rows) and infection status (columns)
#      - mu = mutation rate
#      - sdH = standard deviation of mutations in horizontal transmission probability
#      - sdV = standard deviation of mutations in vertical transmission probability
#      - c = number of potentially infectious contacts a newborn host has
#      - tmax = number of time steps (hosts births) to run simulation for
#      - record = population state will be recorded when time is a multiple of record
#      - burnin = time steps to wait before starting to record population state
#  - Output: simHistory, an array of hosts present at the time points given by record and burnin (columns indicate different time points)

function simSymbCntrl(pop::Vector{host}, d::Float64, fec::Array{Float64}, est::Array{Float64},
    mort::Array{Float64}, mu::Float64, sdH::Float64, sdV::Float64,
    c::Int, tmax::Int, record::Int, burnin::Int=0)

    # Calculate number of patches and make array for storing states of the population
    patches = maximum([pop[i].patch for i=1:length(pop)]);
    simHistory = Array{host}(length(pop), Int(floor((tmax - burnin)/record)));

    # Define functions
    birth = makeBirth(fec);
    mut = makeSymbMut(mu, sdH, sdV);
    dispersal = makeDispersal(d, patches);
    transmission = makeTransmission(c);
    establishment = makeEstablishment(est, mort);
    # (transmission does not depend on any external factors or mutation rate)



    for t=1:tmax

        #birth
        parent = birth(pop);
        newborn = host(0, 0, NaN, NaN);

        # dispersal
        newborn = dispersal(parent, newborn);

        #infection
        newborn = transmission(pop, parent, newborn);

        # mutation (of symbionts)
        newborn = mut(newborn);

        # establishment/death
        pop = establishment(pop, newborn);

        #store data
        if (t - burnin > 0) && (mod(t - burnin, record) == 0)
            simHistory[:,Int((t - burnin)/record)] = pop;
        end
    end
    return simHistory
end
