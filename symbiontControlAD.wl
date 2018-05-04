(* ::Package:: *)

(* ::Text:: *)
(*This script contains functions for calculating the direction of transmission mode evolution when symbionts control transmission.*)


BeginPackage["symbiontControlAD`", {"ecologicalEquilibrium`", "adaptiveDynamicsFunctions`"}]


dMutGrowthSymbiont::usage = "Symbiont control version: dMutGrowthSymbiont[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] estimates the partial derivative of the mutant symbiont growth rate in terms of the horizontal and vertical transmission rates evaluated at the resident's transmission rates. It returns a list of the partial derivative in terms of the horizontal transmission rate (1st entry) and the partial derivative in terms of the vertical transmission rate (2nd entry). The function takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the resident's horizontal and vertical transmission rates (h and v), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds).";
dMutGrowthTableSymbiont::usage = "Symbiont control version: dMutGrowthTableSymbiont[equil, coordEquil, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] produces a table of derivatives of mutant symbiont growth rates, evaluated at the points given in coordEquil. It takes as input a table of resident equilibria (equil), a table of resident horizontal and vertical transmission rates at which those equilibria occur (coordEquil), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds). If a pair of transmission rates have multiple stable ecological equilibria, dMutGrowthTable will return a list of derivatives (e.g. {{ddh_equil1, ddv_equil1}, {ddh_equil2, ddv_equil2}}).";


Begin["`Private`"]


(* ::Text:: *)
(*These functions give the probability that a host will arrive in a given patch (P or Q) with a given infection status (U for uninfected, I for infected) prior to horizontal transmission. This is helpful for calculating a mutant's ability to horizontally infect hosts, and it is also useful for calculating death rates. All four functions take as input the fraction of infected hosts in both patches (iP and iQ), the resident's vertical transmission rate (vres), the dispersal rate (d), and fecundities.*)


arrivalPU::usage = "arrivalPU[iP, iQ, vres, d, fPU, fPI, fQU, fQI] gives the probability an uninfected newborn arrives in patch P. It takes as input the fraction of infected hosts in each patch (iP and iQ), the resident's vertical transmission rate (vres), the host dispersal rate (d), and the fecundities (fPU, fPI, fQU, fQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
arrivalPU[iP_, iQ_, vres_, d_, fPU_, fPI_, fQU_, fQI_] :=
(1/fbar[iP, iQ, fPU, fPI, fQU, fQI] )* ((1 - d) * (fPU * (1 - iP) + (1 - vres) * fPI * iP) + d * (fQU * (1 - iQ) + (1 - vres) * fQI * iQ))


arrivalPI::usage = "arrivalPI[iP, iQ, vres, d, fPU, fPI, fQU, fQI] gives the probability an infected newborn arrives in patch P. It takes as input the fraction of infected hosts in each patch (iP and iQ), the resident's vertical transmission rate (vres), the host dispersal rate (d), and the fecundities (fPU, fPI, fQU, fQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
arrivalPI[iP_, iQ_,  vres_, d_, fPU_, fPI_, fQU_, fQI_] :=
(1/fbar[iP, iQ, fPU, fPI, fQU, fQI] )* ((1 - d) * vres * fPI * iP + d * vres * fQI * iQ)


arrivalQU::usage = "arrivalPI[iP, iQ, vres, d, fPU, fPI, fQU, fQI] gives the probability an uninfected newborn arrives in patch Q. It takes as input the fraction of infected hosts in each patch (iP and iQ), the resident's vertical transmission rate (vres), the host dispersal rate (d), and the fecundities (fPU, fPI, fQU, fQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
arrivalQU[iP_, iQ_, vres_, d_, fPU_, fPI_, fQU_, fQI_] :=
(1/fbar[iP, iQ, fPU, fPI, fQU, fQI] )* (d * (fPU * (1 - iP) + (1 - vres) * fPI * iP) + (1 - d) * (fQU * (1 - iQ) + (1 - vres) * fQI * iQ))


arrivalQU::usage = "arrivalPI[iP, iQ, vres, d, fPU, fPI, fQU, fQI] gives the probability an infected newborn arrives in patch Q. It takes as input the fraction of infected hosts in each patch (iP and iQ), the resident's vertical transmission rate (vres), the host dispersal rate (d), and the fecundities (fPU, fPI, fQU, fQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
arrivalQI[iP_, iQ_,  vres_, d_, fPU_, fPI_, fQU_, fQI_] :=
(1/fbar[iP, iQ, fPU, fPI, fQU, fQI] )* (d * vres * fPI * iP + (1 -d) *vres * fQI * iQ)


(* ::Text:: *)
(*birthSymbiont and deathSymbiont return matrices of mutant birth and death rates, where time is measured in units of (t N), where t is time in units of host births. The rows correspond to the location of the mutant being born or dying (the 1st row indicates establishment/death in patch P, the 2nd in patch Q). The columns correspond to the location of the mutant giving birth or dying (the death matrix only has entries along the diagonal). The columns are in the same order as the rows. Only newborns whose hosts successfully establish are counted as births.*)
(*These functions take as input the fraction of infected hosts at the resident's equilibrium (iP and iQ), (for birth only) the mutant's horizontal and vertical transmission rates (hmut and vmut), the resident's horizontal and vertical transmission rates (hres and vres), and the dispersal rate, number of infectious contacts, fecundities, and establishment probabilities (and mortalities, for death).*)


birthSymbiont::usage = "Symbiont control version: birthSymbiont[iP, iQ, hmut, vmut, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] returns a matrix of mutant symbiont death rates. It takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the mutant's horizontal and vertical transmission rates (hmut and vmut), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI) and establishment probabilities (sPU, sPI, sQU, sQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q." ;
birthSymbiont[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_] :=
{{((1 - d) * vmut * fPI/fbar[iP, iQ, fPU, fPI, fQU, fQI]  + arrivalPU[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * hmut * ( 1 + Sum[(1 - hres iP)^j, {j, c-1}])) * sPI,
d * vmut * fQI * sPI/fbar[iP, iQ, fPU, fPI, fQU, fQI]},
{d * vmut * fPI * sQI/fbar[iP, iQ, fPU, fPI, fQU, fQI],
((1- d) * vmut * fQI/fbar[iP, iQ, fPU, fPI, fQU, fQI] + arrivalQU[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * hmut * (1 + Sum[(1 - hres iQ)^j, {j, c-1}])) * sQI}}


deathSymbiont::usage = "Symbiont control version: deathSymbiont[iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] returns a matrix of mutant symbiont death rates. It takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the resident's horizontal and vertical transmission rates (hres and vres), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
deathSymbiont[iP_, iQ_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=
{{(arrivalPU[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * ((1 - hres * iP) ^ c * sPU + (1 - (1 - hres * iP)^c) * sPI) + arrivalPI[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * sPI) * mPI/mbar[iP, mPU, mPI],
0},
{0,
(arrivalQU[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * ((1 - hres * iQ) ^ c * sQU + (1 - (1 - hres * iQ)^c) * sQI) + arrivalQI[iP, iQ, vres, d, fPU, fPI, fQU, fQI] * sQI) * mQU/mbar[iQ, mQU, mQI]}}


(* ::Text:: *)
(*mutXSymbiont returns a matrix of the expected number of mutants produced (time in units of (t N)) by a mutant host when it is at low frequency, taking into account that production may be negative if the host dies. If the leading eigenvalue of mutXSymbiont is > 1, the mutant should invade. The rows correspond to the location and infection status of the mutant being produced, and the columns correspond to the location and infection status of the individual producing the new mutant. The 1st row and column correspond to an uninfected host in patch P,  the 2nd to an infected host in patch P, the 3rd to an uinfected host in patch Q, and the fourth to an infected host in patch Q. mutX takes as input the fraction of infected hosts at the resident's equilibrium (iP and iQ), the mutant and resident horizontal and vertical transmission rates (hmut and vmut for the mutant, hres and vres for the resident), the dispersal rate, number of infectious contacts, fecundities, establishment probabilities, and mortalities.*)


mutXSymbiont::usage = "Symbiont control version: mutXSymbiont[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, mPU, mPI, mQU, mQI] returns a matrix of the expected number of mutants produced by a mutant host when it is at low frequency. It takes as input the fraction of infected hosts in each patch (iP and iQ), the mutant and resident horizontal and vertical transmission rates (hmut and vmut for the mutant, hres and vres for the resident), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
mutXSymbiont[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=MatrixExp[birthSymbiont[iP, iQ , hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI] - deathSymbiont[iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]


dMutGrowthSymbiont[iP_, iQ_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_,mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] := adaptiveDynamicsFunctions`Private`dMutGrowth[mutXSymbiont, iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta, timeOut]


dMutGrowthTableSymbiont[equil_,coordEquil_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] := adaptiveDynamicsFunctions`Private`dMutGrowthTable[mutXSymbiont, equil,coordEquil, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta, timeOut]


End[]
EndPackage[]
