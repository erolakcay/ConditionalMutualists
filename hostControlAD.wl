(* ::Package:: *)

(* ::Text:: *)
(*This file is `hostControlAD.wl`.*)


(* ::Text:: *)
(*This script contains functions for calculating the direction of transmission mode evolution when hosts control transmission.*)


BeginPackage["hostControlAD`", {"ecologicalEquilibrium`", "adaptiveDynamicsFunctions`"}]


dMutGrowthHost::usage = "Host control version: dMutGrowthHost[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] estimates the partial derivative of the mutant host growth rate in terms of the horizontal and vertical transmission rates evaluated at the resident's transmission rates. It returns a list of the partial derivative in terms of the horizontal transmission rate (1st entry) and the partial derivative in terms of the vertical transmission rate (2nd entry). The function takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the resident's horizontal and vertical transmission rates (h and v), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds).";
dMutGrowthTableHost::usage = "Host control version: dMutGrowthTableHost[equil, coordEquil, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] produces a table of derivatives of mutant host growth rates, evaluated at the points given in coordEquil. It takes as input a table of resident equilibria (equil), a table of resident horizontal and vertical transmission rates at which those equilibria occur (coordEquil), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds). If a pair of transmission rates have multiple stable ecological equilibria, dMutGrowthTable will return a list of derivatives (e.g. {{ddh_equil1, ddv_equil1}, {ddh_equil2, ddv_equil2}}).";


Begin["`Private`"]


(* ::Text:: *)
(*birthHost and deathHost return matrices of mutant birth and death rates, where time is measured in units of (t N), where t is time in units of host births. The rows correspond to the location and infection status of the mutant being born or dying (the 1st row indicates establishment/death of an uninfected host in patch P, the 2nd an infected host in P, the 3rd an uninfected host in Q, and the 4th an infected host in Q). The columns correspond to the location and infection status of the mutant giving birth or dying (the death matrix only has entries along the diagonal). The columns are in the same order as the rows. Only newborns who successfully establish are counted as births.*)
(*These functions take as input the fraction of infected hosts at the resident's equilibrium (iP and iQ), the mutant's (hmut and vmut, for birth) or resident's (hres and vres, for death) horizontal and vertical transmission rates, and the dispersal rate, number of infectious contacts, fecundities, and establishment probabilities (and mortalities, for death).*)


birthHost::usage = "Host control: birthHost[iP, iQ, hmut, vmut, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI] returns a matrix of mutant host birth rates. It takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the mutant's horizontal and vertical transmission rates (hmut and vmut), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI) and establishment probabilities (sPU, sPI, sQU, sQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q." ;
birthHost[iP_, iQ_, hmut_, vmut_, d_,c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_] := (1/ fbar[iP, iQ, fPU, fPI, fQU, fQI]) * 
{{fPU * (1 - d) * (1 - hmut*iP)^c * sPU ,
fPI * (1 - d) * (1 - vmut) * (1 - hmut*iP)^c * sPU,
fQU * d * (1 - hmut*iP)^c * sPU,
fQI * d * (1 - vmut) * (1 - hmut*iP)^c * sPU},
{fPU * (1 - d) * (1 - (1 - hmut*iP)^c ) * sPI,
fPI * (1- d) * (vmut + (1 - vmut)*(1 - (1 - hmut*iP)^c)) * sPI,
fQU * d * (1 - (1 - hmut*iP)^c) * sPI,
fQI * d * (vmut + (1 - vmut)* (1 - (1 - hmut*iP)^c)) * sPI},
{fPU * d * (1 - hmut*iQ)^c * sQU,
fPI * d * (1 - vmut) * (1 - hmut*iQ)^c * sQU,
fQU * (1 - d) * (1 - hmut*iQ)^c * sQU,
fQI * (1 - d) * (1 - vmut) * (1 - hmut*iQ)^c * sQU},
{fPU * d * (1 - (1 - hmut*iQ)^c) * sQI,
fPI * d * (vmut + (1 - vmut) * (1 - (1 - hmut*iQ)^c)) * sQI,
fQU * (1 - d) * (1 - (1 - hmut*iQ)^c)*sQI,
fQI * (1 - d) * (vmut + (1 - vmut) * (1 - (1 - hmut*iQ)^c)) * sQI}}


deathHost::usage = "Host control version: death[iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] returns a matrix of mutant host death rates. It takes as input the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the resident's horizontal and vertical transmission rates (hres and vres), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
deathHost[iP_, iQ_ , hres_, vres_, d_,c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_]:=  
DiagonalMatrix[Flatten[{{1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 1, 1}} . birthHost[iP, iQ, hres, vres, d,c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI] . {{1 - iP}, {iP}, {1 - iQ}, {iQ}}]] . DiagonalMatrix[{mPU/mbar[iP, mPU, mPI], mPI/mbar[iP, mPU, mPI], mQU/mbar[iQ, mQU, mQI], mQI/mbar[iQ, mQU, mQI]}]


(* ::Text:: *)
(*mutXHost returns a matrix of the expected number of mutants produced (time in units of (t N)) by a mutant host when it is at low frequency, taking into account that production may be negative if the host dies. If the leading eigenvalue of mutXHost is > 1, the mutant should invade. The rows correspond to the location and infection status of the mutant being produced, and the columns correspond to the location and infection status of the individual producing the new mutant. The 1st row and column correspond to an uninfected host in patch P,  the 2nd to an infected host in patch P, the 3rd to an uinfected host in patch Q, and the fourth to an infected host in patch Q. mutXHost takes as input the fraction of infected hosts at the resident's equilibrium (iP and iQ), the mutant and resident horizontal and vertical transmission rates (hmut and vmut for the mutant, hres and vres for the resident), the dispersal rate, number of infectious contacts, fecundities, establishment probabilities, and mortalities.*)


mutXHost::usage = "Host control version: mutX[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, mPU, mPI, mQU, mQI] returns a matrix of the expected number of mutants produced by a mutant host when it is at low frequency. It takes as input the fraction of infected hosts in each patch (iP and iQ), the mutant and resident horizontal and vertical transmission rates (hmut and vmut for the mutant, hres and vres for the resident), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q.";
mutXHost[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=MatrixExp[birthHost[iP, iQ , hmut, vmut, d,c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI] - deathHost[iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]


dMutGrowthHost[iP_, iQ_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_,mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] := adaptiveDynamicsFunctions`Private`dMutGrowth[mutXHost, iP, iQ, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta, timeOut]


dMutGrowthTableHost[equil_,coordEquil_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] := adaptiveDynamicsFunctions`Private`dMutGrowthTable[mutXHost, equil,coordEquil, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta, timeOut]


End[]
EndPackage[]
