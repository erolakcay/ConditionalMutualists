#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script contains functions for finding the ecological equilibrium fraction of infected hosts. These functions work for both host and symbiont control of transmission evolution, since the ecological equilibrium for a monomorphic population is not affected by which partner controls transmission evolution.*)


BeginPackage[ "ecologicalEquilibrium`"]


fbar::usage = "fbar[iP, iQ, fPU, fPI, fQU, fQI] finds the average host fecundity given the fraction of infected hosts in each patch (iP and iQ) and the fecundities of all possible hosts types (fPU, fPI, fQU, fQI, where P and Q indicate patch and U and I indicate uninfected or infected hosts).";
mbar::usage = "mbar[iP, mPU, mPI] finds the average mortality in a patch given the fraction of infected hosts (iP) and the mortality of uninfected (mPU) and infected (mPI) hosts.";
diPdt::usage = "diPdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] gives the change in the fraction of infected hosts in patch P. It takes as input the fraction of infected hosts in each patch (iP and iQ), the horizontal and vertical transmission rates (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q.";
diQdt::usage = "diQdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] gives the change in the fraction of infected hosts in patch Q. It takes as input the fraction of infected hosts in each patch (iP and iQ), the horizontal and vertical transmission rates (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q.";
stableEquil::usage = "stableEquil[h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, precision: MachinePrecision] finds stable equilibrium fractions of infected hosts in each patch. It takes as input the resident horizontal (h) and vertical (v) transmission rates, the dispersal rate (d), number of infection contacts (c), as well as the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q. The precision with which to solve for the equilibrium is an optional input (precision) that defaults to MachinePrecision. It returns a list of ecological equilibria in the form {{iP -> iPvalue1, iQ -> iQvalue1}, {iP -> iPvalue2, iQ -> iQvalue2}, ...}.";
stableEquilTable::usage = "stableEquilTable[hlist, vlist, d, c, fPU, fPI, fQU, fQI, sPU, spI, sQU, sQI, mPU, mPI, mQU, mQI, precision: MachinePrecision, epsilon:1/10000] produces a list containing (1) a table of stable equilibria at all combinations of the horizontal and vertical transmission rates given in hlist and vlist and (2) a table of the horizontal and vertical transmission rates producing those equilibria. It takes as input lists of horizontal and vertical transmission rates at which to find the equilibrium (hlist and vlist, respectively), the dispersal rate (d), number of infection contacts (c), as well as the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q. It takes two optional inputs: the precision with which to solve for the equilibrium (precision, defaults to MachinePrecision), and the amount by which to try shifting the coordinates if a stable equilibrium is not found at first (epsilon, defaults to 1/10000).";


Begin["`Private`"]


(* ::Section:: *)
(*Helper functions*)


(* ::Text:: *)
(*The functions fbar and mbar finds the average host fecundity and mortality, respectively. fbar takes as input the fraction of infected hosts in each patch (iP and iQ) as well as the fecundity of uninfected and infected hosts in each patch (fPU, fPI, fQU, fQI), where P or Q indicate the patch, and U and I indicate uninfected or infected hosts, respectively. mbar takes as input the fraction of infected hosts in one patch (iP for a generic patch) and the mortality of uninfected and infected hosts in that patch (mPU and mPI).*)


fbar[iP_, iQ_, fPU_, fPI_, fQU_, fQI_] := fPU * (1 - iP) + fPI * iP + fQU * (1 - iQ) + fQI * iQ


mbar[iP_, mPU_, mPI_] := mPU * (1 - iP) + mPI * iP


(* ::Text:: *)
(*The functions diPtdtFor0s and diQdtFor0s give the change in the fraction of infected hosts in each patch per time step (time measured in units of host births) multiplied by the average fecundity. This is more efficient for solving these functions for 0, and holds as long as all fecundities and mortalities are positive (which we assume). These functions take as input the fraction of infected hosts in each patch (iP and iQ), the resident's horizontal (h) and vertical (v) transmission rates, the dispersal rate of newborn hosts (d, the probability a newborn leaves its natal patch), the number of infectious contacts a newborn has (c, which we assume to be 1 in all analyses), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, sQI) of uninfected and infected hosts. P and Q indicate the patch, and U and I indicated uninfected or infected hosts, respectively.*)


diPdtForOs::usage = "diPdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] gives the change in the fraction of infected hosts in patch P up to a factor of 1/(fbar * mbar_P), where fbar is the average host fecundity and mbar_P the average host mortality in patch P. It takes as input the fraction of infected hosts in each patch (iP and iQ), the horizontal and vertical transmission rates (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q.";
diPdtFor0s[iP_, iQ_, h_, v_,d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_]:=
((v * (1 - d) * fPI * iP     +    v * d * fQI * iQ    +
((1 - v) * (1 - d) * fPI * iP    +    (1 - v) * d * fQI * iQ    +    (1 - d) * fPU * (1 - iP)    +   d * fQU * (1 - iQ)) *
( 1 - (1 - h * iP)^c)) * sPI * (1 - iP) * mPU 
- ((1 - v) * (1 - d) * fPI * iP + (1 - v) * d * fQI * iQ + (1 - d) * fPU * (1 - iP) + d * fQU * (1 - iQ))
* (1 - h * iP)^c * sPU * iP * mPI )


diQdtForOs::usage = "diQdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] gives the change in the fraction of infected hosts in patch P up to a factor of 1/(fbar * mbar_Q), where fbar is the average host fecundity and mbar_Q the average host mortality in patch Q. It takes as input the fraction of infected hosts in each patch (iP and iQ), the horizontal and vertical transmission rates (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q.";
diQdtFor0s[iP_, iQ_, h_, v_,d_,c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_]:=
((v * (1 - d) * fQI * iQ     +    v * d * fPI * iP    +
((1 - v) * (1 - d) * fQI * iQ    +    (1 - v) * d * fPI * iP    +    (1 - d) * fQU * (1 - iQ)    +   d * fPU * (1 - iP)) *
( 1 - (1 - h * iQ)^c)) * sQI * (1 - iQ) * mQU
- ((1 - v) * (1 - d) * fQI * iQ + (1 - v) * d * fPI * iP + (1 - d) * fQU * (1 - iQ) + d * fPU * (1 - iP))
* (1 - h * iQ)^c * sQU * iQ * mQI )


(* ::Text:: *)
(*The functions diPdt and diQdt give the change in the fraction of infected hosts in each patch per time step (time in units of host births).*)


diPdt[iP_, iQ_, h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] := (1/fbar[iP, iQ, fPU, fPI, fQU, fQI]) *
  diPdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] * (1/mbar[iP, mPU, mPI])


diQdt[iP_, iQ_, h_, v_,d_,c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_]:= (1/fbar[iP, iQ, fPU, fPI, fQU, fQI]) *
diQdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] * (1/ mbar[iQ, mQU, mQI])


(* ::Text:: *)
(*myJacobian finds the Jacobian of diPdt and diQdt, for use in stability analysis.*)


myJacobian::usage = "myJacobian finds the Jacobian of diPdt and diQdt (the change in the number of infected hosts in patches P and Q). It takes as input the fraction of infected hosts in each patch (iP and iQ), the horizontal and vertical transmission rates (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), and the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q.";
myJacobian[iP_, iQ_, h_,v_,d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_]:=
 {{D[diPdt[iP, iQ, h, v, d, c,  fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], iP], 
D[diPdt[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], iQ]}, 
{D[diQdt[iP, iQ, h, v, d, c,  fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], iP], 
D[diQdt[iP, iQ, h, v, d, c,  fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], iQ]}}


(* ::Section:: *)
(*Finding Stable Equilibria*)


(* ::Text:: *)
(*equilI finds the equilibrium fraction of infected hosts (iP and iQ), in patches P and Q. It takes as input the resident horizontal (h) and vertical (v) transmission rates, the dispersal rate (d), number of infection contacts (c), fecundities, establishment probabilities, and mortalities. The precision with which to solve for the equilibrium is an optional input that defaults to MachinePrecision. It returns a list of ecological equilibria in the form {{iP -> iPvalue1, iQ -> iQvalue1}, {iP -> iPvalue2, iQ -> iQvalue2}, ...}.*)
(*Because this function solves for iP and iQ numerically, it is possible that some solutions that are really 0 or 1 are found to be slightly below or above their true value. To avoid discarding these solutions, the function checks for any solutions that are pretty close to 0 or 1, keeps these solutions, and rounds them to 0 or 1.*)


equilI::usage = "equilI[h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, precision: MachinePrecision] finds the equilibrium fraction of infected hosts (iP and iQ), in patches P and Q. It takes as input the resident horizontal (h) and vertical (v) transmission rates, the dispersal rate (d), number of infection contacts (c), as well as the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q. The precision with which to solve for the equilibrium is an optional input (precision) that defaults to MachinePrecision. It returns a list of ecological equilibria in the form {{iP -> iPvalue1, iQ -> iQvalue1}, {iP -> iPvalue2, iQ -> iQvalue2}, ...}.";
equilI[h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, precision_: MachinePrecision, chop_: 10^-10] := Module[{result}, (* result will be the equilibrium values before rounding to 0 or 1 *) 
  result = 
   (* Solve for iP and iQ at equilibrium, keep solutions where iP, iQ are within approximately [0, 1] *)
   Select[NSolve[{diPdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] == 0, diQdtFor0s[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] == 0}, {iP, iQ}, Reals, WorkingPrecision -> precision], Function[x, ((iP /. x) >= 0) && (Chop[(iP /. x) - 1, chop] <= 0) && ((iQ /. x) >= 0) && (Chop[(iQ /. x) - 1, chop] <= 0)]];
  If[Length[result] > 0,(* Round solutions close to 0 or 1, if there are any solutions *)
   Map[If[Chop[#, chop] == 0 || Chop[# - 1, chop] == 0, Round[#], #] &, {iP, iQ} /. result, {2}],
   result]]


(* ::Text:: *)
(*stableEquil finds stable equilibrium fractions of infected hosts in each patch. It first finds equilibria and then selects those equilibria for which the eigenvalues of the Jacobian are negative.*)


stableEquil[h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, precision_:MachinePrecision, chop_:10^-10]:= Select[equilI[h, v,d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, precision, chop], Function[y, Function[x, (x[[1]] <0) && (x[[2]] <0)][ 
Eigenvalues[myJacobian[iP, iQ, h, v,d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]/.{iP -> y[[1]], iQ -> y[[2]]}]]]]


(* ::Text:: *)
(*stableEquilTable produces a list containing (1) a table of stable equilibria at all combinations of the horizontal and vertical transmission rates given in hlist and vlist and (2) a table of the horizontal and vertical transmission rates producing those equilibria. If a point does not have any stable equilibria, the function checks if a slight adjustment in the horizontal and vertical transmission rates produces a stable equilibrium. If it does, it returns that equilibrium instead. If similar horizontal and vertical transmission rates also do not lead to a stable equilibrium, the function returns an empty list for that combination of transmission rates.*)
(**)
(*Output is returned in the form of a table where rows correspond to the horizontal transmission rate and columns to the vertical transmission rate, in the same order as in hlist and vlist. Entries in the table are lists of {iP -> iP_at_equilibrium, iQ -> iQ_at_equilibrium} (where iP and iQ are the fraction of infected hosts in patches P and Q). In some cases there may be more or less than one stable equilibrium for a given horizontal and vertical transmission rate pair, which is why the entries are lists of equilibria.*)
(**)
(*The function takes as input lists of horizontal and vertical transmission rates at which to find the equilibrium (hlist and vlist, respectively), the dispersal rate, number of infectious contacts, fecundites, and establishment probabilities. It also takes two optional inputs: precision, the precision with which to numerically solve for the equilibrium, and epsilon, the distance to shift the transmission rates when looking for a nearby equilibrium for a pair of transmission rates that has no stable equilibrium.*)


stableEquilTable::usage = "stableEquilTable[hlist, vlist, d, c, fPU, fPI, fQU, fQI, sPU, spI, sQU, sQI, mPU, mPI, mQU, mQI, precision: MachinePrecision, epsilon:1/10000] produces a list containing (1) a table of stable equilibria at all combinations of the horizontal and vertical transmission rates given in hlist and vlist and (2) a table of the horizontal and vertical transmission rates producing those equilibria. It takes as input lists of horizontal and vertical transmission rates at which to find the equilibrium (hlist and vlist, respectively), the dispersal rate (d), number of infection contacts (c), as well as the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q. It takes two optional inputs: the precision with which to solve for the equilibrium (precision, defaults to MachinePrecision), and the amount by which to try shifting the coordinates if a stable equilibrium is not found at first (epsilon, defaults to 1/10000).";
stableEquilTable::NoStableEquil = "No stable equilibrium found near point. Try a different epsilon.";
stableEquilTable[hlist_, vlist_, d_, c_, fPU_, fPI_,fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, precision_:MachinePrecision, epsilon_:1/10000, chop_:10^-10] := 
Module[ {equil, nEquil, maxEquil, minEquil, coordEquil}, 

(* Make a table of stable equilibria (equil) and calculate the number of stable equilibria at each entry (nequil) *)
equil=Table[stableEquil[h, v,d, c, fPU, fPI,fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, precision, chop], {h, hlist}, {v, vlist}];
nEquil = Map[Length, equil, {2}];
{minEquil, maxEquil} = {Min[nEquil], Max[nEquil]}; (* Minimum and maximum number of stable equilibria, respectively *)

(* For each (h, v) pair in equil, write the coordinates of the equilibrium/equilibria *)
coordEquil= Table[{h, v}, {h, hlist}, {v, vlist}];

(* Re-write the entries of equil with no stable equilibria with approximations based on similar transmission rates; adjust coodinates accordingly *)
If[minEquil == 0, (* If some points have no stable ecological equilibria, replace them with nearby points' equilibria *)
equil =Table[
If[nEquil[[hind, vind]] == 0, (* If the point does not have a stable equilibrium, try shifting the transmission rates slightly *)
Module[{equilTemp, h = hlist[[hind]], v = vlist[[vind]], hTemp, vTemp},
hTemp = If[h <= epsilon, h+ epsilon, h- epsilon]; (* Nearby horizontal transmission rate; it will be h - epsilon unless that would be \[LessEqual] 0, in which case it will be h + epsilon *)
vTemp = If[v <= epsilon, v + epsilon, v - epsilon]; (* Nearby vertical transmission rate; it will be v - epsilon unless that would be \[LessEqual] 0, in which case it will be v + epsilon *) 
 equilTemp = stableEquil[hTemp, vTemp, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, precision, chop]; (* stable equilibrium for hTemp and vTemp *)
If[Length[equilTemp] == 0, (* If the new points still don't have an equilibrium, try once more, shifting the coordinates by slightly different amounts *)
hTemp = If[ h <= Sqrt[1/2] epsilon, h + Sqrt[1/2] epsilon, h - Sqrt[1/2] epsilon]; (* Nearby horizontal transmission rate *)
vTemp = If[v <= Sqrt[3/2] epsilon, v + Sqrt[3/2] epsilon, v - Sqrt[3/2] epsilon]; (* Nearby vertical transmission rate *)
equilTemp = stableEquil[hTemp, vTemp, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI,precision, chop];]
If[Length[equilTemp] == 0, Message[stableEquilTable::NoStableEquil]]; (* If equilTemp still doesn't have any stable equilibrium, print a warning *)
coordEquil[[hind, vind]] = {hTemp, vTemp}; (* Set hTemp, vTemp as the new coordinates of the equilibrium *)
equilTemp  (* Set equilTemp as the new entry in equil *)],
equil[[hind, vind]]], 
{hind, Length[hlist]}, {vind, Length[vlist]}]];

(* Sort the equilibria for each pair of transmission rates *)
equil = Map[Sort, equil, {2}];
Return[{equil, coordEquil}]]


End[]
EndPackage[]
