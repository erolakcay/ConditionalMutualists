(* ::Package:: *)

BeginPackage["adaptiveDynamicsFunctions`"]


dMutGrowthAvgTable::usage = "dMutGrowthAvgTable[dGrowth] averages multiple growth rate derivatives that are present for a given set of transmission rates. It takes as input a table of growth rate derivatives, with each entry having the form {{ddhValue1, ddvValue1}, ... }, where ddhValue and ddvValue are the derivatives of the growth rate with respect to the horizontal and vertical transmission rates.";
trajectoryEndPoint::usage = "trajectoryEndPoint[hgrowth, vgrowth, dt, hinit, vinit, stoppingThreshold:10^-10] finds the endpoint of an evolutionary trajectory starting at horizontal transmission rate hinit and vertical transmission rate vinit. Besides the initial point, it takes as input the derivatives of the horizontal and vertical transmission rates (hgrowth and vgrowth), which should be given as functions of the resident horizontal and vertical transmission rates. trajectoryEndPoint as takes as input dt, the distance over which to forecast the trajectory. trajectoryEndPoint has a final, optional argument, stoppingThreshold. If two successive steps produce transmission rate pairs which are less than stoppingDistance apart, the function will consider the trajectory to be at its endpoint. stoppingThreshold defaults to 10^-10.";
getEvoTrajectories::usage = "getEvoTrajectories[equil, coord, dgrowth, hStartPoints, vStartPoints, dt, stoppingThreshold: 10^-10, groupingThreshold: 10^-5] finds the endpoints of evolutionary trajectories in a system with thegiven ecological equilibria (equil) located at the transmission rates given by (coord) and with the derivatives of the mutant growth rates at these points givenby (dgrowth). It investigates the trajectories starting at a grid of transmission rates given by hStartPoints and vStartPoints. It also takes as input the distance (dt)over which to forecast changes in transmission rate based on the derivatives. It takes stoppingThreshold and groupingThreshold as optional input. If two successive steps produce transmission rate pairs which are less than stoppingThreshold apart, the function will consider the trajectory to be at its endpoint. stoppingThreshold defaults to 10^-10. Similarly, if the start and endpoints of a trajectory are less than groupingThreshold apart, then the trajectory will be labeled as starting and ending at the same point. groupingThreshold defaults to 10^-5.
For each start point, getEvolTrajectories returns a line a list containing {the start point, the endpoint, the initial fraction of infected hosts, the final fraction of infected hosts, notes on the trajectory}.";


Begin["`Private`"]


(* ::Section:: *)
(*Derivatives of Mutant Growth Rates*)


(* ::Text:: *)
(*This function estimates the partial derivative of the mutant growth rate in terms of the horizontal and vertical transmission rates evaluated at the resident's transmission rates. It returns a list of the partial derivative in terms of the horizontal transmission rate (1st entry) and the partial derivative in terms of the vertical transmission rate (2nd entry). The function takes as input a function giving a matrix of the expected number of mutants of each type produced by each type of mutant (mutX[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]),  the equilibrium fraction of infected hosts at the resident's growth rate, the resident's transmission rates, the dispersal rate, number of infectious contacts, fecundities, establishment probabilities, mortalities, the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds).*)


dMutGrowth::usage = "dMutGrowth[iP, iQ, h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] estimates the partial derivative of the mutant growth rate in terms of the horizontal and vertical transmission rates evaluated at the resident's transmission rates. It returns a list of the partial derivative in terms of the horizontal transmission rate (1st entry) and the partial derivative in terms of the vertical transmission rate (2nd entry). The function takes as input a function giving a matrix of the expected number of mutants of each type produced by each type of mutant (mutX[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]), the equilibrium fraction of infected hosts at the resident's growth rate (iP and iQ for patches P and Q), the resident's horizontal and vertical transmission rates (h and v), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds).";
dMutGrowth::TimedOut = "Eigenvalue calculation timed out. Used numerical approximation of mutX to get eigenvalues.";
dMutGrowth[mutX_, iP_, iQ_, hres_, vres_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_,mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] :=
Module[{dGdh, dGdv}, (* dGdh and dGdv will be the derivative of mutant growth rate in terms of the horizontal and vertical transmission rates, respectively *)
If[(* If no hosts are infected in either patch, h & v are selectively neutral. Growth rate does not change with changes in h or v. *)
iP == 0 && iQ == 0, {dGdh, dGdv} = {0, 0},
(* Otherwise, numerically calculate the growth rate. *)
If[(* If hres - delta \[LessEqual] 0, calculate the partial derivative of growth rate in terms of h using the mutant growth at hres and at hres + delta *)
hres <= delta, dGdh = (TimeConstrained[Eigenvalues[mutX[iP, iQ, hres + delta, vres, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], 1], timeOut,Message[dMutGrowth::TimedOut]; Eigenvalues[N[mutX[iP, iQ, hres + delta, vres, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]], 1]][[1]] - 1)/delta,
(* Otherwise, calculate the partial derivative using the mutant growth rate at hres and hres - delta *)
dGdh = (1 - TimeConstrained[Eigenvalues[mutX[iP, iQ, hres - delta, vres, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], 1], timeOut,Message[dMutGrowth::TimedOut]; Eigenvalues[N[mutX[iP, iQ, hres - delta, vres, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]], 1]] [[1]])/delta];
If[(* If vres - delta \[LessEqual] 0, calculate the partial derivative of growth rate in terms of v using the mutant growth at vres and at vres + delta *)
vres <= delta,  dGdv = (TimeConstrained[Eigenvalues[mutX[iP, iQ, hres, vres + delta, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], 1], timeOut, Message[dMutGrowth::TimedOut]; Eigenvalues[N[mutX[iP, iQ, hres, vres + delta, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]], 1]][[1]] - 1)/delta,
(* Otherwise, calculate the partial derivative using the mutant growth rate at vres and vres - delta *)
dGdv = (1 - TimeConstrained[Eigenvalues[mutX[iP, iQ, hres, vres - delta, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], 1], timeOut, Message[dMutGrowth::TimedOut]; Eigenvalues[mutX[iP, iQ, hres, vres - delta, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], 1]][[1]])/delta];];
Return[{dGdh, dGdv}]]


(* ::Text:: *)
(*dMutGrowthTable produces a table of derivatives of mutant growth rates, evaluated at the points given in coordEquil. It takes as input a function giving a matrix of the expected number of mutants of each type produced by each type of mutant (mutX[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]), a  table of resident equilibria (equil), a table of resident horizontal and vertical transmission rates at which those equilibria occur (coordEquil), the dispersal rate, number of infectious contacts, fecundites, establishment probabilities, mortalities, the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds). If a pair of transmission rates have multiple stable ecological equilibria, dMutGrowthTable will return a list of derivatives (e.g. {{ddh_equil1, ddv_equil1}, {ddh_equil2, ddv_equil2}}). *)


dMutGrowthTable::usage = "dMutGrowthTable[mutX, equil, coordEquil, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta:1/1000, timeOut:5] produces a table of derivatives of mutant growth rates, evaluated at the points given in coordEquil. It takes as input a function giving a matrix of the expected number of mutants of each type produced by each type of mutant (mutX[iP, iQ, hmut, vmut, hres, vres, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]), a table of resident equilibria (equil), a table of resident horizontal and vertical transmission rates at which those equilibria occur (coordEquil), the dispersal rate (d), number of infectious contacts (c), fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) hosts in patches P and Q. It also takes as optional inputs the change in transmission rate over which to calculate the derivative (delta, defaults to 1/1000), and the length of time (in seconds) to let the calculation go before switching to a less accurate one (timeOut, defaults to 5 seconds). If a pair of transmission rates have multiple stable ecological equilibria, dMutGrowthTable will return a list of derivatives (e.g. {{ddh_equil1, ddv_equil1}, {ddh_equil2, ddv_equil2}}).";
dMutGrowthTable[mutX_, equil_,coordEquil_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, delta_:1/1000, timeOut_:5] := 
Table[Map[dMutGrowth[mutX, #[[1]], #[[2]], coordEquil[[hind, vind, 1]], coordEquil[[hind, vind, 2]], d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, delta, timeOut]&, equil[[hind, vind]]], {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}]


(* ::Text:: *)
(*dMutGrowthAvgTable averages multiple growth rate derivatives that are present for a given set of transmission rates. It takes as input a table of growth rate derivatives, with each entry having the form {{ddhValue1, ddvValue1}, ... }, where ddhValue and ddvValue are the derivatives of the growth rate with respect to the horizontal and vertical transmission rates. dMutGrowthAvgTable returns a table of derivatives where, the ddhValues and ddvValues at each entry have been average to produce {{ddhAvg,  ddvAvg}}.*)


dMutGrowthAvgTable[dGrowth_]:= Map[If[Length[#] >= 1, {{Mean[Map[First,#]], Mean[Map[Last, #]]}}, {}]&, dGrowth, {2}]


(* ::Section:: *)
(*Basins of Attraction*)


(* ::Text:: *)
(*trajectoryShortDist and trajectoryEndPoint forecast an evolutionary trajectory a short distance (trajectoryShortDist) or until it stops moving (trajectoryEndPoint). They take as input functions giving the derivatives of the horizontal and vertical transmission rates (hgrowth and vgrowth) as functions of the resident transmission rate, the distance over which to forecast the trajectory (dt, trajectoryShortDist stops after one forecast, trajectoryEndPoint then performs another forecast based on the new resident transmission rates), and the initial horizontal and vertical transmission rates (hinit and vinit). trajectoryEndPoint also takes an optional argument, stoppingThreshold, the difference between the successive forecasts below which the trajectory is considered to have reached its endpoint. stoppingThreshold defaults to 10^-10.*)


trajectoryShortDist::usage = "trajectoryShortDist[hgrowth, vgrowth, dt, hinit, vinit] forecasts transmission mode evolution a short distance, starting from initial horizontal and vertical transmission rates hinit and vinit. It takes as input the derivatives of the horizontal and vertical transmission rates (hgrowth and vgrowth,  these should be functions of the resident horizontal and vertical transmission rates), the distance over which to forecast the trajectory (dt), and the initial horizontal and vertical transmission rates (hinit and vinit).";
trajectoryShortDist[hgrowth_, vgrowth_, dt_, hinit_, vinit_] := {Max[Min[hgrowth[hinit, vinit] * dt + hinit, 1], 0], Max[Min[vgrowth[hinit, vinit] * dt+ vinit, 1], 0]}


trajectoryEndPoint[hgrowth_, vgrowth_, dt_, hinit_, vinit_, stoppingThreshold_:10^-10] := 
FixedPoint[trajectoryShortDist[hgrowth, vgrowth, dt, #[[1]], #[[2]]]&, {hinit, vinit}, SameTest->(EuclideanDistance[#1, #2] <stoppingThreshold &)]


(* ::Text:: *)
(*getEvoTrajectories finds the endpoints of evolutionary trajectories in a system with the given ecological equilibria (equil) located at the transmission rates given by (coord) and with the derivatives of the mutant growth rates at these points givenby (dgrowth). It investigates the trajectories starting at a grid of transmission rates given by hStartPoints and vStartPoints. It also takes as input the distance (dt)over which to forecast changes in transmission rate based on the derivatives. It takes stoppingThreshold and groupingThreshold as optional input. If two successive steps produce transmission rate pairs which are less than stoppingThreshold apart, the function will consider the trajectory to be at its endpoint. stoppingThreshold defaults to 10^-10. Similarly, if the start and endpoints of a trajectory are less than groupingThreshold apart, then the trajectory will be labeled as starting and ending at the same point. groupingThreshold defaults to 10^-5. For each start point, getEvolTrajectories returns a line a list containing {the start point, the endpoint, the initial fraction of infected hosts, the final fraction of infected hosts, notes on the trajectory}.*)


getEvoTrajectories[equil_, coord_, dgrowth_, hStartPoints_, vStartPoints_, dt_ , stoppingThreshold_: 10^-10, groupingThreshold_:10^-5] := Module[{nEquil, minEquil, maxEquil, analysisOrder},
  (* Determine the number of analyses needed to investigate all equilibria and produce an ordering of those analyses *)
  nEquil = Map[Length, equil, {2}]; (* Number of equilibria at each horizontal transmission rate/vertical transmission rate pair *)
  {minEquil, maxEquil} = {Min[nEquil], Max[nEquil]}; (* Minimum (minEquil) and maximum (maxEquil) equilibria any point has *)
  analysisOrder = Tuples[Table[Range[n], {n, minEquil, maxEquil}]]; (* Order in which to calculate evolutionary trajectories. 
  First level (rows) of the list corresponds to the analysis. Second level (columns) corresponds to the number of equilibria a point has. 
  plotOrder[[row, column]] gives the equilibrium for the rowth analysis for all points with column # of equilibria. *)
  
  (* Find evolutionary trajectories for each combination of equilibria *)
  Table[Module[{tempEquil, tempdGrowth, iPInterpolation, iQInterpolation, dhInterpolation, dvInterpolation, evoEndPoints},
    (* Select the equilibria (tempEquil) and growth rates (tempdGrowth) to find evolutionary trajectories for *)
    tempEquil = Map[#[[analysisOrder[[analysis, Length[#] - minEquil + 1]]]] &, equil, {2}]; tempdGrowth = Map[#[[analysisOrder[[analysis, Length[#]]]]] &, dgrowth, {2}];
    
    (* Interpolate the equilibria and growth rates *)
    iPInterpolation = Interpolation[Flatten[Table[{{N[coord[[i, j, 1]]], N[coord[[i, j, 2]]]}, N[tempEquil[[i, j, 1]]]}, {i, Dimensions[coord][[1]]}, {j, Dimensions[coord][[2]]}], 1], InterpolationOrder -> 1];
    iQInterpolation = Interpolation[Flatten[Table[{{N[coord[[i, j, 1]]], N[coord[[i, j, 2]]]}, N[tempEquil[[i, j, 2]]]},  {i, Dimensions[coord][[1]]}, {j, Dimensions[coord][[2]]}], 1], InterpolationOrder -> 1];
    dhInterpolation = Interpolation[Flatten[Table[{{N[coord[[i, j, 1]]], N[coord[[i, j, 2]]]}, N[tempdGrowth[[i, j, 1]]]}, {i, Dimensions[coord][[1]]}, {j, Dimensions[coord][[2]]}], 1], InterpolationOrder -> 1];
    dvInterpolation = Interpolation[Flatten[Table[{{N[coord[[i, j, 1]]], N[coord[[i, j, 2]]]}, N[tempdGrowth[[i, j, 2]]]}, {i, Dimensions[coord][[1]]}, {j, Dimensions[coord][[2]]}], 1], InterpolationOrder -> 1];
     
    (* Find trajectory endpoints for trajectories starting at the points given by hStartPoints and vStartPoints. Do not find trajectories for start points where there are no infected hosts. These are all stable equilibria that will be marked as NA. *)
    evoEndPoints = Table[Module[{tempEndpt}, 
       If[iPInterpolation[hinit, vinit] == 0 && iQInterpolation[hinit, vinit] == 0, 
        {{hinit, vinit}, {hinit, vinit}, {0, 0}, {0, 0}, "No initial infection"},
        
        tempEndpt = trajectoryEndPoint[dhInterpolation, dvInterpolation, dt, hinit, vinit, stoppingThreshold];
        {{hinit, vinit}, tempEndpt, {iPInterpolation[hinit, vinit], iQInterpolation[hinit, vinit]}, {iPInterpolation[tempEndpt[[1]], tempEndpt[[2]]], iQInterpolation[tempEndpt[[1]], tempEndpt[[2]]]}, If[EuclideanDistance[tempEndpt, {hinit, vinit}] < groupingThreshold, "Trajectory stayed at start point", ""]}]], 
      
      {hinit, hStartPoints}, {vinit, vStartPoints}];
    Join[{{"Startpoint", "Endpoint", "Start infection", "End infection", "Info"}}, Flatten[evoEndPoints, 1]]],
   {analysis, Length[analysisOrder]}]]


End[]
EndPackage[]
