(* ::Package:: *)

(* ::Text:: *)
(*This script contains functions for plotting and summarizing data.*)


BeginPackage["plottingFunctions`"]


myFontStyle;
myEquilLegend;
mySimLegend;
plotAD::usage = "plotAD[equil, coordEquil, dGrowth] plots the equilibrium fraction of infected hosts in each patch, overlaid with a stream plot showing the direction of transmission mode evolution. It takes as input matrices of equilibria (equil, entries in the format {iP ->  iPvalue, iQ -> iQvalue}), coordinates of equilibria (coordEquil, entries in the format {hres -> hresvalue, vres -> vresvalue}), and mutant growth rates (dGrowth, entries in the format {ddh -> ddhvalue, ddv -> ddvvalue}). Because a single point may have multiple equilibria, this function produces plots for each equilibrium, assuming that points with the same number of equilibria have them ordered so that point1's equilibrium 1 should be plotted with point2's equilibrium 1, and point1's equilibrium 2 should be plotted with point2's equilibrium 2.";
plotBasinsAttraction::usage = "plotBasinsAttraction[essTable] plots the size of basins of attraction vs. dispersal rate. It takes as input essTable, a table of the fraction of trajectories leading to each ESS. The first column of essTable should have the names of the ESSs in the 2nd through last rows. The first row should have the dispersal rates (as numbers) in the second through last columns. Otherwise, entry[[i, j]] should give the number of trajectories leading to essTable[[i, 1]] when the dispersal rate = essTable[[1, j]]. This function also takes as optional input any option to ListLogLinearPlot.";
plotSimMulti::usage = "plotSimMulti[filesPH, filesMH, filesPV, filesMV] plots a heatmap of the data from multiple simulations with the same parameter values. It takes as input four lists of file names: filesPH and filesMH (files of average horizontal transmission rates in the patch(es) where the symbiont is a parasite and a mutualist, respectively) and filesPV and filesMV (the same as filesPH and filesMH except for vertical transmission rates).The output is a pair of heat maps showing the distribution of transmission rates at the last time point and a statement of when the last time point is for the first simulation in filesPH.";
plotSimInfectionMulti::usage = "plotSimInfectionMulti plots the fraction of infected hosts at the final horizontal and vertical transmission of the input simulation files. It takes as input six lists of file names:  filesPH and filesMH (files of average horizontal transmission rates in the patch(es) where the symbiont is a parasite and a mutualist, respectively), filesPV and filesMV (the same as filesPH and filesMH except for vertical transmission rates), and filesPInf and filesMInf (files of the fraction of infected hosts in the patch(es) where the symbiont is a parasite and a mutualist, respectively). The output is a pair of scatter plots, one for each patch, where the coordinates indicate the transmission rates  and the colors of the points indicate the fraction of infected hosts.";
myTally::usage = "myTally[list, tallyBy, nameBy] tallies the contents of a list using the function tallyBy to determine if two entries are equal. It performs the function nameBy on the elements corresponding to each tally, to allow useful descriptions of the items tallied.";
makeESStable::usage = "makeESStable[dlist] makes a table ready to be filled with counts of trajectories leading to evolutionary stable strategies. It takes as input a list of dispersal rates where transmission evolution was investigated.";
addESS::usage = "addESS[essTable, d, ess, count] adds a count of trajectories leading an ESSs to an existing table made by makeESStable. It takes as input a table of counts of evolutionary stable strategies at different dispersal rates (essTable), a dispersal rate (d, must be present in essTable), the ESS (ess), and the number of trajectories leading to the ESS (count). If the ESS already has trajectories at that dispersal rate, count will be added to the existing number of trajectories.";
talliesToESS::usage = "talliesToESS[tallies] takes a list of tallies of evolutionary endpoints at different dispersal rates and produces a table of the number of trajectories leading to each endpoint at each dispersal rate. tallies should be a list with entries of the form {d, tal} where d is the dispersal rate and tal is a list of tallies (produced by myTallies) of trajectories. tal may be of length > 1 when there are multiple possible stable ecological equilibria at pair of transmission rates.";


Begin["`Private`"]


(* ::Section:: *)
(*Fonts and Legends*)


(* ::Text:: *)
(*Color functions*)


myLineColors::usage = "A color function for streamlines. Takes as input the norm of {ddh, ddv}, but returns black no matter what.";
myLineColors[f_] := RGBColor[0, 0, 0]


myColors::usage = " A color function for the fraction of infected hosts. Takes as input the fraction of infected hosts and returns a color.";
myColors[f_] := If[f == 0, RGBColor[1, 1, 1], RGBColor[(85/255)-(85/255)f,(221/255)-(85/255)f, (255/255)-(85/255)f]]


myColorsSim::usage = "A color function for the fraction of simulation trajectories ending at a given point. Takes as input a fraction of trajectories and returns a colors. The function gives separate colors for fractions of trajectories < 30%. All points with more than 30% of trajectories ending there will have the same color.";
myColorsSim[f_] :=Module[{f2}, f2 = f/0.3; If[f <= 0.3, RGBColor[(85/255)-(85/255)f2,(221/255)-(85/255)f2, (255/255)-(85/255)f2], RGBColor[0/255, 68/255, 85/255]]]


(* ::Text:: *)
(*Font style*)


myFontStyle = {FontFamily -> "TimesNewRoman", FontColor -> Black, FontSize -> 14};


(* ::Text:: *)
(*Legend for fraction of infected hosts (based on myColors)*)


myEquilLegend = BarLegend[{myColors[#]&,{0,1}}, LegendLayout -> "Row",LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 200];


(* ::Text:: *)
(*Legend for fraction of simulations (based on myColorsSim)*)


mySimLegend = BarLegend[{myColorsSim[#]&,{0,1}}, LegendLayout -> "Row", LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 200];


(* ::Section:: *)
(*Analytical Model Plots*)


(* ::Text:: *)
(*This function plots the equilibrium fraction of infected hosts in each patch, overlaid with a stream plot showing the direction of transmission mode evolution. It takes as input matrices of equilibria (equil, entries in the format {iP ->  iPvalue, iQ -> iQvalue}), coordinates of equilibria (coordEquil, entries in the format {hres -> hresvalue, vres -> vresvalue}), and mutant growth rates (dGrowth, entries in the format {ddh -> ddhvalue, ddv -> ddvvalue}). Because a single point may have multiple equilibria, this function produces plots for each equilibrium, assuming that points with the same number of equilibria have them ordered so that point1's equilibrium 1 should be plotted with point2's equilibrium 1, and point1's equilibrium 2 should be plotted with point2's equilibrium 2.*)


plotAD[equil_, coordEquil_, dGrowth_] := 
Module[{nEquil, minEquil, maxEquil, plotOrder},
(* Determine the number of plots needed to represent all equilibria and produce an ordering of those plots *)
nEquil = Map[Length, equil, {2}]; (* Number of equilibria at each horizontal transmission rate/vertical transmission rate pair *)
{minEquil, maxEquil} = {Min[nEquil], Max[nEquil]}; (* Minimum (minEquil) and maximum (maxEquil) equilibria any point has *)
plotOrder = Tuples[Table[Range[n], {n, minEquil, maxEquil}]]; (* Order in which to plot equilibria. First level (rows) of the list corresponds to the plot. Second level (columns) corresponds to the number of equilibria a point has. plotOrder[[row, column]] gives the equilibrium to plot for the rowth plot for all points with column # of equilibria. *)

(* Make plots for each combination of equilibria *)
Table[Module[{tempEquil, tempdGrowth, iPplot, iQplot, dPlot, d0h, d0v, d0Plot},
(* Select the equilibria (tempEquil) and growth rates (tempdGrowth) to plot *)
tempEquil = Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, equil, {2}] ;
tempdGrowth = Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, dGrowth, {2}] ;

(* Plot the equilibrium fraction of infected hosts in patches P (iPplot) and Q (iQplot) *)
iPplot = ListDensityPlot[Flatten[Table[{(coordEquil[[hind, vind, 1]]), (coordEquil[[hind, vind, 2]]), (tempEquil[[hind, vind, 1]])}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}], 1],PlotRange -> {{0, 1},{0,1}},  LabelStyle -> myFontStyle, ColorFunction -> myColors, ColorFunctionScaling -> False, ClippingStyle -> {myColors[0], myColors[1]},FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}, RotateLabel -> True];(*PlotLabel \[Rule] "Patch P", PlotLegends \[Rule] Automatic, FrameLabel \[Rule] {"Horiz. Trans. Rate", "Vert. Trans. Rate"}, , RotateLabel \[Rule] True*)
iQplot = ListDensityPlot[Flatten[Table[{(coordEquil[[hind, vind, 1]]), (coordEquil[[hind, vind, 2]]), (tempEquil[[hind, vind, 2]])}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}], 1], PlotRange -> {{0, 1},{0,1}},  LabelStyle -> myFontStyle, ColorFunction -> myColors, ColorFunctionScaling -> False,  ClippingStyle -> {myColors[0], myColors[1]},FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}, RotateLabel -> True];
(*PlotLabel \[Rule] "Patch Q",  PlotLegends \[Rule] Automatic,  *)
(* Plot points with derivatives of growth rate = 0 as dots *)
d0v = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempdGrowth[[hind, vind, 1]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
d0h = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempdGrowth[[hind, vind, 2]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
d0Plot = ListPlot[Select[Flatten[Table[{h, v, d0h[h, v], d0v[h, v]}, {h, Range[0, 1, 1/10]}, {v, Range[0, 1, 1/10]}], 1], #[[3]] ==0 && #[[4]] == 0&][[All, 1;; 2]], PlotStyle -> {Black, PointSize [0.01]}];
(* Plot the derivative of the growth rate as a stream plot *)
dPlot = ListStreamPlot[Flatten[Table[{coordEquil[[hind, vind]],tempdGrowth[[hind, vind]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], StreamColorFunction -> myLineColors, StreamScale -> {Large, Automatic, 0.04}, StreamPoints -> Coarse, RegionFunction -> (((d0v[#1, #2] == 0) && (d0h[#1, #2] == 0)) == False&)];

(* Overlay the derivative of the growth rate on the plots of the equilibria in each patch *)
{Show[iPplot, dPlot, d0Plot],  Show[iQplot, dPlot, d0Plot]}], {plot, Length[plotOrder]}]]


(* ::Text:: *)
(*makeADfigure produces files with plots from the analytical model. It takes as input two strings (mainFile and suppFile) giving the files where the main analytical results and the results from any extra equilibria should be saved. It also takes as the initial horizontal and vertical transmission rates (hlist and vlist), 3 dispersal rates (d1, d2, and d3) to calculate transmission evolution at, the number of infectious contacts (c), fecundities of uninfected and infected individuals in patch P (fPU and fPI) and patch Q (fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI). It also takes 5 optional parameters : the names of the patches (patchNames, a list of 2 strings for the names to use for two patches when plotting), the precision to use when finding the ecological equilibrium (precision), the amount to shift the transmission rates to find a stable equilibrium if one is not found at a pair of initial transmission rates (epsilon), the difference in transmission rates used to calculate the derivative of the mutant growth rate (delta), and the maximum time to try to find a growth rate derivative analytically before switching to a numerical approximation (timeOut).*)


(* ::Text:: *)
(*plotBasinsAttraction plots the size of basins of attraction vs. dispersal rate. It takes as input essTable, a table of the fraction of trajectories leading to each ESS. The first column of essTable should have the names of the ESSs in the 2nd through last rows. The first row should have the dispersal rates (as numbers) in the second through last columns. Otherwise, entry[[i, j]] should give the number of trajectories leading to essTable[[i, 1]] when the dispersal rate = essTable[[1, j]]. This function also takes as optional input any option to ListLogLinearPlot.*)


plotBasinsAttraction[essTable_, opts : OptionsPattern[]] :=
 ListLogLinearPlot[Table[{essTable[[1, j]], essTable[[i, j]]/Total[essTable[[2 ;; Dimensions[essTable][[1]], j]]]}, 
   {i, 2,  Dimensions[essTable][[1]]}, {j, 2, Dimensions[essTable][[2]]}], Evaluate[FilterRules[{opts}, Options[ListLogLinearPlot]]], 
  FrameLabel -> {{"Fraction of Initial\nTransmission Rates\nThat Lead to ESS", ""}, {"Dispersal Rate", ""}}, 
  BaseStyle -> myFontStyle, RotateLabel -> False, Frame -> True, 
  ImageSize -> 400, Joined -> True, PlotLegends -> Map[Style @@ Join[ToString[#], myFontStyle] &, essTable[[2 ;; Length[essTable], 1]]]]


(* ::Section:: *)
(*Simulation Model Plots*)


(* ::Text:: *)
(*plotSimMulti plots a heatmap of the data from multiple simulations with the same parameter values. It takes as input four lists of file names: filesPH and filesMH (files of average horizontal transmission rates in the patch(es) where the symbiont is a parasite and a mutualist, respectively) and filesPV and filesMV (the same as filesPH and filesMH except for vertical transmission rates).The output is a pair of heat maps showing the distribution of transmission rates at the last time point and a statement of when the last time point is for the first simulation in filesPH.*)


plotSimMulti[filesPH_, filesMH_, filesPV_, filesMV_] :=
Module[{sims, simPlotP, simPlotM, lastTime = Import[filesPH[[1]], {"Data", -1}][[2]],
dataPH = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesPH}]],
dataMH = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesMH}]],
dataPV = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesPV}]],
dataMV = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesMV}]]},
sims = Length[dataPH]; (* Number of simulations run *)
(* Plot simulated data as a heat map showing the number of simulations with different transmission rates at the last timepoint. *)
(* Plot of transmission rates in patch(es) where symbiont is parasitic *)
simPlotP = DensityHistogram[Table[{ dataPH[[i]], dataPV[[i]]}, {i, sims}], {0,1, 1/5}, "Probability", ColorFunction -> myColorsSim, ImageSize -> 200, PlotRange -> {{0, 1}, {0, 1}}, ColorFunctionScaling -> False, LabelStyle -> myFontStyle,  FrameLabel -> {"Horiz. Transmission Rate", "Vert. Transmission Rate"}];
(* Plot of transmission rates in patch(es) where symbiont is mutualistic *)
simPlotM = DensityHistogram[Table[{ dataMH[[i]], dataMV[[i]]}, {i, sims}], {0,1, 1/5}, "Probability", ColorFunction -> myColorsSim, ImageSize -> 200, PlotRange -> {{0, 1}, {0, 1}}, ColorFunctionScaling -> False, LabelStyle -> myFontStyle,  FrameLabel -> {"Horiz. Transmission Rate", "Vert. Transmission Rate"}];
{{simPlotP, simPlotM}, StringJoin["Transmission rates at time ", ToString[lastTime, FormatType -> TraditionalForm]]}]


(* ::Text:: *)
(*plotSimInfectionMulti plots the fraction of infected hosts at the final horizontal and vertical transmission of the input simulation files. It takes as input six lists of file names:  filesPH and filesMH (files of average horizontal transmission rates in the patch(es) where the symbiont is a parasite and a mutualist, respectively), filesPV and filesMV (the same as filesPH and filesMH except for vertical transmission rates), and filesPInf and filesMInf (files of the fraction of infected hosts in the patch(es) where the symbiont is a parasite and a mutualist, respectively). The output is a pair of scatter plots, one for each patch, where the coordinates indicate the transmission rates  and the colors of the points indicate the fraction of infected hosts.*)


plotSimInfectionMulti[filesPH_, filesMH_, filesPV_, filesMV_, filesPInf_, filesMInf_] :=
Module[{sims, simPlotP, simPlotM, outlinePlotP, outlinePlotM, lastTime = Import[filesPH[[1]], {"Data", -1}][[2]],
dataPH = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesPH}]],
dataMH = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesMH}]],
dataPV = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesPV}]],
dataMV = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesMV}]],
dataPInf = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesPInf}]],
dataMInf = Catenate[Table[Import[file, {"Data", -1}][[3;;-1]], {file, filesMInf}]]},
sims = Length[dataPH]; (* Number of simulations run *)

(* Plot of fraction of infected hosts (color) and transmission rates (coordinates) in patch(es) where symbiont is parasitic *)
simPlotP = ListPlot[Table[Style[{dataPH[[i]], dataPV[[i]]}, myColors[dataPInf[[i]]]], {i, sims}], PlotStyle -> PointSize[0.05], AspectRatio -> 1, PlotRange -> {{0, 1}, {0, 1}}, ImageSize -> 200, LabelStyle -> myFontStyle,  Frame -> True, FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}];
outlinePlotP = ListPlot[Select[Table[{dataPH[[i]], dataPV[[i]], dataPInf[[i]]}, {i, sims}], (#[[3]] == 0)&][[All, 1;;2]], PlotMarkers -> {Graphics[{Thin, Black, Circle[]}], 0.05}]; (* outlines the points with no infection in black *)
(* Plot of fraction of infected hosts (color) and transmission rates (coordinates) in patch(es) where symbiont is mutualistic *)
simPlotM = ListPlot[Table[Style[{dataMH[[i]], dataMV[[i]]}, myColors[dataMInf[[i]]]], {i, sims}], PlotStyle -> PointSize[0.05], AspectRatio -> 1, PlotRange -> {{0, 1}, {0, 1}}, ImageSize -> 200, LabelStyle -> myFontStyle,  Frame -> True, FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}];
outlinePlotM = ListPlot[Select[Table[{dataMH[[i]], dataMV[[i]], dataMInf[[i]]}, {i, sims}], (#[[3]] == 0)&][[All, 1;;2]], PlotMarkers -> {Graphics[{Thin, Black, Circle[]}], 0.05}]; (* outlines the points with no infection in black *)
{{Show[simPlotP, outlinePlotP], Show[simPlotM, outlinePlotM]}, StringJoin["Transmission rates at time ", ToString[lastTime, FormatType -> TraditionalForm]]}]


(* ::Section:: *)
(*Summarizing Data*)


(* ::Text:: *)
(*myTally tallies the contents of a list using the input function tallyBy to determine if two entries are equal. It performs the input function nameBy on the elements corresponding to each tally, to allow useful descriptions of the items tallied.*)


myTally[list_, tallyBy_, nameBy_] := 
 Module[{tempList}, 
  tempList = Tally[list, tallyBy];
  Flatten[Join[{Map[nameBy, tempList[[All, 1]]]}, {tempList[[All, 2]]}], {2}]]


(* ::Text:: *)
(*makeESStable makes a table ready to be filled with counts of trajectories leading to evolutionary stable strategies. It takes as input a list of dispersal rates where transmission evolution was investigated.*)


makeESStable[dlist_] := 
{Prepend[dlist, "ESS"]}


(* ::Text:: *)
(*addESS adds a count of trajectories leading an ESSs to an existing table made by makeESStable. It takes as input a table of counts of evolutionary stable strategies at different dispersal rates (essTable), a dispersal rate (d, must be present in essTable), the ESS (ess), and the number of trajectories leading to the ESS (count). If the ESS already has trajectories at that dispersal rate, count will be added to the existing number of trajectories.*)


addESS::tableError = "d not present in essTable.";
addESS[essTable_, d_, ess_, count_] :=
Module[{essTableTemp = essTable},
If[!MemberQ[essTable[[1]], d], Message[addESS::tableError]; Return[essTable]];
If[!MemberQ[essTable[[All, 1]], ess], essTableTemp = Append[essTable, Join[{ess}, Table[0, {i, Dimensions[essTable][[2]] - 1}]]]];
essTableTemp[[Position[essTableTemp, ess][[1, 1]], Position[essTableTemp, d][[1, 2]]]] += count;
Return[essTableTemp]]


(* ::Text:: *)
(*talliesToESS takes a list of tallies of evolutionary endpoints at different dispersal rates and produces a table of the number of trajectories leading to each endpoint at each dispersal rate. tallies should be a list with entries of the form {d, tal} where d is the dispersal rate and tal is a list of tallies (produced by myTallies) of trajectories. tal may be of length > 1 when there are multiple possible stable ecological equilibria at pair of transmission rates.*)


talliesToESS[tallies_] :=Module[{nTallies, minTallies, maxTallies, tableOrder, essTables},
nTallies = Map[Length, tallies[[All,2]]]; (* Number of sets of trajectories for each dispersal rate *)
{minTallies, maxTallies} = {Min[nTallies], Max[nTallies]}; (* Minimum (minTallies) and maximum (maxTallies) sets of trajectories any point has *)
tableOrder = Tuples[Table[Range[n], {n, minTallies, maxTallies}]]; (* Order in which to create essTables based on the tallied trajectories. First level (rows) of the list corresponds to the essTable. Second level (columns) corresponds to the number of sets of trajectories a dispersal rate has. tableOrder[[row, column]] gives the tally of trajectories to put in the essTable for the rowth essTable for all points with column # of lists of trajectories. *)
essTables = Table[makeESStable[tallies[[All,1]]], {t, Length[tableOrder]}];
Table[Module[{tempTallies = tallies[[d, 2, tableOrder[[t, Length[tallies[[d, 2]]] - minTallies + 1]]]]}, Table[essTables[[t]] = addESS[essTables[[t]],tallies[[d, 1]],  tempTallies[[i, 1]], tempTallies[[i, 2]]];, {i, Length[tempTallies]}];], {t, Length[tableOrder]}, {d, Length[tallies]} ];
Return[essTables]]



End[]
EndPackage[]
