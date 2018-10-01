(* ::Package:: *)

(* ::Text:: *)
(*This file is `plottingFunctions.wl`.*)


(* ::Text:: *)
(*This script contains functions for plotting and summarizing data.*)


BeginPackage["plottingFunctions`"]


myFontStyle;
myEquilLegend;
mySimLegend;
plotAD::usage = "plotAD[equil, coordEquil, dGrowth] plots the equilibrium fraction of infected hosts in each patch, overlaid with a stream plot showing the direction of transmission mode evolution. It takes as input matrices of equilibria (equil, entries in the format {iP ->  iPvalue, iQ -> iQvalue}), coordinates of equilibria (coordEquil, entries in the format {hres -> hresvalue, vres -> vresvalue}), and mutant growth rates (dGrowth, entries in the format {ddh -> ddhvalue, ddv -> ddvvalue}). Because a single point may have multiple equilibria, this function produces plots for each equilibrium, assuming that points with the same number of equilibria have them ordered so that point1's equilibrium 1 should be plotted with point2's equilibrium 1, and point1's equilibrium 2 should be plotted with point2's equilibrium 2.";
plotBasinsAttraction::usage = "plotBasinsAttraction[essTable] plots the size of basins of attraction vs. dispersal rate. It takes as input essTable, a table of the fraction of trajectories leading to each ESS. The first column of essTable should have the names of the ESSs in the 2nd through last rows. The first row should have the dispersal rates (as numbers) in the second through last columns. Otherwise, entry[[i, j]] should give the number of trajectories leading to essTable[[i, 1]] when the dispersal rate = essTable[[1, j]]. This function also takes as optional input any option to ListLogLinearPlot.";
plotR0::Usage = "plotR0[equil, coordEquil, r0, colors, fontStyle, colorFunctionScaling] plots the basic reproductive number. It takes as input the equilibria (equil), their coordinates (coordEquil) and the R_0 values for each equilibrium (R_0). It will also take as optional input a colorfunction and fontstyle to use, as well as whether to scale the colorfunction (default True).";
plotSimTransmission::usage = "plotSimTransmissionAndInfection[file, patchTypes] plots a heatmap of the average transmission rates at the end of several simulations. Takes as input file, a csv file listing the transmission probabilities and fraction of infected hosts for several simulations (columns) at different time points (rows). Also takes as input patchTypes, the number of types of patches in thesummary file.\nFor the input file, the data for each time point should be listed on several rows, first by patch (name of patch in 2nd column) and then by horizontal transmission probability, fraction infected, and vertical transmission probability. The simulated datashould begin in the third column.\nExample:\ntime\t patch\t info\t sim1\t sim2\t ...\n1\t M\t Avg HT\t 0.5  \t 0.53 \n1\t M\t Infected\t 0.2  \t 0.21 \n1\t M\t Avg VT\t 0.8  \t 0.9  \n1\t P\t Avg HT\t 0.44 \t 0.4  \n...........";
plotSimTransmissionAndInfection::usage = "plotSimTransmissionAndInfection[file, patchTypes] plots a point for each simulation in file. Points show the ending transmission probabilities and are colored to indicate the fraction of infected hosts. Takes as input file, a csv file listing the transmission probabilities and fraction of infected hosts for several simulations (columns) at different time points (rows). Also takes as input patchTypes, the number of types of patches in thesummary file.\nFor the input file, the data for each time point should be listed on several rows, first by patch (name of patch in 2nd column) and then by horizontal transmission probability, fraction infected, and vertical transmission probability. The simulated datashould begin in the third column. Example\ntime\t patch\t info\t sim1\t sim2...\n1\t M\t Avg HT\t 0.5\t 0.53 \n1\t M\t Infected\t 0.2\t 0.21 \n1\t M\t Avg VT\t 0.8\t 0.9 \n1\t P\t Avg HT\t 0.44 \t 0.4 \n...";
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


myColors::usage = "A color function for the fraction of infected hosts. Takes as input the fraction of infected hosts and returns a color.";
myColors[f_] := ColorData["Candy"][(4/20) + (15f/20)];


noInfectionColor::usage = "Color that indicates that infection cannot be sustained in either patch. Used in plotAD function.";
noInfectionColor = White;


(* ::Text:: *)
(*Font style*)


myFontStyle = {FontFamily -> "TimesNewRoman", FontColor -> Black, FontSize -> 14};


(* ::Text:: *)
(*Legend for fraction of infected hosts (based on myColors)*)


myEquilLegend = BarLegend[{myColors[#]&,{0,1}}, LegendLayout -> "Row",LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 200];


(* ::Text:: *)
(*Legend for fraction of simulations (based on myColorsSim)*)


mySimLegend = BarLegend[{myColors[#]&,{0,1}}, ScalingFunctions -> {#^(1/4)&, #^(4)&}, Ticks -> {0, 0.001, 0.02, 0.1, 0.5, 1.0}, 
LegendLayout -> "Row", LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 200];


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
Table[Module[{tempEquil, tempdGrowth, iPplot, iQplot, dPlot, d0h, d0v, d0Plot, iP0, iQ0, iP0iQ0Plot},
(* Select the equilibria (tempEquil) and growth rates (tempdGrowth) to plot *)
tempEquil = Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, equil, {2}];
tempdGrowth = Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, dGrowth, {2}];

(* Plot the equilibrium fraction of infected hosts in patches P (iPplot) and Q (iQplot) *)
iPplot = ListDensityPlot[Flatten[Table[{(coordEquil[[hind, vind, 1]]), (coordEquil[[hind, vind, 2]]), (tempEquil[[hind, vind, 1]])}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}], 1],PlotRange -> {{0, 1},{0,1}, {0, 1}},  LabelStyle -> myFontStyle, ColorFunction -> myColors, ColorFunctionScaling -> False, FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}, RotateLabel -> True];
iQplot = ListDensityPlot[Flatten[Table[{(coordEquil[[hind, vind, 1]]), (coordEquil[[hind, vind, 2]]), (tempEquil[[hind, vind, 2]])}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}], 1], PlotRange -> {{0, 1},{0,1}, {0, 1}},  LabelStyle -> myFontStyle, ColorFunction -> myColors, ColorFunctionScaling -> False,  FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}, RotateLabel -> True];

(* Plot points with derivatives of growth rate = 0 as dots *)
d0h = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempdGrowth[[hind, vind, 1]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
d0v = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempdGrowth[[hind, vind, 2]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
d0Plot = ListPlot[Select[Flatten[Table[{h, v, d0h[h, v], d0v[h, v]}, {h, Range[0, 1, 1/10]}, {v, Range[0, 1, 1/10]}], 1], #[[3]] ==0 && #[[4]] == 0&][[All, 1;; 2]], PlotStyle -> {Black, PointSize [0.01]}];
(* Plot the derivative of the growth rate as a stream plot *)
dPlot = ListStreamPlot[Flatten[Table[{coordEquil[[hind, vind]],tempdGrowth[[hind, vind]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], StreamColorFunction -> myLineColors, StreamScale -> {Large, Automatic, 0.04}, StreamPoints -> Coarse, RegionFunction -> (((d0v[#1, #2] == 0) && (d0h[#1, #2] == 0)) == False&)];

(* Find the regions where there is no infection in both patches *)
iP0 = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempEquil[[hind, vind, 1]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
iQ0 = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempEquil[[hind, vind, 2]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
iP0iQ0Plot = RegionPlot[iP0[h, v] <= 0 && iQ0[h, v] <= 0, {h, 0, 1}, {v, 0, 1}, PlotStyle -> noInfectionColor, BoundaryStyle -> None];

(* Overlay the derivative of the growth rate on the plots of the equilibria in each patch *)
{Show[iPplot, iP0iQ0Plot, dPlot, d0Plot],  Show[iQplot, iP0iQ0Plot, dPlot, d0Plot]}], {plot, Length[plotOrder]}]]


(* ::Text:: *)
(*plotBasinsAttraction plots the size of basins of attraction vs. dispersal rate. It takes as input essTable, a table of the fraction of trajectories leading to each ESS. The first column of essTable should have the names of the ESSs in the 2nd through last rows. The first row should have the dispersal rates (as numbers) in the second through last columns. Otherwise, entry[[i, j]] should give the number of trajectories leading to essTable[[i, 1]] when the dispersal rate = essTable[[1, j]]. This function also takes as optional input any option to ListLogLinearPlot.*)


plotBasinsAttraction[essTable_, opts : OptionsPattern[]] :=
 ListLogLinearPlot[Table[{essTable[[1, j]], essTable[[i, j]]/Total[essTable[[2 ;; Dimensions[essTable][[1]], j]]]}, 
   {i, 2,  Dimensions[essTable][[1]]}, {j, 2, Dimensions[essTable][[2]]}], Evaluate[FilterRules[{opts}, Options[ListLogLinearPlot]]], 
  FrameLabel -> {{"Fraction of Initial\nTransmission Rates\nThat Lead to ESS", ""}, {"Dispersal Rate", ""}}, 
  BaseStyle -> myFontStyle, RotateLabel -> False, Frame -> True, 
  ImageSize -> 400, Joined -> True, PlotLegends -> Map[Style @@ Join[ToString[#], myFontStyle] &, essTable[[2 ;; Length[essTable], 1]]]]


(* ::Text:: *)
(*plotR0 plots the basic reproductive number. It takes as input the equilibria (equil), their coordinates (coordEquil) and the R_0 values for each equilibrium (R_0).*)


plotR0[equil_, coordEquil_, r0_, colors_:myColors, fontStyle_:myFontStyle, colorFunctionScaling_:True] := 
Module[{nEquil, minEquil, maxEquil, plotOrder, iP0, iQ0, iP0iQ0Plot},
(* Determine the number of plots needed to represent all equilibria and produce an ordering of those plots *)
nEquil = Map[Length, equil, {2}]; (* Number of equilibria at each horizontal transmission rate/verticaltransmission rate pair *)
{minEquil, maxEquil} = {Min[nEquil], Max[nEquil]}; (* Minimum (minEquil) and maximum (maxEquil) equilibria any point has *)
plotOrder = Tuples[Table[Range[n], {n, minEquil, maxEquil}]]; (* Order in which to plot equilibria. First level (rows) of the list corresponds to the plot. Second level (columns) corresponds to the number of equilibria a point has. plotOrder[[row, column]] gives the equilibrium to plot for the rowth plot for all points with column # of equilibria. *)

(* Make plots for each combination of equilibria *)
Table[Module[{tempEquil,tempR0, r0plot},
(* Select the equilibria (tempEquil) to plot *)
tempEquil = Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, equil, {2}];
tempR0= Map[#[[plotOrder[[plot, Length[#] - minEquil + 1]]]]&, r0, {2}];

(* Plot r0 *)
r0plot = ListDensityPlot[Flatten[Table[{(coordEquil[[hind, vind, 1]]), (coordEquil[[hind, vind, 2]]), (tempR0[[hind, vind]])}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}], 1],PlotRange -> {{0, 1},{0,1}, {0, Max[tempR0]}},  LabelStyle -> fontStyle, ColorFunction -> colors, ColorFunctionScaling -> colorFunctionScaling, Frame -> True, FrameTicks -> {{{0, 0.5, 1}, None}, {{0, 0.5, 1}, None}}];

(* Find the regions where there is no infection in both patches *)
iP0 = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempEquil[[hind, vind, 1]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
iQ0 = Interpolation[Flatten[Table[{coordEquil[[hind, vind]], tempEquil[[hind, vind, 2]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1], InterpolationOrder -> 1];
iP0iQ0Plot = RegionPlot[iP0[h, v] <= 0 && iQ0[h, v] <= 0, {h, 0, 1}, {v, 0, 1}, PlotStyle -> White, BoundaryStyle -> None, Frame -> True, FrameTicks -> {{{0, 0.5, 1}, None}, {{0, 0.5, 1}, None}}]; 

(* Overlay the derivative of the growth rate on the plots of the equilibria in each patch *)
Show[r0plot, iP0iQ0Plot]], {plot, Length[plotOrder]}]]


(* ::Section::Closed:: *)
(*Simulation Model Plots*)


(* ::Text:: *)
(*plotSimTransmission produces a heatmap indicating the relative frequencies of transmission rates at the endpoints of simulations for the simulations and patches given in a summary file. *)


plotSimTransmission[file_, patchTypes_] :=
Module[{results = Import[file, {"Data", patchTypes * -3;;All}]},
Table[{results[[3*patch, 2]],
DensityHistogram[Transpose[results[[{3*patch-2, 3*patch}, 4;;All]]], {0, 1, 1/7}, "Probability", ColorFunction -> myColors, ColorFunctionScaling -> False, AspectRatio -> 1, 
ScalingFunctions -> {#^(1/4)&, #^(4)&},
PlotRange -> {{0, 1}, {0, 1}}, ImageSize -> 200, LabelStyle -> myFontStyle,  Frame -> True, FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}]}, {patch, patchTypes}]]


(* ::Text:: *)
(*plotSimTransmissionAndInfection produces a scatter plot indicating the fraction of infected hosts (color) and transmission probabilities (coordinates) for each simulation and patch given in a summary file.*)


plotSimTransmissionAndInfection[file_, patchTypes_] :=
Module[{results = Import[file, {"Data", patchTypes * -3;;All}]},
Table[{results[[3*patch, 2]],
Show[ListPlot[Map[Style[{#[[1]], #[[3]]}, myColors[#[[2]]]]&, Transpose[results[[3*patch-2;;3*patch, 4;;All]]]],
PlotStyle -> PointSize[0.025], AspectRatio -> 1, PlotRange -> {{0, 1}, {0, 1}}, ImageSize -> 200, LabelStyle -> myFontStyle,  Frame -> True, FrameLabel -> {"Horiz. Trans. Prob.", "Vert. Trans. Prob."}],
ListPlot[Transpose[results[[{3*patch-2,3*patch}, 4;;All]]], PlotMarkers -> {Graphics[{Thin, Black, Opacity[0.15], Circle[]}], 0.025}]]}, {patch, patchTypes}]]


(* ::Section::Closed:: *)
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
