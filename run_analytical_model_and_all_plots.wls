#!/usr/bin/env wolframscript
(* ::Package:: *)

If[!Quiet[FailureQ[NotebookDirectory[]]], SetDirectory[NotebookDirectory[]]];
If[!DirectoryQ["Figures"], CreateDirectory["Figures"]];
If[!DirectoryQ["wdx"], CreateDirectory["wdx"]];


Get["hostControlAD`"]
Get["symbiontControlAD`"]
Get["plottingFunctions`"]
LaunchKernels[3]
ParallelNeeds["hostControlAD`"]
ParallelNeeds["symbiontControlAD`"]
ParallelNeeds["plottingFunctions`"]


(* ::Section:: *)
(*Plots of simulation results*)


(* ::Subsection:: *)
(*Host control, population size 200, transmission rates (main text)*)


hostControlFilesFec = Table["Simulation results/Host_Control_rep_" <> letter <> "/" <> i <> 
", f = 0.5, s = 1.0, m = 1.0, N = 200, patches = 2, "<> j <> ".csv", 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, 
{i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}, 
{letter, {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}}];
hostControlPlotsFec = ParallelMap[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, hostControlFilesFec];


hostControlFilesEst = Table["Simulation results/Host_Control_rep_" <> date <> "/" <> i <> 
", f = 1.0, s = 0.5, m = 1.0, N = 200, patches = 2, "<> j <> ".csv", 
 {j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, 
 {i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}, 
{letter, {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}}];
 hostControlPlotsEst = ParallelMap[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, hostControlFilesEst];


(* ::Text:: *)
(*Distribution of average transmission probabilities at simulation endpoints (Main text figure 4)*)


Export["Figures/HostControl_MainText.pdf", 
Join[Transpose[hostControlPlotsEst[[All, 1]]][[{2, 1}]], Transpose[hostControlPlotsFec[[All, 1]]][[{2, 1}]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Host Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Host Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of simulations ending\nat given average transmission rates", mySimLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Subsection:: *)
(*Host control, population size 200, fraction of infected hosts (supplement)*)


hostControlFilesFecInfection = Table["Simulation results/Host_Control_rep_" <> letter <> "/" <> i <> 
", f = 0.5, s = 1.0, m = 1.0, N = 200, patches = 2, "<> j <> ".csv", 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}},
{i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M", "Fraction infected hosts patch P", "Fraction infected hosts patch M"}},  
{letter, {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}}];
 hostControlPlotsFecInfection = Map[plotSimInfectionMulti[#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], #[[6]]]&, hostControlFilesFecInfection];


hostControlFilesEstInfection = Table["Simulation results/Host_Control_rep_" <> letter <> "/" <> i <> 
", f = 1.0, s = 0.5, m = 1.0, N = 200, patches = 2, "<> j <> ".csv", 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}},
{i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M", "Fraction infected hosts patch P", "Fraction infected hosts patch M"}},  
{letter, {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"}}];
 hostControlPlotsEstInfection = ParallelMap[plotSimInfectionMulti[#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], #[[6]]]&, hostControlFilesEstInfection];


hostControlPlotsFecInfectionRaster = Map[Rasterize[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
	ImageResolution -> 300]&, hostControlPlotsFecInfection, {3}];
hostControlPlotsEstInfectionRaster = Map[Rasterize[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
	ImageResolution -> 300]&, hostControlPlotsEstInfection, {3}];


(* ::Text:: *)
(*Fraction of infected hosts at different simulation endpoints (Supplement figure 8)*)


Export["Figures/HostControl_FractionInfected.pdf", 
Join[Transpose[hostControlPlotsEstInfectionRaster[[All, 1]]][[{2, 1}]], Transpose[hostControlPlotsFecInfectionRaster[[All, 1]]][[{2, 1}]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Host Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Host Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Subsection:: *)
(*Host control, population size 1000 (supplement)*)


largePopFilesFec = Table["Simulation results/Host_Control_rep_" <> letter <> "/" <> i <> 
", f = 0.5, s = 1.0, m = 1.0, N = 1000, patches = 2, "<> j <> ".csv", 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, 
{i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}, 
{letter, {"a", "b", "c", "d", "e"}}];
largePopPlotsFec = Map[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, largePopFilesFec];


largePopFilesEst = Table["Simulation results/Host_Control_rep_" <> letter <> "/" <> i <> 
", f = 1.0, s = 0.5, m = 1.0, N = 1000, patches = 2, "<> j <> ".csv", 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, 
{i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}, 
{letter, {"a", "b", "c", "d", "e"}}];
largePopPlotsEst = Map[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, largePopFilesEst];


(* ::Text:: *)
(*Simulation of host evolution in larger population (Supplement figure 9)*)


Export["Figures/HostControl_LargePop.pdf", 
Join[Transpose[largePopPlotsEst[[All, 1]]][[{2, 1}]], Transpose[largePopPlotsFec[[All, 1]]][[{2, 1}]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Host Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Host Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of simulations ending\nat given average transmission rates", mySimLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Subsection:: *)
(*Symbiont control, population size 200 (supplement)*)


symbControlFilesFec = Table[{"Simulation results/Symbiont_Control/" <> i <> ", f = 0.5, s = 1.0, m = 1.0, N = 200, patches = 2, "<> j <> ".csv"}, 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, {i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}];
symbControlPlotsFec = Map[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, symbControlFilesFec];


symbControlFilesEst = Table[{"Simulation results/Symbiont_Control/" <> i <> ", f = 1.0, s = 0.5, m = 1.0, N = 200, patches = 2, "<> j <> ".csv"}, 
{j, {"d = 0.005", "d = 0.05", "d = 0.5"}}, {i, {"Horiz trans evo patch P", "Horiz trans evo patch M", "Vert trans evo patch P", "Vert trans evo patch M"}}];
symbControlPlotsEst = Map[plotSimMulti[#[[1]], #[[2]], #[[3]], #[[4]]]&, symbControlFilesEst];


(* ::Text:: *)
(*Simulation of symbiont evolution (Supplement figure 11)*)


Export["Figures/SymbiontControlSimulations.pdf", 
Join[Transpose[symbControlPlotsEst[[All, 1]]][[{2, 1}]], Transpose[symbControlPlotsFec[[All, 1]]][[{2, 1}]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Host Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Host Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob.", "       Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of simulations ending\nat given average transmission rates", mySimLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Section::Closed:: *)
(*Symbiont affects one fitness component symmetrically*)


(* ::Subsection::Closed:: *)
(*Calculations*)


(* ::Text:: *)
(*Symbiont affects host fecundity*)


resultsFecAlone = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecAlone.wdx", resultsFecAlone];


(* ::Text:: *)
(*Symbiont affects host establishment probability (lifespan)*)


resultsEstAlone = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1, fPI -> 1, fQU -> 1, fQI -> 1, sPU -> 1/2, sPI -> 1, sQU -> 1, sQI -> 1/2, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsEstAlone.wdx", resultsEstAlone];


(* ::Text:: *)
(*Symbiont affects adult host mortality (lifespan)*)


resultsMortAlone = ParallelTable[Module[{eq, coord, gHost, gSymb, 
(* Mortality effects are 2 instead of 1/2 because a lower mortality is good, while lower fecundity/establishment is bad *)
params = {c  -> 1, fPU -> 1, fPI -> 1, fQU -> 1, fQI -> 1, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 2, mPI -> 1, mQU -> 1, mQI -> 2}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsMortAlone.wdx", resultsMortAlone];


(* ::Subsection::Closed:: *)
(*Figures*)


resultsFecAlone = Import["wdx/resultsFecAlone.wdx"];
resultsEstAlone = Import["wdx/resultsEstAlone.wdx"];
resultsMortAlone = Import["wdx/resultsMortAlone.wdx"];


plotsFecAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAlone];
plotsFecAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsFecAlone];
plotsEstAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstAlone];
plotsEstAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsEstAlone];
plotsMortAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsMortAlone];
plotsMoreAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsMortAlone];


(* ::Text:: *)
(*Host evolution when symbiont affects host establishment probability or fecundity (Main text figure 3)*)


Export["Figures/HostControl_MainScenario.pdf",
Join[Transpose[First /@ plotsEstAloneHost], Transpose[First /@ plotsFecAloneHost]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Text:: *)
(*Host evolution at alternate ecological equilibria when symbiont affects establishment probability (Supplement figure 1)*)


Export["Figures/HostControlEstAlone_extraEquil.pdf", Transpose[Map[Last, plotsEstAloneHost[[1;;2]]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]]}}, 
Join[{{Item[Rotate["Symbiont Affects Lifespan", 90Degree], Background -> Lighter[Gray]],
Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}}, 
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts\nthat are infected", myEquilLegend}}], SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Text:: *)
(*Host evolution when symbiont affects host mortality (Supplement figure 2)*)


Export["Figures/HostControlMortAlone.pdf", Join[Transpose[plotsMortAloneHost[[All, 1]]], Transpose[plotsMortAloneHost[[1;;2, 2]]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}}, 
Join[{{Item[Rotate["Ecological Equilibria - Set 1", 90Degree], Background -> Lighter[Gray]],
Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]},
{Item[Rotate["Ecological Equilibria - Set 2", 90Degree], Background -> LightGray],
Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}}, 
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Text:: *)
(*Symbiont evolution when symbiont affects host establishment probability or fecundity (Supplement figure 10)*)


Export["Figures/SymbiontControl_MainScenario.pdf",
Join[Transpose[First /@ plotsEstAloneSymb], Transpose[First /@ plotsFecAloneSymb]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Lifespan", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Fecundity", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Section::Closed:: *)
(*Host-symbiont conflict*)


(* ::Subsection::Closed:: *)
(*Calculations*)


resultsFecAlone = Import["wdx/resultsFecAlone.wdx"];
resultsEstAlone = Import["wdx/resultsEstAlone.wdx"];


conflictFec = Map[MapThread[{#1[[1]], #1[[2]], If[(#2[[1]] == 0) && (#2[[2]] == 0), ,
If[(Norm[#3] == 0) || (Norm[#4] == 0), 0, (#3 . #4)/(Norm[#3] * Norm[#4]), -0.99]]}&,
{Flatten[(coordEquil/. #), 1], Chop[Flatten[Map[First, (equil/. #), {2}], 1]], Chop[Flatten[Map[First, (dGrowthHost/. #), {2}], 1]], 
Chop[Flatten[Map[First, (dGrowthSymb/. #), {2}], 1]]}]&, resultsFecAlone];


conflictEst = Map[MapThread[{#1[[1]], #1[[2]], If[(#2[[1]] == 0) && (#2[[2]] == 0), ,
If[(Norm[#3] == 0) || (Norm[#4] == 0), 0, (#3 . #4)/(Norm[#3] * Norm[#4]), -0.99]]}&,
{Flatten[(coordEquil/. #), 1], Chop[Flatten[Map[First, (equil/. #), {2}], 1]], Chop[Flatten[Map[First, (dGrowthHost/. #), {2}], 1]], 
Chop[Flatten[Map[First, (dGrowthSymb/. #), {2}], 1]]}]&, resultsEstAlone];


(* ::Subsection::Closed:: *)
(*Figure*)


myConflictColors = ColorData["Candy"][(11/20) + (9#/20)]&;
myConflictLegend = BarLegend[{myConflictColors, {-1, 1}}, {-1, -0.6, -0.3, 0, 0.3, 0.6, 1}, 
LegendLayout -> "Row",LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 300];


conflictPlots = Map[ListContourPlot[#, BaseStyle -> myFontStyle, ImageSize -> 150, PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ClippingStyle -> White,
ColorFunctionScaling -> False, ColorFunction -> myConflictColors, Contours -> {-0.6, -0.3, 0, 0.3, 0.6}, ContourStyle -> {Black, Black, {Dashed, Thick}, Black, Black}]&, 
{conflictEst, conflictFec}, {2}];


conflictPlotsRaster = Map[Rasterize[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}},
	FrameStyle -> ReleaseHold[Hold[Directive[FontColor, FontSize, FontFamily]]/. myFontStyle]], ImageResolution -> 300]&, conflictPlots, {2}];


(* ::Text:: *)
(*Host-symbiont conflict (Main text figure 5)*)


Export["Figures/Host_Symbiont_Conflict.pdf",
conflictPlotsRaster//Grid[Join[{{"", ""},{Item[Rotate["Symbiont Affects Lifespan", 90 Degree], Background -> Lighter[Gray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Fecundity", 90 Degree], Background -> LightGray], Rotate["Vert. Trans. Prob.", 90 Degree]}, {"", ""}},
Join[{{Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}}, 
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130], Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 
{{"         Horiz. Trans. Prob.", "         Horiz. Trans. Prob.", "         Horiz. Trans. Prob."},
{Grid[{{"Host-symbiont conflict (negative values)\nand concordance (positive values)", myConflictLegend}}], SpanFromLeft, SpanFromLeft}}], 2], 
BaseStyle -> myFontStyle, Alignment -> {Center, Center}]&, ImageResolution -> 300];


(* ::Section::Closed:: *)
(*Symbiont affects fecundity asymmetrically*)


(* ::Subsection::Closed:: *)
(*Calculations*)


resultsFecAloneMoreCosts = ParallelTable[Module[{eq, coord, gHost, gSymb, params = {c  -> 1, fPU -> 3/4, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100],disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecAloneMoreCosts.wdx", resultsFecAloneMoreCosts];


resultsFecAloneMoreBen = ParallelTable[Module[{eq, coord, gHost, gSymb, params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 3/4, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100],disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecAloneMoreBen.wdx", resultsFecAloneMoreBen];


(* ::Subsection::Closed:: *)
(*Figure*)


resultsFecAloneMoreCosts = Import["wdx/resultsFecAloneMoreCosts.wdx"];
resultsFecAloneMoreBen = Import["wdx/resultsFecAloneMoreBen.wdx"];


plotsFecAloneMoreCosts = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAloneMoreCosts];
plotsFecAloneMoreBen = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAloneMoreBen];


(* ::Text:: *)
(*Symbiont affects fecundity asymmetrically (supplement figure 3)*)


Export["Figures/HostControl_FecAlone_MoreCostsOrBen.pdf",
Join[Transpose[First /@ plotsFecAloneMoreCosts], Transpose[First /@ plotsFecAloneMoreBen]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont costs in Patch P\noutweight benefits in Patch M", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont benefits in Patch M\noutweight costs in Patch P", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Section::Closed:: *)
(*Symbiont has weak conditional effect on a second fitness component*)


(* ::Subsection::Closed:: *)
(*Calculations*)


(* ::Text:: *)
(*Symbiont affects host fecundity and weakly affects host lifespan (establishment probability) in the same direction*)


resultsFecWeakEstSameDirection = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 99/100, sPI -> 1, sQU -> 1, sQI -> 99/100, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecWeakEstSameDirection.wdx", resultsFecWeakEstSameDirection];


(* ::Text:: *)
(*Symbiont affects host fecundity and weakly affects host lifespan (establishment probability) in the opposite direction*)


resultsFecWeakEstOppositeDirection = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 1, sPI -> 99/100, sQU -> 99/100, sQI -> 1, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecWeakEstOppositeDirection.wdx", resultsFecWeakEstOppositeDirection];


(* ::Text:: *)
(*Symbiont affects host lifespan (establishment probability) and weakly affects host fecundity in the same direction*)


resultsEstWeakFecSameDirection = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 99/100, fPI -> 1, fQU -> 1, fQI -> 99/100, sPU -> 1/2, sPI -> 1, sQU -> 1, sQI -> 1/2, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsEstWeakFecSameDirection.wdx", resultsEstWeakFecSameDirection];


(* ::Text:: *)
(*Symbiont affects host lifespan (establishment probability) and weakly affects host lifespan in the opposite direction*)


resultsEstWeakFecOppositeDirection = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1, fPI -> 99/100, fQU -> 99/100, fQI -> 1, sPU -> 1/2, sPI -> 1, sQU -> 1, sQI -> 1/2, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsEstWeakFecOppositeDirection.wdx", resultsEstWeakFecOppositeDirection];


(* ::Subsection::Closed:: *)
(*Figures*)


resultsFecWeakEstSameDirection = Import["wdx/resultsFecWeakEstSameDirection.wdx"];
resultsFecWeakEstOppositeDirection = Import["wdx/resultsFecWeakEstOppositeDirection.wdx"];
resultsEstWeakFecSameDirection = Import["wdx/resultsEstWeakFecSameDirection.wdx"];
resultsEstWeakFecOppositeDirection = Import["wdx/resultsEstWeakFecOppositeDirection.wdx"];


(equil/.resultsFecWeakEstOppositeDirection[[3]])[[101, 2]] == (equil/.resultsFecWeakEstSameDirection[[3]])[[101, 2]]


plotsFecWeakEstSameDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecWeakEstSameDirection];
plotsFecWeakEstOppositeDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecWeakEstOppositeDirection];
plotsEstWeakFecSameDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstWeakFecSameDirection];
plotsEstWeakFecOppositeDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstWeakFecOppositeDirection];


(* ::Text:: *)
(*Conditional effects act in same direction (supplement figure 4)*)


Export["Figures/HostControl_WeakConditional_SameDirection.pdf",
Join[Transpose[First /@ plotsEstWeakFecSameDirection], Transpose[First /@ plotsFecWeakEstSameDirection]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Lifespan\nMore Strongly than Fecundity", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Fecundity\nMore Strongly than Lifespan", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Text:: *)
(*Conditional effects act in opposite directions (supplement figure 5)*)


Export["Figures/HostControl_WeakConditional_OppositeDirection.pdf",
Join[Transpose[First /@ plotsEstWeakFecOppositeDirection], Transpose[First /@ plotsFecWeakEstOppositeDirection]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects Lifespan\nMore Strongly than Fecundity", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Affects Fecundity\nMore Strongly than Lifespan", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Section::Closed:: *)
(*Symbiont affects fecundity conditionally, has unconditional lifespan cost*)


(* ::Subsection::Closed:: *)
(*Calculations*)


(* ::Text:: *)
(*Symbiont affects host fecundity conditionally and has a small unconditional lifespan (establishment) cost*)


resultsFecConditEstCost = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 1, sPI -> 7/8, sQU -> 1, sQI -> 7/8, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsFecConditEstCost.wdx", resultsFecConditEstCost];


(* ::Subsection::Closed:: *)
(*Figures*)


resultsFecConditEstCost = Import["wdx/resultsFecConditEstCost.wdx"];


plotsFecConditEstCost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecConditEstCost];


(* ::Text:: *)
(*Symbiont affects fecundity conditionally, has unconditional lifespan cost (Supplement figure 6)*)


Export["Figures/HostControl_FecCondit_EstCost.pdf", Transpose[Map[First, plotsFecConditEstCost]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{ {Item[Rotate["Symbiont Conditionally Affects Fecundity\nUnconditionally Decreases Lifespan", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Section::Closed:: *)
(*Symbiont asymmetrically affects lifespan, has unconditional fecundity cost*)


(* ::Subsection::Closed:: *)
(*Calculations*)


resultsEstConditMoreBen = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1, fPI -> 1, fQU -> 1, fQI -> 1, sPU -> 1/2, sPI -> 1, sQU -> 1, sQI -> 3/4, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100],disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsEstConditMoreBen.wdx", resultsEstConditMoreBen];


(* ::Text:: *)
(*Symbiont increases lifespan (establishment probability) in Patch M more than it decreases it in Patch P and has an unconditional cost on host fecundity*)


resultsEstConditMoreBenFecCost = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1, fPI -> 7/8, fQU -> 1, fQI ->7/8, sPU -> 1/2, sPI -> 1, sQU -> 1, sQI -> 3/4, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsEstConditMoreBenFecCost.wdx", resultsEstConditMoreBenFecCost];


(* ::Subsection::Closed:: *)
(*Figures*)


resultsEstConditMoreBen = Import["wdx/resultsEstConditMoreBen.wdx"];
resultsEstConditMoreBenFecCost = Import["wdx/resultsEstConditMoreBenFecCost.wdx"];


plotsEstConditMoreBen = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstConditMoreBen];
plotsEstConditMoreBenFecCost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstConditMoreBenFecCost];


(* ::Text:: *)
(*Symbiont increases lifespan in Patch M more than it decreases it in Patch P, with and without unconditional fecundity cost (Supplement figure 7)*)


Export["Figures/HostControl_EstConditional_MoreBen_FecCost.pdf", Join[Transpose[Map[First, plotsEstConditMoreBen]],Transpose[Map[First, plotsEstConditMoreBenFecCost]]]//
Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Without Uncondition Fecundity Cost", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["With Unconditional Fecundity Cost", 90Degree], Background -> LightGray], 
	Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];
