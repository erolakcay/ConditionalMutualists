#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This file is `run_analytical_model_and_all_plots.wls`.*)


(* ::Text:: *)
(*This script runs the analytical model and produces plots for the analytical and simulation results.*)


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


plotsMortSimsHost = ParallelMap[plotSimTransmission[#, 2]&, {"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-005_summary.csv",
"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-05_summary.csv",
"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-5_summary.csv"}];


plotsFecSimsHost = ParallelMap[plotSimTransmission[#, 2]&, {"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


Export["Figures/HostControl_N200_AvgTM_MainText.pdf", 
Join[Transpose[Map[Last, plotsMortSimsHost, {2}]], Transpose[Map[Last, plotsFecSimsHost, {2}]]]//
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


plotsMortSimsHostInfected = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-005_summary.csv",
"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-05_summary.csv",
"Simulation results/Host_Control_f_1-0_s_1-0_m_0-5_d_0-5_summary.csv"}];


plotsFecSimsHostInfected = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


plotsMortSimsHostInfectedRaster = Map[{#[[1]], Rasterize[Show[#[[2]], FrameLabel -> None, FrameTicks -> None], RasterSize -> 300]}&, plotsMortSimsHostInfected, {2}];


plotsFecSimsHostInfectedRaster = Map[{#[[1]], Rasterize[Show[#[[2]], FrameLabel -> None, FrameTicks -> None], RasterSize -> 300]}&, plotsFecSimsHostInfected, {2}];


(* ::Text:: *)
(*Fraction of infected hosts at different simulation endpoints (Supplement figure 8)*)


Export["Figures/HostControl_N200_Infection_Supplement.pdf", 
Join[Transpose[Map[Last, plotsMortSimsHostInfectedRaster, {2}]], Transpose[Map[Last, plotsFecSimsHostInfectedRaster, {2}]]]//
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
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Subsection:: *)
(*Host control, population size 1000 (supplement)*)


plotsEstSimsHostInfectedN1000 = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_N_1000_f_1-0_s_0-5_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_N_1000_f_1-0_s_0-5_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_N_1000_f_1-0_s_0-5_m_1-0_d_0-5_summary.csv"}];


plotsFecSimsHostInfectedN1000 = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_N_1000_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_N_1000_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_N_1000_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


plotsEstSimsHostInfectedN1000Raster = Map[{#[[1]], Rasterize[Show[#[[2]], FrameLabel -> None, FrameTicks -> None], RasterSize -> 300]}&, plotsEstSimsHostInfectedN1000, {2}];


plotsFecSimsHostInfectedN1000Raster = Map[{#[[1]], Rasterize[Show[#[[2]], FrameLabel -> None, FrameTicks -> None], RasterSize -> 300]}&, plotsFecSimsHostInfectedN1000, {2}];


(* ::Text:: *)
(*Simulation of host evolution in larger population (Supplement figure 9)*)


Export["Figures/HostControl_N1000_Infection_Supplement.pdf", 
Join[Transpose[Map[Last, plotsEstSimsHostInfectedN1000Raster, {2}]], Transpose[Map[Last, plotsFecSimsHostInfectedN1000Raster, {2}]]]//
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
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Subsection:: *)
(*Symbiont control, population size 200 (supplement)*)


plotsEstSimsSymb = ParallelMap[plotSimTransmission[#, 2]&, {"Simulation results/Symbiont_Control_f_1-0_s_0-5_m_1-0_d_0-005_summary.csv",
"Simulation results/Symbiont_Control_f_1-0_s_0-5_m_1-0_d_0-05_summary.csv",
"Simulation results/Symbiont_Control_f_1-0_s_0-5_m_1-0_d_0-5_summary.csv"}];


plotsFecSimsSymb = ParallelMap[plotSimTransmission[#, 2]&, {"Simulation results/Symbiont_Control_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Symbiont_Control_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Symbiont_Control_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


(* ::Text:: *)
(*Simulation of symbiont evolution (Supplement figure 11)*)


Export["Figures/SymbiontControlSimulations.pdf", 
Join[Transpose[Map[Last, plotsEstSimsSymb, {2}]], Transpose[Map[Last, plotsFecSimsSymb, {2}]]]//
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
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Subsection:: *)
(*Host control, unequal numbers of M and P patches*)


plotsEstSimsHostInfectedMoreM = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_More_M-patches_f_1-0_s_0-5_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_More_M-patches_f_1-0_s_0-5_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_More_M-patches_f_1-0_s_0-5_m_1-0_d_0-5_summary.csv"}];


plotsFecSimsHostInfectedMoreM = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_More_M-patches_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_More_M-patches_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_More_M-patches_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


plotsEstSimsHostInfectedMoreP = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_More_P-patches_f_1-0_s_0-5_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_More_P-patches_f_1-0_s_0-5_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_More_P-patches_f_1-0_s_0-5_m_1-0_d_0-5_summary.csv"}];


plotsFecSimsHostInfectedMoreP = ParallelMap[plotSimTransmissionAndInfection[#, 2]&, {"Simulation results/Host_Control_More_P-patches_f_0-5_s_1-0_m_1-0_d_0-005_summary.csv",
"Simulation results/Host_Control_More_P-patches_f_0-5_s_1-0_m_1-0_d_0-05_summary.csv",
"Simulation results/Host_Control_More_P-patches_f_0-5_s_1-0_m_1-0_d_0-5_summary.csv"}];


(* ::Text:: *)
(*Unequal M and P patches when symbiont affects lifespan (Supplement figure 14)*)


Export["MoreMorPpatches_establishment.pdf",
Join[Transpose[plotsEstSimsHostInfectedMoreM], Transpose[plotsEstSimsHostInfectedMoreP]] //
Grid[Join[{{ "", "", "", Item["d = 0.0067", Background -> Lighter[Gray]], Item["d = 0.067", Background -> Lighter[Lighter[Gray]]], Item["d = 0.67", Background -> LightGray]}},
Join[{{Item[Rotate["M-patches more common", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["M-Patches", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["P-Patch", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["P-patches more common", 90Degree], Background -> LightGray], 
	Item[Rotate["M-Patch", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["P Patches", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Text:: *)
(*Unequal M and P patches when symbiont affects fecundity (Supplement figure 15)*)


Export["MoreMorPpatches_fecundity.pdf",
Join[Transpose[plotsFecSimsHostInfectedMoreM], Transpose[plotFecSimsHostInfectedMoreP]] //
Grid[Join[{{ "", "", "", Item["d = 0.0067", Background -> Lighter[Gray]], Item["d = 0.067", Background -> Lighter[Lighter[Gray]]], Item["d = 0.67", Background -> LightGray]}},
Join[{{Item[Rotate["M-patches more common", 90Degree], Background -> Lighter[Gray]],  
	Item[Rotate["M-Patches", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["P-Patch", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["P-patches more common", 90Degree], Background -> LightGray], 
	Item[Rotate["M-Patch", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["P Patches", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 300];


(* ::Section::Closed:: *)
(*Symbiont affects one fitness component symmetrically*)


(* ::Subsection:: *)
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


resultsFecAloneOtherDispersalRates = ParallelTable[Module[{eq, coord, gHost, gSymb, 
params = {c  -> 1, fPU -> 1/2, fPI -> 1, fQU -> 1, fQI -> 1/2, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 1, mPI -> 1, mQU -> 1, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/10, 25/100, 4/10}}];
Export["wdx/resultsFecAloneOtherDispersalRates.wdx", resultsFecAloneOtherDispersalRates];


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
params = {c  -> 1, fPU -> 1, fPI -> 1, fQU -> 1, fQI -> 1, sPU -> 1, sPI -> 1, sQU -> 1, sQI -> 1, mPU -> 1, mPI -> 1/2, mQU -> 1/2, mQI -> 1}}, 
{eq, coord} = ReleaseHold[Hold[stableEquilTable[Range[0, 1, 1/100], Range[0, 1, 1/100], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, 50]]/.params];
gHost = ReleaseHold[Hold[dMutGrowthTableHost[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
gSymb = ReleaseHold[Hold[dMutGrowthTableSymbiont[N[eq], N[coord], disp, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI]]/.params];
Join[Prepend[params, d -> disp], {equil  -> eq, coordEquil -> coord, dGrowthHost -> gHost, dGrowthSymb -> gSymb}]], {disp, {1/200, 1/20, 1/2}}];
Export["wdx/resultsMortAlone.wdx", resultsMortAlone];


(* ::Subsection:: *)
(*Figures*)


resultsFecAlone = Import["wdx/resultsFecAlone.wdx"];
resultsEstAlone = Import["wdx/resultsEstAlone.wdx"];
resultsMortAlone = Import["wdx/resultsMortAlone.wdx"];
resultsFecAloneOtherDispersalRates = Import["wdx/resultsFecAloneOtherDispersalRates.wdx"];


plotsFecAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAlone];
plotsFecAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsFecAlone];
plotsEstAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstAlone];
plotsEstAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsEstAlone];
plotsMortAloneHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsMortAlone];
plotsMortAloneSymb = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthSymb]]]]/.#]&, resultsMortAlone];
plotsFecAloneOtherDispersalRatesHost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAloneOtherDispersalRates];


(* ::Text:: *)
(*Add marks for equilibria to the plots*)


markedPlotsEstAloneHost = MapIndexed[Show[#1,
ListPlot[{{1.01, -0.01}, If[#2[[1]] != 3, {-0.01, 1.01}]}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 12]]],
PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]&, 
plotsEstAloneHost, {3}];


markedPlotsMortAloneHost = MapIndexed[Show[#1,
ListPlot[{{1.01, -0.01}, If[#2[[1]] != 3, {-0.01, 1.01}]}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 12]]],
PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]&, 
plotsMortAloneHost, {3}];


getHstar[coordEquil_, dGrowthHost_, equilNumber_:1, estimate_:0.7] := Module[{vDerivs, vDerivsFun, vCutoff, growthRates = Map[#[[equilNumber]]&, dGrowthHost, {2}]},
vDerivs = Flatten[Table[{coordEquil[[hind, vind]], growthRates[[hind, vind, 2]]}, {hind, Dimensions[coordEquil][[1]]}, {vind, Dimensions[coordEquil][[2]]}],1];
vDerivsFun = Interpolation[vDerivs, InterpolationOrder -> 1];
vCutoff = h/.FindRoot[vDerivsFun[h, 1], {h, estimate}];
Return[vCutoff]]


markedPlotsFecAloneHost = MapIndexed[Module[{tmpPlot},
If[#2[[1]] == 3, 
tmpPlot = #1,
tmpPlot = Module[{hstar = getHstar[coordEquil/. resultsFecAlone[[#2[[1]]]], dGrowthHost/.resultsFecAlone[[#2[[1]]]]]},
     Show[#1, ListPlot[{{-0.01, 1.01}}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 12]]], 
     ParametricPlot[{1.02, v}, {v, 0, 1.02}, PlotStyle -> Directive[Thickness[0.025], Black]],
     ParametricPlot[{h, 1.02}, {h, hstar, 1}, PlotStyle -> Directive[Thickness[0.025], Black]]]]]; 
Show[tmpPlot, PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]]&, 
plotsFecAloneHost, {3}];


markedPlotsFecAloneOtherDispersalRatesHost = MapIndexed[Module[{tmpPlot},
tmpPlot = Module[{hstar = getHstar[coordEquil/. resultsFecAlone[[#2[[1]]]], dGrowthHost/.resultsFecAlone[[#2[[1]]]]]},
     Show[#1, ListPlot[{{-0.01, 1.01}}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 12]]], 
     ParametricPlot[{1.02, v}, {v, 0, 1.02}, PlotStyle -> Directive[Thickness[0.025], Black]],
     ParametricPlot[{h, 1.02}, {h, hstar, 1}, PlotStyle -> Directive[Thickness[0.025], Black]]]]; 
Show[tmpPlot, PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]]&, 
plotsFecAloneOtherDispersalRatesHost, {3}];


markedPlotsFecAloneSymb = MapIndexed[Show[#1,
ParametricPlot[{h, 1.02}, {h, 0, 1.02}, PlotStyle -> Directive[Thickness[0.025], Black]],
PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]&, 
plotsFecAloneSymb, {3}];


markedPlotsEstAloneSymb = MapIndexed[Show[#1,
ParametricPlot[{h, 1.02}, {h, 0, 1.02}, PlotStyle -> Directive[Thickness[0.025], Black]],
PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]&, 
plotsEstAloneSymb, {3}];


markedPlotsMortAloneSymb = MapIndexed[Show[#1,
ParametricPlot[{h, 1.02}, {h, 0, 1.02}, PlotStyle -> Directive[Thickness[0.025], Black]],
PlotRangeClipping -> False, ImagePadding -> {{All, 3}, {All, 3}}]&, 
plotsMortAloneSymb, {3}];


(* ::Text:: *)
(*Host evolution when symbiont affects host mortality or fecundity (Main text figure 3)*)


Export["Figures/HostControl_MainScenario.pdf",
Join[Transpose[Last /@ markedPlotsMortAloneHost], Transpose[Last /@ markedPlotsFecAloneHost]]//
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
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 600];


(* ::Text:: *)
(*Host evolution at alternate ecological equilibria when symbiont affects mortality (Supplement figure 2)*)


Export["Figures/HostControlMortAlone_extraEquil.pdf", Transpose[Map[First, markedPlotsMortAloneHost[[1;;2]]]]//
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
(*Host evolution when symbiont affects host establishment probability (Supplement figure 3)*)


Export["Figures/HostControlEstAlone.pdf", Join[Transpose[Map[Last, markedPlotsEstAloneHost]], Transpose[markedPlotsEstAloneHost[[1;;2, 1]]]]//
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
(*Host evolution when symbiont affects host fecundity for more dispersal rates (Supplement figure 4)*)


Export["Figures/HostControlFecAloneOtherDispersalRates.pdf", Transpose[markedPlotsFecAloneOtherDispersalRatesHost[[All, 1]]] //
Grid[Join[{{ "", "", "", Item["d = 0.1", Background -> Lighter[Gray]], Item["d = 0.25", Background -> Lighter[Lighter[Gray]]], Item["d = 0.4", Background -> LightGray]}}, 
Join[{{Item[Rotate["Symbiont Affects Fecundity", 90Degree], Background -> Lighter[Gray]],
Item[Rotate["Patch M", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Patch P", 90 Degree], Background -> Lighter[LightGray]], Rotate["Vert. Trans. Prob.", 90 Degree]}}, 
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Fraction of hosts that are infected", myEquilLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];


(* ::Text:: *)
(*Symbiont evolution when symbiont affects host mortality or fecundity (Supplement figure 12)*)


Export["Figures/SymbiontControl_MainScenario.pdf",
Join[Transpose[First /@ markedPlotsMortAloneSymb], Transpose[First /@ markedPlotsFecAloneSymb]]//
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
resultsMortAlone = Import["wdx/resultsMortAlone.wdx"];


conflictFec = Map[MapThread[{#1[[1]], #1[[2]], If[(#2[[1]] == 0) && (#2[[2]] == 0), ,
If[(Norm[#3] == 0) || (Norm[#4] == 0), 0, (#3 . #4)/(Norm[#3] * Norm[#4]), -0.99]]}&,
{Flatten[(coordEquil/. #), 1], Chop[Flatten[Map[First, (equil/. #), {2}], 1]], Chop[Flatten[Map[First, (dGrowthHost/. #), {2}], 1]], 
Chop[Flatten[Map[First, (dGrowthSymb/. #), {2}], 1]]}]&, resultsFecAlone];


conflictEst = Map[MapThread[{#1[[1]], #1[[2]], If[(#2[[1]] == 0) && (#2[[2]] == 0), ,
If[(Norm[#3] == 0) || (Norm[#4] == 0), 0, (#3 . #4)/(Norm[#3] * Norm[#4]), -0.99]]}&,
{Flatten[(coordEquil/. #), 1], Chop[Flatten[Map[First, (equil/. #), {2}], 1]], Chop[Flatten[Map[First, (dGrowthHost/. #), {2}], 1]], 
Chop[Flatten[Map[First, (dGrowthSymb/. #), {2}], 1]]}]&, resultsEstAlone];


conflictMort = Map[MapThread[{#1[[1]], #1[[2]], If[(#2[[1]] == 0) && (#2[[2]] == 0), ,
If[(Norm[#3] == 0) || (Norm[#4] == 0), 0, (#3 . #4)/(Norm[#3] * Norm[#4]), -0.99]]}&,
{Flatten[(coordEquil/. #), 1], Chop[Flatten[Map[First, (equil/. #), {2}], 1]], Chop[Flatten[Map[First, (dGrowthHost/. #), {2}], 1]], 
Chop[Flatten[Map[First, (dGrowthSymb/. #), {2}], 1]]}]&, resultsMortAlone];


(* ::Subsection:: *)
(*Figure*)


myConflictColors = ColorData["Candy"][(11/20) + (9#/20)]&;
myConflictLegend = BarLegend[{myConflictColors, {-1, 1}}, {-1, -0.6, -0.3, 0, 0.3, 0.6, 1}, 
LegendLayout -> "Row",LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 300];


conflictPlots = Map[ListContourPlot[#, BaseStyle -> myFontStyle, ImageSize -> 150, PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ClippingStyle -> White,
ColorFunctionScaling -> False, ColorFunction -> myConflictColors, Contours -> {-0.6, -0.3, 0, 0.3, 0.6}, ContourStyle -> {Black, Black, {Dashed, Thick}, Black, Black}]&, 
{conflictMort, conflictFec}, {2}];


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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Figure*)


resultsFecAloneMoreCosts = Import["wdx/resultsFecAloneMoreCosts.wdx"];
resultsFecAloneMoreBen = Import["wdx/resultsFecAloneMoreBen.wdx"];


plotsFecAloneMoreCosts = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAloneMoreCosts];
plotsFecAloneMoreBen = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecAloneMoreBen];


(* ::Text:: *)
(*Symbiont affects fecundity asymmetrically (supplement figure 5)*)


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


(* ::Subsection:: *)
(*Figures*)


resultsFecWeakEstSameDirection = Import["wdx/resultsFecWeakEstSameDirection.wdx"];
resultsFecWeakEstOppositeDirection = Import["wdx/resultsFecWeakEstOppositeDirection.wdx"];
resultsEstWeakFecSameDirection = Import["wdx/resultsEstWeakFecSameDirection.wdx"];
resultsEstWeakFecOppositeDirection = Import["wdx/resultsEstWeakFecOppositeDirection.wdx"];


plotsFecWeakEstSameDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecWeakEstSameDirection];
plotsFecWeakEstOppositeDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsFecWeakEstOppositeDirection];
plotsEstWeakFecSameDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstWeakFecSameDirection];
plotsEstWeakFecOppositeDirection = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstWeakFecOppositeDirection];


(* ::Text:: *)
(*Conditional effects act in same direction (supplement figure 6)*)


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
(*Conditional effects act in opposite directions (supplement figure 7)*)


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
(*Symbiont affects fecundity conditionally, has unconditional lifespan cost (Supplement figure 8)*)


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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Figures*)


resultsEstConditMoreBen = Import["wdx/resultsEstConditMoreBen.wdx"];
resultsEstConditMoreBenFecCost = Import["wdx/resultsEstConditMoreBenFecCost.wdx"];


plotsEstConditMoreBen = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstConditMoreBen];
plotsEstConditMoreBenFecCost = ParallelMap[ReleaseHold[Hold[plotAD[N[equil], N[coordEquil], Chop[N[dGrowthHost]]]]/.#]&, resultsEstConditMoreBenFecCost];


(* ::Text:: *)
(*Symbiont increases lifespan in Patch M more than it decreases it in Patch P, with and without unconditional fecundity cost (Supplement figure 9)*)


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


(* ::Section:: *)
(*Basic reproductive number*)


(* ::Subsection::Closed:: *)
(*Calculations*)


resultsFecAlone = Import["wdx/resultsFecAlone.wdx"];
resultsEstAlone = Import["wdx/resultsEstAlone.wdx"];
resultsMortAlone = Import["wdx/resultsMortAlone.wdx"];


(* ::Text:: *)
(*Functions for calculating R0 for an infection in each patch*)


individualR0P[h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=
(mPU/(2 mPI ((1-d) fPU + d fQU) sPU)) (v fPI ((1-d) sPI + d sQI) + h c ((1 - d) fPU + d fQU) sPI)


individualR0Q[h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=
(mQU/(2 mQI ((1-d) fQU + d fPU) sQU)) (v fQI ((1-d) sQI + d sPI) + h c ((1 - d) fQU + d fPU) sQI)


(* ::Text:: *)
(*Function for calculating the expected R0, where R0P and R0Q are weighted by the fraction of the total infections that occur in their patch*)


expectedR0[iP_, iQ_,h_, v_, d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_] :=
If[(iP > 0) && (iQ > 0),
(iP/(iP + iQ)) individualR0P[h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI] +
(iQ/(iP + iQ)) individualR0Q[h, v, d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI],
0]


(* ::Text:: *)
(*Calculate R0 for each symbiont effect when it is alone and symmetric*)


r0FecAlone = Table[ReleaseHold[Hold[Table[expectedR0[equil[[hind, vind,eqind, 1]],equil[[ hind, vind,eqind, 2]], coordEquil[[hind, vind, 1]], coordEquil[[hind, vind, 2]], d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], {hind, Length[equil]}, {vind, Length[equil[[hind]]]},{eqind, Length[equil[[hind, vind]]]}]]/. resultsFecAlone[[i]]], {i, Length[resultsFecAlone]}];
r0EstAlone = Table[ReleaseHold[Hold[Table[expectedR0[equil[[hind, vind,eqind, 1]],equil[[ hind, vind,eqind, 2]], coordEquil[[hind, vind, 1]], coordEquil[[hind, vind, 2]], d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], {hind, Length[equil]}, {vind, Length[equil[[hind]]]},{eqind, Length[equil[[hind, vind]]]}]]/. resultsEstAlone[[i]]], {i, Length[resultsEstAlone]}];
r0MortAlone = Table[ReleaseHold[Hold[Table[expectedR0[equil[[hind, vind,eqind, 1]],equil[[ hind, vind,eqind, 2]], coordEquil[[hind, vind, 1]], coordEquil[[hind, vind, 2]], d, c, fPU, fPI, fQU, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI], {hind, Length[equil]}, {vind, Length[equil[[hind]]]},{eqind, Length[equil[[hind, vind]]]}]]/. resultsMortAlone[[i]]], {i, Length[resultsMortAlone]}];


(* ::Subsection:: *)
(*Plots*)


myColorsR0[f_] := ColorData["Candy"][f/Max[{Max[r0FecAlone], Max[r0EstAlone], Max[r0MortAlone]}]];


myR0Legend = BarLegend[{myColorsR0[#]&,{0,Max[{Max[r0FecAlone], Max[r0EstAlone], Max[r0MortAlone]}]}}, LegendLayout -> "Row",LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 400];


r0FecAlonePlots = Table[plotR0altColor[(equil/. resultsFecAlone[[i]]),(coordEquil/. resultsFecAlone[[i]]), r0FecAlone[[i]], myColorsR0, myFontStyle, False], {i, Length[resultsFecAlone]}];


r0EstAlonePlots = Table[plotR0altColor[(equil/. resultsEstAlone[[i]]),(coordEquil/. resultsEstAlone[[i]]), r0EstAlone[[i]], myColorsR0, myFontStyle, False], {i, Length[resultsEstAlone]}];


r0MortAlonePlots = Table[plotR0altColor[(equil/. resultsMortAlone[[i]]),(coordEquil/. resultsMortAlone[[i]]), r0MortAlone[[i]], myColorsR0, myFontStyle, False], {i, Length[resultsMortAlone]}];


Export["r0.pdf", {r0EstAlonePlots,r0MortAlonePlots,r0FecAlonePlots}//Grid[Join[{{ "", "", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Symbiont Affects...", 90Degree], Background -> Lighter[LightGray]],  
	Item[Rotate["Establishment", 90 Degree], Background -> Lighter[Gray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, Item[Rotate["Mortality", 90 Degree], Background -> Lighter[Lighter[Gray]]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{SpanFromAbove, 
	Item[Rotate["Fecundity", 90 Degree], Background -> LightGray], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "", "", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{Subscript["R", 0], myR0Legend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&, ImageResolution -> 600];


(* ::Section:: *)
(*Comparing establishment and mortality effects*)


(* ::Text:: *)
(*Find the magnitude of the growth rate derivative for establishment and mortality effects; calculate establishment/mortality to compare the magnitudes*)


resultsEstAlone = Import["wdx/resultsEstAlone.wdx"];
resultsMortAlone = Import["wdx/resultsMortAlone.wdx"];


estVsMortMagnitudeHost = Table[ListDensityPlot[Flatten[MapThread[{#1[[1]], #1[[2]], If[Norm[First[#2]] == 0, , If[Norm[Last[#4]] == 0, 1, Norm[Last[#3]]/Norm[Last[#4]]]]}&, 
{(coordEquil/. resultsEstAlone[[index]]), (equil/.resultsEstAlone[[index]]),
(dGrowthHost/.resultsEstAlone[[index]]), (dGrowthHost/.resultsMortAlone[[index]])}, 2], 1], 
ColorFunction -> plottingFunctions`Private`myColors, PlotRange -> {{0, 1}, {0, 1}, {0, 1}}, LabelStyle -> myFontStyle,
ColorFunctionScaling -> False, ScalingFunctions -> {#^4&, #^(1/4)&}], {index, Length[resultsEstAlone]}];


estVsMortMagnitudeSymb = Table[ListDensityPlot[Flatten[MapThread[{#1[[1]], #1[[2]], If[Norm[First[#2]] == 0, , If[Norm[Last[#4]] == 0, 1, Norm[Last[#3]]/Norm[Last[#4]]]]}&, 
{(coordEquil/. resultsEstAlone[[index]]), (equil/.resultsEstAlone[[index]]),
(dGrowthSymb/.resultsEstAlone[[index]]), (dGrowthSymb/.resultsMortAlone[[index]])}, 2], 1], 
ColorFunction -> plottingFunctions`Private`myColors, PlotRange -> {{0, 1}, {0, 1}, {0, 1}}, LabelStyle -> myFontStyle,
ColorFunctionScaling -> False, ScalingFunctions -> {#^4&, #^(1/4)&}], {index, Length[resultsEstAlone]}];


estVsMortMagnitudeLegend = BarLegend[{plottingFunctions`Private`myColors[#]&, {0, 1}}, ScalingFunctions -> {#^4&, #^(1/4)&}, Ticks -> {0, 0.5, 0.7, 0.8, 0.9, 1},
LabelStyle -> myFontStyle, Method -> {AxesStyle -> Black, TicksStyle -> Black}, LegendMarkerSize -> 300, LegendLayout -> "Row"];


Export["Figures/EstVsMort.pdf",
Join[{estVsMortMagnitudeHost, estVsMortMagnitudeSymb}]//
Grid[Join[{{"", "", Item["d = 0.005", Background -> Lighter[Gray]], Item["d = 0.05", Background -> Lighter[Lighter[Gray]]], 
Item["d = 0.5", Background -> LightGray]}},
Join[{{Item[Rotate["Host Control", 90 Degree], Background -> Lighter[Gray]], Rotate["Vert. Trans. Prob.", 90 Degree]}, 
{Item[Rotate["Symbiont Control", 90 Degree], Background -> LightGray], Rotate["Vert. Trans. Prob.", 90 Degree]}},  
Module[{pos = 1}, Map[Labeled[Show[#, ImageSize -> 130, FrameLabel -> None, FrameTicks -> {{{0, 0.2, 0.4, 0.6, 0.8, 1}, None}, {{0, 0.2, 0.4, 0.6, 0.8, 1}, None}}], 
Alphabet[][[pos++]] <> ")", {{Left, Top}}]&, #, {2}]], 2], 
{{"", "",  "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob.", "          Horiz. Trans. Prob."}},
{{"", "", Grid[{{"Ratio of derivatives of growth rate\n(establishment/mortality)", 
estVsMortMagnitudeLegend}}], SpanFromLeft, SpanFromLeft}}], 
Alignment -> {Center, Center}, BaseStyle -> myFontStyle]&,  ImageResolution -> 600];
