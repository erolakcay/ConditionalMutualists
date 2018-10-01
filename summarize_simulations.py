# This file is "summarize_simulations.py"
import pandas
import math
import time
import numpy
import re
import csv
import os

# helper function for summarizing simulations
# Input:
#   - simfile: string giving the path to a file with simulation results
#   - patchreplace: function for renaming patches. Renamed patches with the same name will be grouped for
#                   before calculating summary statistics
#   - name: name of the output file for the summary results
# Output:
#   - CSV file with average transmission and infection at each time step for each patch/group of patches
def singlesummary(simfile, patchreplace, name=None):
    simdata = pandas.read_csv(simfile, skipinitialspace=True, na_values="NaN ");
    simdata.patch = map(patchreplace, simdata.patch);
    simsummary = pandas.pivot_table(simdata, values=['symb', 'h', 'v '], index=['time', 'patch'],
    aggfunc={'symb': lambda x: numpy.mean(x) - 1, 'h':lambda x: None if pandas.isnull(x).all() else numpy.mean(x),
     'v ': lambda x: None if pandas.isnull(x).all() else numpy.mean(x)}, dropna=False).stack(dropna=False);
    simsummary = simsummary.rename(name);
    simsummary = simsummary.rename({"v ": "v"})
    return simsummary

# function to summarize simulations
# Input:
#   - directories: list of directories containing simulation files to be analyzed; files may be in subdirectories,
#                  provided all directories and subdirectories contain only simulation results and/or other
#                  subdirectories
#   - outputfile: string giving the name of the file the output should be saved to
#   - patchnamefun: (optional) a function that takes patch numbers and converts them to a name. Patches will
#                   grouped for analysis based on their new names
#   - simnamefun: (optional) a function that takes the full path of a file as input and returns a string. Columns
#                  representing simulation results will be named using the function
# Output:
#    - a CSV file (outputfile) with the average see transmission rates and fraction of infected hosts in each patch
#      or type of patch at each time point
def simulationsummary(directories, outputfile, patchnamefun=None, simnamefun=None):
    # set renaming functions to defaults if not given as input
    if patchnamefun == None:
        # default: patches keep their original names
        patchnamefun = lambda x:x
    if simnamefun == None:
        # default: simulations are not named
        simnamefun = lambda x:None

    # body of the function
    pandas.concat([singlesummary(dirname + "/" + filename, patchnamefun, simnamefun(dirname + "/" + filename))
                   for rootdir in directories for dirname, subdirs, files in os.walk(rootdir)
                   for filename in files], axis=1).to_csv(outputfile)


# for naming patches in simulations with 2 patches
def patchreplaceTwoPatches(patchID):
    if patchID == 1:
        return "P"
    elif patchID == 2:
         return "M"
    return patchID

#### for naming patches in simulations with 3 patches
def patchreplace(regime, patchID):
        if regime == "P-patch 1, M-patches 2 & 3":
            if patchID == 1:
                return "P"
            elif (patchID == 2) or (patchID == 3):
                return "M"
        elif regime == "P-patches 1 & 2, M-patch 3":
            if (patchID == 1) or (patchID == 2):
                return "P"
            elif patchID == 3:
                return "M"
        return patchID


### Example:
# for fx in ("f = 0.5, s = 1.0, m = 1.0", "f = 1.0, s = 0.5, m = 1.0", "f = 1.0, s = 1.0, m = 0.5"):
#     for d in ("0.005", "0.05", "0.5"):
#         simulationsummary(["Host_Control_rep" + rep + "/"
#                    "/" + fx + "/d = " + d + "/N = 200, patches = 2" for
#                      rep in ("a", "b", "c", "d",
#                    "e", "f", "g", "h", "i", "j")],
#                    "Host_Control_" + fx.replace(" = ", "_").replace(".", "-").replace(", ", "_") + "_d_" +
#                   d.replace(".", "-") + "_summary.csv",
#                  patchreplaceTwoPatches, lambda x: "\"" + x[-20:-16] + x[-8:-4] + "\"")
#
