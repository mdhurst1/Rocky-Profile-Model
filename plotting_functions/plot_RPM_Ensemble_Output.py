# -*- coding: utf-8 -*-
"""
Script to plot the results of RPM and RockyCoastCRN Ensembles for sensitivity testing

Martin Hurst,
September 2nd 2019

@author: mhurst
"""

#import modules
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from RPM_CRN_Figure import *

# define workspace
ResultsFolder = Path(r'C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Results\ModelOutput')
PlotsFolder = Path(r'C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Results\Plots')
    
#set parameter values explored
Gradients = [0.5, 1, 0]
SLRs = [-0.001, 0, 0.001]
TidalRanges = [1, 4, 8]
WeatheringRates = [0.01, 0.1, 0.5]
SubtifalEfficacies = [0.001, 0.01, 0.1]
Resistances = [0.1, 1, 10.]
WaveHeights = [1, 2, 3]
WaveAttenuationConst = [0.01, 0.1, 1]

# set up a figure
GradientsFigure = RPM_CRN_Figure(FigWidth_Inches=11.)

# make plot for gradient variation
for i in range(0,len(Gradients)):

    # makes plots
    ProfileName = ("EnsembleShoreProfile_G" + str(Gradients[i])
                            + "_S_" + str(SLRs[1])
                            + "_T_" + str(TidalRanges[1])
                            + "_W_" + str(WeatheringRates[1])
                            + "_Ws_" + str(SubtifalEfficacies[1])
                            + "_R_" + str(Resistances[1])
                            + "_H_" + str(WaveHeights[1])
                            + "_A_" + str(WaveAttenuationConst[1])
                            + ".xz")

    ConcentrationsName = "EnsembleConcentrations_"+ProfileName.lstrip("EnsembleShoreProfile_").rstrip("xz")+"xn"
    
    # Load profile and concentration data
    FigureFile = PlotsFolder / "Gradients.png"
    
    # label
    Label="Initial Gradient " + str(Gradients[i])
    
    # add final results to existing plot
    GradientsFigure.PlotProfileAndConcentrationFigure(ResultsFolder / ProfileName, 
                                                      ResultsFolder / ConcentrationsName, 
                                                      Label=Label, Legend=True)

# save results
GradientsFigure.SaveFig(FigureFile)

"""
for b in range(0,len(TidalRanges)):
    
    
    for c in range(0,len(WaveHeights)):
    
    
    
    for d in range(0,len(WeatheringRates)):
    
    
    for e in range(0,len(Resistances)):
    
    
    for f in range(0,len(BreakingCoefficients)):
    
    
    for g in range(0,len(BrokenCoefficients)):
    
    
    Parameters = [Gradients[a],TidalRanges[b],WaveHeights[c],WeatheringRates[d],Resistances[e],BreakingCoefficients[f],BrokenCoefficients[g]]
                                FileName    = Path + ("ShoreProfile_G" + str(Gradients[a])
                                            + "_T_" + str(TidalRanges[b])
                                            + "_H_" + str(WaveHeights[c])
                                            + "_W_" + str(WeatheringRates[d])
                                            + "_R_" + str(Resistances[e])
                                            + "_Br_" + str(BreakingCoefficients[f])
                                            + "_Bo_" + str(BrokenCoefficients[g])
                                            + ".xz")
                                make_plot(FileName, Parameters)
                                #sys.exit("Temp break")
"""