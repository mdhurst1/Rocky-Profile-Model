# -*- coding: utf-8 -*-
"""
Script to plot the results of RPM and RockyCoastCRN Ensembles for sensitivity testing

Martin Hurst,
September 2nd 2019

@author: mhurst
"""

#import modules
from pathlib import Path
from RPM_CRN_Figure import *

# define workspace
ResultsFolder = Path(r'C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Rocky-Profile-Model\driver_files\Ensemble_Drivers')
    
#set parameter values explored
Gradients = [0.175, 1, 0]
SLRs = [-0.001, 0, 0.001]
TidalRanges = [1, 4, 8]
WeatheringRates = [0.001, 0.01, 0.1]
SubtidalEfficacies = [0.001, 0.01, 0.1]
Resistances = [10, 100, 1000]
WaveHeights = [1, 2, 3]
WaveAttenuationConst = [0.01, 0.1, 1]

for i in range(0,len(Gradients)):
    
    # set up a figure
    SLRFigure = RPM_CRN_Figure(FigWidth_Inches=11.)
    
    # makes plots
    ProfileName = ("EnsembleShoreProfile_G" + str(Gradients[i])
                            + "_S_" + str(SLRs[1])
                            + "_T_" + str(TidalRanges[1])
                            + "_W_" + str(WeatheringRates[1])
                            + "_Ws_" + str(SubtidalEfficacies[1])
                            + "_R_" + str(Resistances[1])
                            + "_H_" + str(WaveHeights[1])
                            + "_A_" + str(WaveAttenuationConst[1])
                            + ".xz")
    
    ConcentrationsName = "EnsembleConcentrations_"+ProfileName.lstrip("EnsembleShoreProfile_").rstrip("xz")+"xn"
    
    # add final results to existing plot
    SLRFigure.PlotProfileEvolutionFigure(ResultsFolder / ProfileName)
    
    # save results
    FigureFile = PlotsFolder / "Gradients" + str(SLRs[i]) + ".png"
    SLRFigure.SaveFig(FigureFile)


    