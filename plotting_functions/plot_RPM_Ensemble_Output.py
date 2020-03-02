# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RPM and RockyCoastCRN

Martin Hurst,
September 2nd 2019

@author: mhurst
"""

#import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
#rc('text', usetex=True)
padding = 5

def make_plot(FileName, Parameters):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,3.3))

    #First plot the morphology through time
    # declare the file and the axis
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    EndTime = float(Lines[-1].strip().split(" ")[0])

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    CliffHeight = Header[0]
    dz = Header[1]
    
    # Only plot every 1 000 years
    PlotTime = 0
    PlotInterval = 100
    EndTime = 1001
    
    ax1 = plt.subplot(111)

    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
        Line = (Lines[j].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[1:],dtype="float64")
        NValues = len(X)
        Z = np.linspace(CliffHeight,-CliffHeight, NValues)
        
        if ((Time >= PlotTime) and (Time < EndTime)):
            ax1.plot(X,Z,'-',lw=1.5,color=cm.coolwarm(Time/EndTime))
            ax1.text(X[0],Z[0]-1,str(np.int(Time))+" years",rotation=-90)
            PlotTime += PlotInterval
    
    plt.text(X[-1]+5,Z[-1]+1,(  "Tr = "+str(Parameters[1])+" m\n"+
                    "H = "+str(Parameters[2])+" m\n"+
                    "W = "+str(Parameters[3])+"\n"+
                    "R = "+str(Parameters[4])+"\n"+
                    "Br = "+str(Parameters[5])+"\n"+
                    "Bo = "+str(Parameters[6])),fontsize=8)

    # tweak the plot
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    xmin, xmax = ax1.get_xlim()
    plt.xlim(xmin-10,xmax+10)
    plt.ylim(-CliffHeight,CliffHeight)
    plt.tight_layout()
    plt.savefig(FileName.rstrip("xz")+"png")
    plt.clf()

if __name__ == "__main__":
    Path = "../results/Hiro_Ensemble_Runs/"
    
    #set parameter values explored
    Gradients = [0, 0.1, 1.]
    SLR = [-0.001, 0, 0.001]
    TidalRanges = [1., 4., 8.]
    WeatheringRates = [0.005, 0.05, 0.5]
    SubtifalEfficacy = [0.001, 0.01, 0.1]
    Resistances = [10., 100., 1000.]
    WaveAttenuationConst = [0.01, 0.1, 1.]

    #Loop across parameter space
    for a in range(0,len(Gradients)):
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