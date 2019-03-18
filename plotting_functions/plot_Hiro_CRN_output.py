# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

def make_plot(FileName,ColourMap):
    
    #create blank figure
    plt.figure(1,figsize=(8,4))

    #First load the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"_ShoreProfile.xz"
    f = open(ProfileName,'r')
    MorphLines = f.readlines()
    NoLines = len(MorphLines)
    StartTime = float(MorphLines[2].strip().split(" ")[0])
    EndTime = float(MorphLines[-1].strip().split(" ")[0])
    f.close()
    
    # Get Z Values
    HeaderLine = MorphLines[0].strip().split(" ")
    CliffHeight = float(HeaderLine[0])
    MinElev = float(HeaderLine[1])
    dZ = float(HeaderLine[2])
    NValues = (int)((CliffHeight-MinElev)/dZ+1)
    Z = np.linspace(CliffHeight, MinElev, NValues)
    
    #Second load CRN concentrations through time
    # declare the file and the axis
    ProfileName = FileName+"Concentrations.xn"
    f = open(ProfileName,'r')
    NLines = f.readlines()
    #NNoLines = len(NLines)
    #NEndTime = float(NLines[-1].strip().split(" ")[0])
    f.close()
    
    # Only plot every 1 000 years
    PlotTime = 9900
    PlotInterval = 100

    
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    
    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
        MorphLine = (MorphLines[j].strip().split(" "))
        
        Time = int(MorphLine[0])
        RSL = float(MorphLine[1])
        
        #Read morphology
        X = np.array(MorphLine[2:],dtype="float64")
                
        if (j == 1):
            ax1.plot(X,Z,'--',lw=1.5,color=ColourMap((Time)/(StartTime)))
        if (Time == PlotTime):
            print (Time,StartTime)
            ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time)/(StartTime)))
            PlotTime -= PlotInterval
    
    N10Line = (NLines[-2].strip().split(" "))
    N10 = np.array(N10Line[1:],dtype="float64")
    X2 = np.arange(0,len(N10))*0.1
    mask = [N10!=N10[-1]]
    N10 = N10[mask]
    X2 = X2[mask]
           
    ax2.plot(X2,N10/1000.,'k-',lw=1.5)
        
    # tweak the plot
    #ax1.set_xticklabels([])
    
    ax1.set_ylabel("Elevation (m)")
    ax2.set_ylabel("Concentration (x 10$^3$ a g$^-1$)")
    #xmin, xmax = ax1.get_xlim()
    #ax1.set_xlim(xmin-10,xmax+10)
    #ax2.set_xlim(xmin-10,xmax+10)
    #ax1.set_ylim(-10,10)
    
    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('CB_Test_2.png',dpi=300)

if __name__ == "__main__":
    FileName = "../../RPM_JRS/CB_Test_2" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        