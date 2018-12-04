# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

@author: mhurst
"""

# import plotting tools and set the back end for running on server
#import matplotlib
#matplotlib.use('Agg')

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

def make_plot(FileName,ColourMap):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,3.3))

    #First plot the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    StartTime = float(Lines[1].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    CliffHeight = Header[0]
    dz = Header[1]
    
    # Only plot every so many years
    StartTime = 9000
    PlotTime = StartTime
    PlotInterval = -1000
    EndTime = 0
    
    ax1 = plt.subplot(111)
    plt.axis('equal')
    
    #Colourmap
    ColourMap = cm.bone_r
    
    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
        Line = (Lines[j].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[1:],dtype="float64")
        NValues = len(X)
        Z = np.linspace(CliffHeight,-CliffHeight, NValues)
        
        if Time >= EndTime:
            ax1.plot(X,Z,'k-',lw=1.,zorder=10,label="Final Profile")
            PlotTime += PlotInterval
            break
        elif Time == StartTime:
            print(Time)
            ax1.plot(X,Z,'k--',lw=1.,zorder=10,label="Initial Profile")
            PlotTime += PlotInterval
        elif (Time <= PlotTime):
            print(Time)
            colour =(Time-StartTime)/(EndTime-StartTime)
            ax1.plot(X,Z,'-',lw=1.,color=ColourMap(colour))
            PlotTime += PlotInterval

    
    ax1.plot(X,Z,'k-',lw=1.5)
            
        
    # tweak the plot
    #ax1.set_xticklabels([])
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    #plt.xlim(np.min(X),np.max(X))
    plt.ylim(-CliffHeight/2,CliffHeight/2)
    #plt.ylim(-30,30)
    plt.tight_layout()
    plt.savefig('RPM_test.png',dpi=300)

if __name__ == "__main__":
    FileName = "../driver_files/"
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        
