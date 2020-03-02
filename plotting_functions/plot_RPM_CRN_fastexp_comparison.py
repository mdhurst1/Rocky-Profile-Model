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
rc('ytick.major',pad=1)
rc('xtick.major',pad=1)
padding = 1

#create blank figure
plt.figure(1,figsize=(6.6,9.9))

# spatial plot
ax1 = plt.subplot(311)
plt.axis('equal')

# concentration plot
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)

    
ProfileFileList = ["../driver_files/" + "SlowExp_ShoreProfile.xz", "../driver_files/" + "FastExp_ShoreProfile.xz"]
ConcentrationFileList = ["../driver_files/" + "SlowExp_Concentrations.xn", "../driver_files/" + "FastExp_Concentrations.xn"]
ColourMaps = [cm.bone_r,cm.Reds]

for i in range(0,len(ProfileFileList)):

    #First plot the morphology through time
    # declare the file and the axis
    f = open(ProfileFileList[i],'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    StartTime = float(Lines[1].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])
    f.close()

    f = open(ConcentrationFileList[i],'r')
    ConcLines = f.readlines()
    NoLines = len(ConcLines)
    StartTime = float(Lines[1].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])
    f.close()

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    CliffHeight = Header[0]
    print(CliffHeight)
    dz = Header[2]

    # Only plot every so many years
    StartTime = 10000
    PlotTime = StartTime
    PlotInterval = -1000
    EndTime = 0
    
    #Get header info and setup X coord
    for j in range(2,NoLines-1):
        
        Line = (Lines[j].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[2:],dtype="float64")
        NValues = len(X)
        Z = np.arange(CliffHeight,CliffHeight-NValues*dz,-dz)
                
        if (Time == StartTime) or (Time == -9999):
            print(Time)
            ax1.plot(X+i,Z,'k--',lw=1.,zorder=10,label="Initial Profile")
            PlotTime += PlotInterval
        elif Time <= EndTime:
            ax1.plot(X+i,Z,'k-',lw=1.,zorder=10,label="Final Profile")
            PlotTime += PlotInterval
            break
        elif (Time <= PlotTime):
            print(Time)
            colour =(Time-StartTime)/(EndTime-StartTime)
            ax1.plot(X+i,Z,'-',lw=1.,color=ColourMaps[i]((Time)/(StartTime)))
            PlotTime += PlotInterval

    ax1.plot(X,Z,'k-',lw=1.5)

    N10Line = (ConcLines[-1].strip().split(" "))   #reading line -2?   
    N10 = np.array(N10Line[1:],dtype="float64")
    X2 = np.arange(0,len(N10))*0.1
    mask = [N10!=N10[-1]]
    N10 = N10[mask]
    X2 = X2[mask]
           
    ax2.plot(X2,N10,'-',lw=1.5, color=ColourMaps[i](1.))
    
    if i == 0:
        N10old = N10
    
ax3.plot(X2,N10old-N10,'k-')

# tweak the plot
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xlabel("Distance (m)")
ax1.set_ylabel("Elevation (m)")
ax2.set_ylabel("$^{10}$Be concentration (at g$^{-1}$)")
ax3.set_ylabel("Difference in $^{10}$Be concentration (at g$^{-1}$)")
ax1.set_xlim(0,175)
ax2.set_xlim(0,175)
ax3.set_xlim(0,175)

#plt.ylim(-30,30)
#plt.tight_layout()
plt.savefig('RPM_test_fast_exp.png',dpi=300)

        
