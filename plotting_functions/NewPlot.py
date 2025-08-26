#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 12:19:01 2025

@author: mhurst
"""
# import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc

def ReadShoreProfile(Folder, Filename):

    """
    Function to read shore profile data from file

    MDH, Aug 2025

    """

    # load the profile file
    f = open(Folder+Filename,'r')
    Lines = f.readlines()
    f.close()

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=float)
    dz = Header[0]
    NTimes = (int)((len(Lines)-1)/2)
    print(NTimes)
    
    #Get header info and setup X coord
    Times = []
    
    # loop through lines and append X and Z data to array of vectors
    for i in range(1,len(Lines)-2, 2):
        print(i)
        XLine = Lines[i]
        ZLine = Lines[i+1]
        
        
        XSplitLine = XLine.strip().split(" ")
        ZSplitLine = ZLine.strip().split(" ")
        
        Times.append(float(XSplitLine[1]))

        if i == 1:
            X = np.array(XSplitLine[3:],dtype="float64")
            Z = np.array(ZSplitLine[3:],dtype="float64")
        else:
            X = np.vstack((X,np.array(XSplitLine[3:],dtype="float64")))
            Z = np.vstack((Z,np.array(ZSplitLine[3:],dtype="float64")))

    Times = np.array(Times)
    return Times, X, Z

def PlotShoreProfile(Folder, Filename, PlotInterval=1000.):
    
    """
    Function to plot shore profile data from file

    MDH, Aug 2025
    
    """
    
    # read the data from file
    Times, X, Z = ReadShoreProfile(Folder, Filename)
    
    
    #create blank figure
    fig = plt.figure(1,figsize=(12.,6.))
    ax1 = plt.subplot(111)
    
    # Only plot every so many years
    PlotTime = Times[0]
    StartTime = Times[0]
    EndTime = Times[-1]
    
    #Colourmap
    ColourMap = cm.bone_r
    
    #Loop through times and plot at time interval
    while PlotTime >= EndTime:
        
        # get index
        Index = np.argmin(np.abs(Times-PlotTime))
                    
        if (PlotTime == StartTime):
            ax1.plot(X[Index], Z[Index], 'k--', lw=1., zorder=10, label="Initial Profile")
            plt.show()
            
        elif PlotTime == EndTime:
            ax1.plot(X[-1], Z[-1], 'r-',lw=2., zorder=10, label="Final Profile")
            break
        
        else:
            colour = (PlotTime-StartTime)/(EndTime-StartTime)
            ax1.plot(X[Index], Z[Index], '-', lw=1., color=ColourMap(colour))
        
        PlotTime -= PlotInterval

    # tweak the plot
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    plt.xlim(np.min(X[-1]),np.max(X[-1]))
    plt.ylim(np.min(Z),np.max(Z))
    
    #add the colourbar
    #cbar = fig.colorbar(self.Times)
    #cbar.set_label('Time (years)')

    # add the legend
    plt.tight_layout()
    
    # save fig
    plt.savefig(Folder+Filename.rstrip("xz")+"png",dpi=300)
    
if __name__ == "__main__":
    Folder = '../test/'
    Filename = 'TestProject_ShoreProfile.xz'
    PlotShoreProfile(Folder,Filename)