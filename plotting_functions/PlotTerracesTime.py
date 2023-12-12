# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 18:41:47 2023

@author: mh322u
"""

#import modules
from Find_Terraces import *
from RPM_CRN_Analysis_Functions import *
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import numpy as np
import sys, os
from tqdm import tqdm
import pdb

# dont plot in spyder pls
%matplotlib qt5
#plt.ion()

# figure properties
# Set up fonts for plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 16

def AnimateTerraceTimeseries(Folder, RunID, MinWidth=3.):
    
    # read profile file
    ResultsFolder = Folder + "Results/"
    PlotsFolder = Folder + "Plots/"
    
    # load the profile data
    ProfileFileName = str(RunID) + "_ShoreProfile.xz"
    Times, SeaLevels, Z, X = ReadShoreProfile(ResultsFolder+ProfileFileName)
    FinalX = X[-1]
    dZ = np.round(Z[1]-Z[0],1)
    
    StartTime = Times[0]
    EndTime = Times[-1]
    Time = StartTime
    TimeInterval = 1000
    OldIndex = -9999
    
    if StartTime > EndTime:
        TimeInterval *= -1
    
    # Set up the figure and axis
    fig = plt.figure(1,(16,9))
    ax = fig.add_subplot(111)
    ax.set_ylim(Z[-1],Z[0])
    ax.set_xlim(FinalX[-1], FinalX[0])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    plt.tight_layout()
    
    # Create an empty line to be animated
    Profile, = ax.plot([], [], 'k-', lw=2)
    Segments = LineCollection([], colors='r', lw=2)
    ax.add_collection(Segments)
    
    # define an init function to set up the beginning of each frame for rendering
    def init():
        Profile.set_xdata([])
        Segments.set_segments([])
        return [Profile, Segments]

    # Function to update the plot in each animation frame
    def update(frame):
        
        # Update the profile line with new X
        ThisX = X[frame]
        Profile.set_data(ThisX, Z)
        
        # find the terraces
        TerracesDF = FindTerracesRaw(ThisX, Z)
        
        # plot the terraces
        TerraceSegments = []
        for i,Terrace in TerracesDF.iterrows():
            
            # make sure indices are ints
            Ind1 = Terrace.StartIndex
            Ind2 = Terrace.EndIndex
            X_Segment = ThisX[Ind1:Ind2]
            Z_Segment = Z[Ind1:Ind2]
            
            TerraceSegments.append(np.column_stack((X_Segment, Z_Segment)))
        
        # set the line collection
        Segments.set_segments(TerraceSegments)
        
        return [Profile, Segments]
   
    # Create the animation
    animation = FuncAnimation(fig, update, frames=len(Times), init_func=init, repeat=False)
    
    # Save the animated plot as a GIF
    animation.save(PlotsFolder+'animated_plot.gif', writer='pillow', fps=5)
    
    plt.show()       
    
if __name__ == "__main__":
    
    # setup workspace
    #Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
    Folder = "C:/Users/mh322u/OneDrive/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
    
    # RunID = Record.RunID.values[0]
    RunID = 2198
    AnimateTerraceTimeseries(Folder, RunID)