# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 01:09:55 2023

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

%matplotlib qt5
#def TrackTerraces(Folder, RunID, MinWidth=3.):

# setup workspace
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"

# RunID = Record.RunID.values[0]
RunID = 2198
    
# get all the profile data
ResultsFolder = Folder + "Results/"
PlotsFolder = Folder + "Plots/"
AnalysisFolder = Folder + "Analysis/"

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

# similarity threshold as a function of mean elevation 
# do i somehow need to factor in RSL change here too, might help with autogenic!
SimilarThresh = 2.

TerracesTSList = []
print(TerracesTSList)
Columns = ["TerraceID","StartIndex","EndIndex","Width","MeanElev","ElevChange","Slope","Time"]
BlankTerraceDF = pd.DataFrame(columns = Columns)

# loop through time
for i, Time in enumerate(Times):

    # identify a terrace and track its progress
    TerracesDF = FindTerracesRaw(X[i],Z)
    
    print(Time)
    
    if Time == 4900:
        
        ThisX = X[i]
        plt.plot(ThisX,Z,'k-')
        
        for i, Terrace in TerracesDF.iterrows():
            # make sure indices are ints
            Ind1 = Terrace.StartIndex
            Ind2 = Terrace.EndIndex
            plt.plot(ThisX[Ind1:Ind2],Z[Ind1:Ind2],'r-', lw=2, zorder=9)
        
        plt.show()
        
        pdb.set_trace()
        
    if (len(TerracesDF) == 0):
        continue
    
    # first time out no old terraces to compare to
    # or if no terraces have previously been found
    if (len(TerracesTSList) == 0):
        # loop through terraces, probs only one
        for i, Terrace in TerracesDF.iterrows():
            
            # add time to DF
            Terrace["Time"] = Time
            
            # create the new DF for this new terrace and add the terrace to it
            TerracesTSList.append(BlankTerraceDF)
            TerracesTSList[Terrace["TerraceID"]] = pd.concat([TerracesTSList[Terrace["TerraceID"]], Terrace.to_frame().T], ignore_index=True)
        
        # assign old terraces DF for later comparison and skip out
        OldTerracesDF = TerracesDF
        
        continue
    
    # otherwise loop through the terraces (if any)
    for i, Terrace in TerracesDF.iterrows():
        
        # find the same terrace in the old DF
        OldIndex = (Terrace["MeanElev"] - OldTerracesDF["MeanElev"]).abs().idxmin()
        
        # assign the previous ID to any terraces at similar elevations
        if (np.abs(Terrace["MeanElev"] - OldTerracesDF["MeanElev"].iloc[OldIndex]) < SimilarThresh):
            Terrace["TerraceID"] = OldTerracesDF.iloc[OldIndex]["TerraceID"]
        else:
           # print("New Terrace found, setting ID to list length")
            Terrace["TerraceID"] = len(TerracesTSList)
        
        # set time
        Terrace["Time"] = Time
        
        # if there is not yet a DF for this particular terrace, add a blank one to populate
        if not (Terrace["TerraceID"] <= len(TerracesTSList)-1):
            
            TerracesTSList.append(BlankTerraceDF)
            print("Terrace #, ", len(TerracesTSList)-1, "found")
        
        # append the terrace to the appropriate list
        TerracesTSList[Terrace["TerraceID"]] = pd.concat([TerracesTSList[Terrace["TerraceID"]], Terrace.to_frame().T], ignore_index=True)
    
    # assign old terraces DF for later comparison and skip out
    OldTerracesDF = TerracesDF
    
# loop through terracesDF list and save a spreadsheet for each terrace through time
Filename = AnalysisFolder + str(RunID) + "_Terraces.xlsx"

with pd.ExcelWriter(Filename) as writer:
    
    for i, TerraceDF in enumerate(TerracesTSList):
        
        TerraceDF.to_excel(writer, sheet_name="Terrace_"+str(i), index=False)
        
    
    
