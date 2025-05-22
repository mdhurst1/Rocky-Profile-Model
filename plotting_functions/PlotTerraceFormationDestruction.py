# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 08:54:35 2023

@author: mh322u
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

%matplotlib qt5

Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
AnalysisFolder = Folder + "Analysis/"
PlotsFolder = Folder + "Plots/"
RunID = 2198

# Read terraces timeseries spreadsheet
Filename = AnalysisFolder + str(RunID) + "_Terraces.xlsx"
TerracesSheets = pd.read_excel(Filename,sheet_name=None)

# create figure and axes
fig1 = plt.figure(1,figsize=(16,9))
ax1 = fig1.add_subplot(111)
ax1.set_ylabel("Width (m)")
ax1.set_xlabel("Time (m)")
ax1.set_xlim(13000,0)

# loop through terraces
for i, (sheet_name, TerraceDF) in enumerate(TerracesSheets.items()):

    print(i)
    print(TerraceDF["TerraceID"][0])
    
    # get time and width
    # Calculate y-values for the shaded region
    WidthBottom = (len(TerracesSheets)-TerraceDF["TerraceID"][0])*200 - 0.5*TerraceDF["Width"]
    WidthTop = (len(TerracesSheets)-TerraceDF["TerraceID"][0])*200 + 0.5*TerraceDF["Width"]

    # Plot the shaded region using fill_between
    ax1.fill_between(TerraceDF["Time"], WidthBottom, WidthTop, color=[1., 0.1, 0.1])
    ax1.text(TerraceDF["Time"][0]+800,(len(TerracesSheets)-i)*200,"Terrace " + str(i),va="center")

plt.savefig(PlotsFolder+"WidthTime.png", dpi=300)
plt.show()