# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 20:48:03 2023

@author: mh322u
"""

# import modules
import pandas as pd
from Find_Terraces import *

# range of values used in simulations
SeaLevels = [1,2,3]
UpliftFreqs = [2,3,4]
Clusterings = [1,2,3]
UpliftMags = [1,2,3]
Subsidences = [1,2,3]
InitSlopes = [1,2,3]
Tides = [1,2,3]
Weatherings = [1,2]
Resistances = [1,2]
Waves = [1,2]

# setup workspace
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"

# locad the records and pick one
RecordFile = "TerraceRuns.xlsx"
RecordDF = pd.read_excel(Folder+RecordFile,sheet_name="Sheet1",header=0)

# set up a table to store the results
Columns = ["Uplift", "RSL-1", "RSL-2", "RSL-3"]
NTerracesDF = pd.DataFrame(columns=Columns)
KappaDF = pd.DataFrame(columns=Columns)
RunDF = pd.DataFrame(columns=Columns)

# parameter to extract
# Param = "Width"

# Define base scenario
BaseDict = { "SeaLevel": 1,
             "UpliftFreq": 3,   # medium (every 1ky)
             "UpliftMag": 2,    # medium (2m)
             "Clustering": 1,   # no clustering
             "Subsidence": 1,   # no subsidence
             "InitSlope": 2,    # 35 degrees
             "Tide": 2,         # medium (2 m range)
             "Weathering": 1,   # low efficacy (0.001*Resistance)
             "Resistance": 1,   # low resistance (10 kg/m2/s)
             "Waves": 2 }       # high efficacy, low gamma 0.01
    
# loop across uplift magitude whilst keeping all else const
for i, UpliftFreq in enumerate(UpliftFreqs):
    
    # create a blank record
    NTRow = ["UpFreq_" + str(i)]
    KappaRow = ["UpFreq_" + str(i)]
    RunRow = ["UpFreq_" + str(i)]
    
    for SeaLevel in SeaLevels:
        
        # get record
        SearchDict = BaseDict.copy()
        SearchDict["SeaLevel"] = SeaLevel
        SearchDict["UpliftFreq"] = UpliftFreq
        
        # Find the record
        Record = RecordDF[((RecordDF["SeaLevel"] == SearchDict["SeaLevel"]) 
                           & (RecordDF["UpliftFreq"] == SearchDict["UpliftFreq"])
                           & (RecordDF["UpliftMag"] == SearchDict["UpliftMag"])
                           & (RecordDF["Clustering"] == SearchDict["Clustering"])
                           & (RecordDF["Subsidence"] == SearchDict["Subsidence"])
                           & (RecordDF["InitSlope"] == SearchDict["InitSlope"]) 
                           & (RecordDF["Tide"] == SearchDict["Tide"]) 
                           & (RecordDF["Weathering"] == SearchDict["Weathering"]) 
                           & (RecordDF["Resistance"] == SearchDict["Resistance"]) 
                           & (RecordDF["Waves"] == SearchDict["Waves"]))]
        
        RunID = Record.RunID.values[0]
        TerracesDF = FindTerraces(Folder, RunID)
        PlotTerraces(Folder, RunID)
        
        N_Terraces = len(TerracesDF)
        N_Uplift = len(pd.read_csv(Folder+"Results/"+str(RunID)+"_episodic_uplift.data", header=None))
        
        NTRow.append(N_Terraces)
        Kappa = np.round((N_Terraces-1)/N_Uplift,2)
        print(RunID,N_Terraces,N_Uplift,Kappa)
        KappaRow.append(Kappa)
        RunRow.append(RunID)
        
    NTerracesDF = pd.concat([NTerracesDF,pd.DataFrame([NTRow], columns=Columns)], ignore_index=True)
    KappaDF = pd.concat([KappaDF,pd.DataFrame([KappaRow], columns=Columns)], ignore_index=True)
    RunDF = pd.concat([RunDF,pd.DataFrame([RunRow], columns=Columns)], ignore_index=True)
        
# loop across uplift freq
for i, UpliftMag in enumerate(UpliftMags):
    
    # create a blank record
    NTRow = ["UpMag_" + str(i)]
    KappaRow = ["UpMag_" + str(i)]
    RunRow = ["UpMag_" + str(i)]
    
    for i, SeaLevel in enumerate(SeaLevels):
        
        # get record
        SearchDict = BaseDict.copy()
        SearchDict["SeaLevel"] = SeaLevel
        SearchDict["UpliftMag"] = UpliftMag
        
        # Find the record
        Record = RecordDF[((RecordDF["SeaLevel"] == SearchDict["SeaLevel"]) 
                           & (RecordDF["UpliftFreq"] == SearchDict["UpliftFreq"])
                           & (RecordDF["UpliftMag"] == SearchDict["UpliftMag"])
                           & (RecordDF["Clustering"] == SearchDict["Clustering"])
                           & (RecordDF["Subsidence"] == SearchDict["Subsidence"])
                           & (RecordDF["InitSlope"] == SearchDict["InitSlope"]) 
                           & (RecordDF["Tide"] == SearchDict["Tide"]) 
                           & (RecordDF["Weathering"] == SearchDict["Weathering"]) 
                           & (RecordDF["Resistance"] == SearchDict["Resistance"]) 
                           & (RecordDF["Waves"] == SearchDict["Waves"]))]
        
        RunID = Record.RunID.values[0]
        TerracesDF = FindTerraces(Folder, RunID)
        PlotTerraces(Folder, RunID)
        
        N_Terraces = len(TerracesDF)
        N_Uplift = len(pd.read_csv(Folder+"Results/"+str(RunID)+"_episodic_uplift.data", header=None))
        
        NTRow.append(N_Terraces)
        Kappa = np.round((N_Terraces-1)/N_Uplift,2)
        print(RunID,N_Terraces,N_Uplift,Kappa)
        KappaRow.append(Kappa)
        RunRow.append(RunID)
    
    NTerracesDF = pd.concat([NTerracesDF,pd.DataFrame([NTRow], columns=Columns)], ignore_index=True)
    KappaDF = pd.concat([KappaDF,pd.DataFrame([KappaRow], columns=Columns)], ignore_index=True)
    RunDF = pd.concat([RunDF,pd.DataFrame([RunRow], columns=Columns)], ignore_index=True)
    
# loop across clustering
for i, Clustering in enumerate(Clusterings):
    
    # create a blank record
    NTRow = ["Cluster_" + str(i)]
    KappaRow = ["Cluster_" + str(i)]
    RunRow = ["Cluster_" + str(i)]
    
    for i, SeaLevel in enumerate(SeaLevels):
        
        # get record
        SearchDict = BaseDict.copy()
        SearchDict["SeaLevel"] = SeaLevel
        SearchDict["Clustering"] = Clustering
        
        # Find the record
        Record = RecordDF[((RecordDF["SeaLevel"] == SearchDict["SeaLevel"]) 
                           & (RecordDF["UpliftFreq"] == SearchDict["UpliftFreq"])
                           & (RecordDF["UpliftMag"] == SearchDict["UpliftMag"])
                           & (RecordDF["Clustering"] == SearchDict["Clustering"])
                           & (RecordDF["Subsidence"] == SearchDict["Subsidence"])
                           & (RecordDF["InitSlope"] == SearchDict["InitSlope"]) 
                           & (RecordDF["Tide"] == SearchDict["Tide"]) 
                           & (RecordDF["Weathering"] == SearchDict["Weathering"]) 
                           & (RecordDF["Resistance"] == SearchDict["Resistance"]) 
                           & (RecordDF["Waves"] == SearchDict["Waves"]))]
    
        RunID = Record.RunID.values[0]
        TerracesDF = FindTerraces(Folder,RunID)
        PlotTerraces(Folder, RunID)
        
        N_Terraces = len(TerracesDF)
        N_Uplift = len(pd.read_csv(Folder+"Results/"+str(RunID)+"_episodic_uplift.data", header=None))
        
        NTRow.append(N_Terraces)
        Kappa = np.round((N_Terraces-1)/N_Uplift,2)
        print(RunID,N_Terraces,N_Uplift,Kappa)
        KappaRow.append(Kappa)
        RunRow.append(RunID)
                
    NTerracesDF = pd.concat([NTerracesDF,pd.DataFrame([NTRow], columns=Columns)], ignore_index=True)
    KappaDF = pd.concat([KappaDF,pd.DataFrame([KappaRow], columns=Columns)], ignore_index=True)
    RunDF = pd.concat([RunDF,pd.DataFrame([RunRow], columns=Columns)], ignore_index=True)
    
NTerracesDF.to_excel(Folder+"RSL_Uplift_NTerraces.xlsx")
KappaDF.to_excel(Folder+"RSL_Uplift_Kappa.xlsx")
RunDF.to_excel(Folder+"RSL_Uplift_Run.xlsx")

# Create a figure and axis
#fig, ax = plt.subplots(figsize=(6, 9))

# Plot the table
#table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')

# Color code cells based on values
#for i in range(1, len(df.columns) + 1):  # Skip the first column (index 0)
#    col_values = df.iloc[:, i - 1]
#    cell_colors = plt.cm.RdYlBu(col_values / col_values.max())  # Adjust the colormap as needed
#    for j in range(len(df)):
#        table[(j + 1, i)].set_facecolor(cell_colors[j])

# Hide the axes
#ax.axis('off')

# Show the plot
#plt.show()


