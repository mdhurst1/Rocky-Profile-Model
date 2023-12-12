# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 20:48:03 2023

@author: mh322u
"""

# import modules
import pandas as pd
from Find_Terraces import *
import numpy as np

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
Folder = "C:/Users/mh322u/OneDrive/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"

# locad the records and pick one
RecordFile = "TerraceRuns.xlsx"
RecordDF = pd.read_excel(Folder+RecordFile,sheet_name="Sheet1",header=0)

# set up a table to store the results
Columns = ["Mag", "Freq", "RSL-1", "RSL-2", "RSL-3"]
TotalRunsDF = pd.DataFrame(columns=Columns)
TerracesRunsDF = pd.DataFrame(columns=Columns)
NTerracesDF = pd.DataFrame(columns=Columns)
MeanKappaDF = pd.DataFrame(columns=Columns)

# parameter to extract
# Param = "Width"
    
# loop across uplift magitude whilst keeping all else const
for i, UpliftMag in enumerate(UpliftMags):
    
    # loop across uplift freq
    for j, UpliftFreq in enumerate(UpliftFreqs):
        
        print(i*j, "/", len(UpliftMags)+len(UpliftFreqs), end="\r")
        
        # create records
        TotalRunsRecord = [UpliftMag,UpliftFreq,]
        TerraceRunsRecord = [UpliftMag,UpliftFreq,]
        NTerracesRecord = [UpliftMag,UpliftFreq,]
        KappaRecord = [UpliftMag,UpliftFreq,]
            
        # loop across sea level scenarios
        for k, SeaLevel in enumerate(SeaLevels):
            
            # now loop across everything else and get all values
            # of number of terraces and kappa for table
            Count = 0
            CountTerracesRuns = 0
            KappaList = []
            NTerracesList = []
            
            # clustering
            for Cluster in Clusterings:
            
                # if clustering and frequency are the same there's no model run
                if (UpliftFreq == 4) and (Cluster == 3):
                    continue
            
                for Subs in Subsidences:
                    
                    for Slope in InitSlopes:
                        
                        for Tide in Tides:
                            
                            for Weathering in Weatherings:
                                
                                for Resistance in Resistances:
                                    
                                    for Wave in Waves:
                                        
                                        # get record
                                        Record = RecordDF[   ((RecordDF["SeaLevel"] == SeaLevel) 
                                                           & (RecordDF["UpliftFreq"] == UpliftFreq) 
                                                           & (RecordDF["Clustering"] == Cluster) 
                                                           & (RecordDF["UpliftMag"] == UpliftMag) 
                                                           & (RecordDF["Subsidence"] == Subs) 
                                                           & (RecordDF["InitSlope"] == Slope) 
                                                           & (RecordDF["Tide"] == Tide) 
                                                           & (RecordDF["Weathering"] == Weathering) 
                                                           & (RecordDF["Resistance"] == Resistance) 
                                                           & (RecordDF["Waves"] == Wave))]
                                        
                                        # find scenario
                                        try:
                                            RunID = Record.RunID.values[0]
                                        except:
                                            print("Scenario doesn't exist! Might need to check this!")
                                            print(BaseDict)
                                            print(SearchDict)
                                            sys.exit()
                
                                        # get EQ terraces
                                        TerracesDF = FindTerraces(Folder, RunID)
                                        PlotTerraces(Folder, RunID) #, EQ_Only = True)
                                        
                                        # record if N_Terraces > 1
                                        N_Terraces = len(TerracesDF)
                                        if N_Terraces > 1: 
                                            CountTerracesRuns += 1
                                        
                                        # calucalte kappa
                                        N_Uplift = len(pd.read_csv(Folder+"Results/"+str(RunID)+"_episodic_uplift.data", header=None))
                                        Kappa = np.round((N_Terraces-1)/N_Uplift,2)
                                        KappaList.append(Kappa)
                                        Count += 1
            
            # Update Records and add to DFs
            TotalRunsRecord.append(Count)
            TerraceRunsRecord.append(CountTerracesRuns)
            NTerracesRecord.append(np.mean(N_Terraces))
            KappaRecord.append(np.mean(KappaList))
            
        # add to DFs
        TotalRunsDF = pd.concat([TotalRunsDF,pd.DataFrame([TotalRunsRecord], columns=Columns)], ignore_index=True)
        TerracesRunsDF = pd.concat([TerracesRunsDF,pd.DataFrame([TerraceRunsRecord], columns=Columns)], ignore_index=True)
        NTerracesDF = pd.concat([NTerracesDF,pd.DataFrame([NTerracesRecord], columns=Columns)], ignore_index=True)
        MeanKappaDF = pd.concat([MeanKappaDF,pd.DataFrame([KappaRecord], columns=Columns)], ignore_index=True)
            
# write to file
TotalRunsDF.to_excel(Folder+"Uplift_RSL_TotalRuns.xlsx")
TerracesRunsDF.to_excel(Folder+"Uplift_RSL_Terraces.xlsx")
NTerracesDF.to_excel(Folder+"Uplift_RSL_MeanNTerraces.xlsx")
MeanKappaDF.to_excel(Folder+"Uplift_RSL_MeanKappa.xlsx")

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


