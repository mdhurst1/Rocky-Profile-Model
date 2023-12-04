# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 03:20:36 2023

@author: mh322u
"""

# import modules
import pandas as pd
from tqdm import tqdm

#define workspace
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
ResultsFolder = Folder + "Results/"

# locad the records and pick one
RecordFile = "TerraceRuns.xlsx"
RecordDF = pd.read_excel(Folder+RecordFile,sheet_name="Sheet1",header=0)

for i in tqdm(range(len(RecordDF.RunID)),desc="Processing Terraces",unit="%"):
        
    RunID = RecordDF.RunID.iloc[i]
    
    # load uplift history
    UpliftDF = pd.read_csv(Folder+str(RunID)+"_episodic_uplift.data", header=None,columns=["Time","Mag"])
        
    # get RSL at time of uplift
    SeaDF = pd.read_csv(Folder + str(RunID) + "_rsl.data", delimiter=" ", header=0)
    
    NewCol = []
    for Earthquake in UpliftDF:
        # Find the index of the nearest value
        Index = (SeaDF['Time'] - Earthquake["Time"]).abs().idxmin()
        NewCol.append(SeaDF["RSL"].iloc[Index])
    
    UpliftDF["RSL"] = NewCol
    
    UpliftDF.to_csv(Folder+str(RunID)+"_uplift_RSL.data", header=None,columns=["Time","Mag","RSL"])
