# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 05:51:51 2023

@author: mh322u
"""

# import modules
from Find_Terraces import *
import pandas as pd
import sys

# define workspace
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
ResultsFolder = Folder + "Results/"

# load the records
RecordFile = "TerraceRuns.xlsx"
RecordDF = pd.read_excel(Folder+RecordFile,sheet_name="Sheet1",header=0)
Columns = RecordDF.columns

# get a list of unique values for each param from the Record File
NoID = RecordDF.iloc[:,2:]
UniqueValuesDF = pd.DataFrame({col: pd.Series(pd.unique(NoID[col])) for col in NoID.columns})

# get a list of parameters for testing, and remove SL
ParamsList = list(Columns)[1:]
ParamsList.remove("SeaLevel")
ParamsList.remove("RunID")

# build a table of delta kappa as a function of RSL
KappaDF = pd.DataFrame(columns=["Scenario","SL1","SL2","SL3"])
DeltaKappaDF = pd.DataFrame(columns=ParamsList)

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

# loop through RSL scenarios
for SeaLevel in UniqueValuesDF["SeaLevel"]:
    
    print("Sea Level Scenario:", SeaLevel)
    NewRecord = []
        
    # loop through each parameter in turn
    for Param in ParamsList:
        
        print("\tParameter:", Param)
        
        Kappas = []
        
        # get range of unique values
        Values = UniqueValuesDF[Param]
        
        # loop across values
        for Value in Values:
            
            if np.isnan(Value):
                continue
            
            # get record
            SearchDict = BaseDict.copy()
            SearchDict["SeaLevel"] = SeaLevel
            SearchDict[Param] = Value
            
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
            
            # calculate kappa
            N_Terraces = len(TerracesDF)
            N_Uplift = len(pd.read_csv(Folder+"Results/"+str(RunID)+"_episodic_uplift.data", header=None))
            Kappa = np.round((N_Terraces-1)/N_Uplift,2)
            Kappas.append(Kappa)
        
        # calculate delta kappa
        DeltaKappa = np.max(Kappas) - np.min(Kappas)
        
        # add value to record
        NewRecord.append(DeltaKappa)
        
    #add record to table
    DeltaKappaDF = pd.concat([DeltaKappaDF,pd.DataFrame([NewRecord], columns=ParamsList)], ignore_index=True)

DeltaKappaDF.transpose().to_excel(Folder+"RSL_Uplift_DeltaKappa.xlsx")

