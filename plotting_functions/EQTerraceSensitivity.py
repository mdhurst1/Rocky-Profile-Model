# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 05:51:51 2023

@author: mh322u
"""

# import modules
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

# build a table of delta kappa as a function of RSL
DeltaKappaDF = pd.DataFrame(columns=ParamsList)

# Define base scenario

# loop through RSL scenarios
for SeaLevel in SeaLevels:
    
    Record = []
    
    # loop through each parameter in turn
    for Param in ParamList:
        
        # get range of unique values
        Values = UniqueValues[Param]
        
        # loop across values
        for Value in Values:
            
            # find scenario
            
            # get EQ terraces
            
            # calculate kappa
            
        # calculate delta kappa
        
        # add value to record
        Record.append(DeltaKappa)
        
    #add record to table
    