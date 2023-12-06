# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 01:09:55 2023

@author: mh322u
"""

def TrackTerraces(Folder, RunID, MinWidth=3.):
    
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
    SimilarThresh = 0;
    
    TerraceTSList = []
    BlankTerraceDF = pd.DataFrame(columns = ["Time","TerraceID","StartIndex","EndIndex","Width","MeanElev","ElevChange","Slope"])
    
    # loop through time
    for i, Time in enumerate(Times):
    
        # identify a terrace and track its progress
        TerracesDF = FindTerracesRaw(X[i],Z)
        
        # first time out no old terraces to compare to
        # or if no terraces have previously been found
        if (i==0) or (len(TerracesTSList) == 0):
            # loop through terraces, probs only one
            for Terrace in TerracesDF:
                
                # add time to DF
                Terrace["Time"] = Time
                
                # create the new DF for this new terrace and add the terrace to it
                TerraceTRList.append(BlankTerraceDF)
                TerraceTSList[Terrace["TerraceID"-1]].append(Terrace)
            
            # assign old terraces DF for later comparison and skip out
            OldTerracesDF = TerracesDF
            continue
        
        # otherwise loop through the terraces (if any)
        for i, Terrace in TerracesDF.iterrows():
            
            # find the same terrace in the old DF
            OldIndex = (Terrace["MeanElev"] - OldTerracesDF).abs().idxmin()
            
            # assign the previous ID to any terraces at similar elevations
            if ((Terrace["MeanElev"] - OldTerracesDF.iloc[OldIndex]).abs() < SimilarThresh):
                Terrace["TerraceID"] = OldTerracesDF.iloc[OldIndex]["TerraceID"]
                Terrace["Time"] = Time
            
            # if there is not yet a DF for this particular terrace, add a blank one to populate
            if not TerraceTSList[Terrace["TerraceID"]]:
                TerraceTRList.append(BlankTerraceDF)
            
            # append the terrace to the appropriate list
            TerraceTSList[Terrace["TerraceID"].append(Terrace)

    # loop through terracesDF list and save a spreadsheet for each terrace through time
    for i, TerracesDF in enumerate(TerracesTSList):
        TerracesDF.to_excel(AnalysisFolder + str(RunID)+"_Terrace_" + str(TerracesDF["TerraceID"] + ".xlsx")
        
if __name__ == "__main__":
    
    # setup workspace
    Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Results/"
    
    # RunID = Record.RunID.values[0]
    RunID = 2198
    TrackTerraces(Folder, RunID)