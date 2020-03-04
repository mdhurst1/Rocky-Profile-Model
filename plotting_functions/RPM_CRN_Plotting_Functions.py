# import modules
import numpy as np

def ReadShoreProfile(ProfileName):

    """
    Function to read shore profile data from file

    MDH, Feb 2020

    """

    # load the profile file
    f = open(ProfileName,'r')
    Lines = f.readlines()
    f.close()

    StartTime = float(Lines[1].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    MaxZ = Header[0]
    MinZ = Header[1]
    dz = Header[2]
    NTimes = len(Lines)-1
    Z = np.arange(MaxZ,MinZ-dz,-dz)        

    #Get header info and setup X coord
    Times = np.zeros(NTimes)
    
    # loop through lines and append X data to array of vectors
    for i, Line in enumerate(Lines[1:]):
    
        SplitLine = Line.strip().split(" ")
        
        Times[i] = float(SplitLine[0])

        if i == 0:
            X = np.array(SplitLine[2:],dtype="float64")
        else:
            X = np.vstack((X,np.array(SplitLine[2:],dtype="float64")))

    return Times, Z, X

def ReadConcentrationData(self):
        
    """
    Function to read CRN concentration data from file
    
    MDH, Feb 2020
    
    """

    # load the profile file
    f = open(ConcentrationsName,'r')
    Lines = f.readlines()
    f.close()

    # get which nuclides
    Nuclides = Lines[0].strip().split(" ")
    NNuclides = len(Nuclides)
    
    dX = float(Lines[1].strip().split(" ")[0])
    StartTime = float(Lines[2].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])

    # loop through lines and append N data to arrays of vectors
    Lines = Lines[2:]

    Times = []

    for i in range(0, len(Lines), NNuclides):
        
        # get chunk of lines for nuclides
        NuclideLines = Lines[i:i+NNuclides]
        NuclidesBool = False*NNuclides

        Times.append(Lines[i].strip().split(" ")[0])

        # loop through each nuclide line 
        for Line in NuclideLines:
            
            SplitLine = Line.strip().split(" ")
            Nuclide = SplitLine[1]

            if Nuclide == "10":
                if i == 0:
                    N10 = [np.array(SplitLine[2:],dtype="float64"),]
                else:
                    N10.append(np.array(SplitLine[2:],dtype="float64"))
            
            elif Nuclide == "14":
                if i == 0:
                    N14 = [np.array(SplitLine[2:],dtype="float64"),]
                else:
                    N14.append(np.array(SplitLine[2:],dtype="float64"))
            
            elif Nuclide == "26":
                if i == 0:
                    N26 = [np.array(SplitLine[2:],dtype="float64"),]
                else:
                    N26.append(np.array(SplitLine[2:],dtype="float64"))
            
            else:
                sys.exit("Nuclide " + Nuclide + " not recognised!")

    ConcentrationsDict = {}
    if Nuclide == "10":
        ConcentrationsDict["10Be"] = N10
    if Nuclide == "14":
        ConcentrationsDict["14C"] = N14
    if Nuclide == "26":
        ConcentrationsDict["26Al"] = N26
    
    return Times, ConcentrationsDict

    