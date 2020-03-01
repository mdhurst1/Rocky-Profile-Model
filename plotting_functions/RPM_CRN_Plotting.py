# import modules
import matplotlib.pyplot as plt
import numpy as np

class RPM_CRN_Run:

    """
    Object to contain results of an RPM_CRN model run and analyse and plot

    """


    self...
    
    def __init__(self, Folder, ProjectName):
        """
        Function to inituate RPM_CRN_Run object from output files

        MDH, Feb 2020
        
        """

        self.Folder = Folder
        self.ProjectName = ProjectName

        self.StartTime = None
        self.EndTime = None

        self.MaxZ = None
        self.MinZ = None
        self.dz = None
        self.Z = None


    def ReadShoreProfile(self):

        """
        Function to read shore profile data from file

        MDH, Feb 2020

        """

        # load the profile file
        f = open(ProfileName,'r')
        Lines = f.readlines()
        NoLines = len(Lines)
        StartTime = float(Lines[1].strip().split(" ")[0])
        EndTime = float(Lines[-1].strip().split(" ")[0])

        # Get info on vertical from header
        Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
        MaxZ = Header[0]
        MinZ = Header[1]
        dz = Header[2]
        Line = (Lines[2].strip().split(" "))
        NValues = len(np.array(Line[2:],dtype="float64")) 
        Z = np.arange(MaxZ,MinZ,-dz)
        
    
        #Get header info and setup X coord
        Times = np.zeros(NValues)
        
        for i, Line in enumerate(Lines[2:]):
        
            SplitLine = Line.strip().split(" ")
            
            Times[i] = float(SplitLine[0])
        
            #Read morphology
            X = np.array(SplitLine[2:],dtype="float64")
            NValues = len(X)
            Z = np.linspace(MaxZ,MinZ, NValues)
        
        if (Time == StartTime) or (Time == -9999):
            print(Time)
            ax1.plot(X,Z,'k--',lw=1.,zorder=10,label="Initial Profile")
            PlotTime += PlotInterval
        elif Time <= EndTime:
            ax1.plot(X,Z,'k-',lw=1.,zorder=10,label="Final Profile")
            PlotTime += PlotInterval
            break
        elif (Time <= PlotTime):
            print(Time)
            colour =(Time-StartTime)/(EndTime-StartTime)
            ax1.plot(X,Z,'-',lw=1.,color=ColourMap(colour))
            PlotTime += PlotInterval

    
    ax1.plot(X,Z,'k-',lw=1.5)



    def ReadConcentrationData(self):
        
        """
        Function to read CRN concentration data from file
        
        MDH, Feb 2020
        
        """

    def Save(self, PickleFile):
        """
        Function to save RPM_CRN_Run object
        
        MDH, Feb 2020
        """
        with open(PickleFile, 'wb') as PFile:
            pickle.dump(self, PFile)


    def PlotProfile(self,Time):
        """
        Function to plot a profile line at a specific time

        MDH, Feb 2020

        """

    def PlotProfiles(self,TimeInterval):
        """
        Function to plot profile lines at specific time intervals
        
        MDH, Feb 2020
        
        """

    def PlotCRNConcentrations(self, Time):

    def Animation(self):
        """
        Function to plot an animation of the profile and concentrations evolving through time

        MDH, Feb 2020

        """
