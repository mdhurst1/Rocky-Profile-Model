# import modules
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=1)
rc('xtick.major',pad=1)
padding = 1

class RPM_CRN_Run:

    """
    Object to contain results of an RPM_CRN model run and analyse and plot

    """

    def __init__(self, Folder, ProjectName):
        """
        Function to inituate RPM_CRN_Run object from output files

        MDH, Feb 2020
        
        """

        self.Folder = Folder
        self.ProjectName = ProjectName
        self.ProfileName = Folder + ProjectName + "_ShoreProfile.xz"
        self.ConcentrationsName = Folder + ProjectName + "_Concentrations.xn"

        self.StartTime = None
        self.EndTime = None
        self.NTimes = None
        self.Times = None

        self.MaxZ = None
        self.MinZ = None
        self.dz = None
        self.Z = None
        self.X = None

        self.ReadShoreProfile()


    def ReadShoreProfile(self):

        """
        Function to read shore profile data from file

        MDH, Feb 2020

        """

        # load the profile file
        f = open(self.ProfileName,'r')
        Lines = f.readlines()

        self.StartTime = float(Lines[1].strip().split(" ")[0])
        i=2
        while self.StartTime == -9999:
            self.StartTime = float(Lines[i].strip().split(" ")[0])
            i += 1

        self.EndTime = float(Lines[-1].strip().split(" ")[0])

        # Get info on vertical from header
        Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
        self.MaxZ = Header[0]
        self.MinZ = Header[1]
        self.dz = Header[2]
        Line = (Lines[2].strip().split(" "))
        self.NTimes = len(np.array(Line[2:],dtype="float64")) 
        self.Z = np.arange(self.MaxZ,self.MinZ-self.dz,-self.dz)        
    
        #Get header info and setup X coord
        self.Times = np.zeros(self.NTimes)
        
        # loop through lines and append X data to array of vectors
        for i, Line in enumerate(Lines[2:]):
        
            SplitLine = Line.strip().split(" ")
            
            self.Times[i] = float(SplitLine[0])

            if i == 0:
                self.X = np.array(SplitLine[2:],dtype="float64")
            else:
                self.X = np.vstack((self.X,np.array(SplitLine[2:],dtype="float64")))

        # clean up
        f.close()

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

        #create blank figure
        fig = plt.figure(1,figsize=(6.6,3.3))
        ax1 = plt.subplot(111)
        plt.axis('equal')

        # Only plot every so many years
        PlotTime = self.StartTime
        
        #Colourmap
        ColourMap = cm.bone_r
        
        #Loop through times and plot at time interval
        while PlotTime >= self.EndTime:
            
            print(PlotTime)

            # get index
            Index = np.argmin(np.abs(self.Times-PlotTime))
            
            if (PlotTime == self.StartTime):
                ax1.plot(self.X[Index], self.Z, 'k--', lw=1., zorder=10, label="Initial Profile")
                
            elif PlotTime == self.EndTime:
                ax1.plot(self.X[-1], self.Z, 'k-',lw=1., zorder=10, label="Final Profile")
                break
            
            else:
                colour = (PlotTime-self.StartTime)/(self.EndTime-self.StartTime)
                ax1.plot(self.X[Index], self.Z, '-', lw=1., color=ColourMap(colour))
            
            PlotTime -= TimeInterval

        # tweak the plot
        plt.xlabel("Distance (m)")
        plt.ylabel("Elevation (m)")
        plt.xlim(np.min(self.X[-1]),np.max(self.X[-1]))
        plt.ylim(self.MinZ,self.MaxZ)
        #plt.ylim(-30,30)

        #add the colourbar
        #cbar = fig.colorbar(self.Times)
        #cbar.set_label('Time (years)')

        # add the legend
        plt.tight_layout()
        
        plt.savefig('RPM_CRN_plotting_test.png',dpi=300)

    def PlotCRNConcentrations(self, Time):

        """
        Function to plot concentrations evolving through time goes here
        """

    def Animation(self):
        """
        Function to plot an animation of the profile and concentrations evolving through time

        MDH, Feb 2020

        """

