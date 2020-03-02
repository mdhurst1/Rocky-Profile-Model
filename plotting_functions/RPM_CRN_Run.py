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

        self.Nuclides = None
        self.NNuclides = None
        self.N10 = None
        self.N14 = None
        self.N26 = None


        self.ReadShoreProfile()
        self.ReadConcentrationData()


    def ReadShoreProfile(self):

        """
        Function to read shore profile data from file

        MDH, Feb 2020

        """

        # load the profile file
        f = open(self.ProfileName,'r')
        Lines = f.readlines()
        f.close()

        self.StartTime = float(Lines[1].strip().split(" ")[0])
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


    def ReadConcentrationData(self):
        
        """
        Function to read CRN concentration data from file
        
        MDH, Feb 2020
        
        """

        # load the profile file
        f = open(self.ConcentrationsName,'r')
        Lines = f.readlines()
        f.close()

        # get which nuclides
        self.Nuclides = Lines[0].strip().split(" ")
        self.NNuclides = len(self.Nuclides)
        
        self.dX = float(Lines[1].strip().split(" ")[0])
        self.StartTime = float(Lines[2].strip().split(" ")[0])
        self.EndTime = float(Lines[-1].strip().split(" ")[0])

        #Get header info and setup X coord
        self.Times = np.zeros(self.NTimes)
        
        # loop through lines and append N data to arrays of vectors
        Lines = Lines[2:]

        for i in range(0, len(Lines), self.NNuclides):
            
            print(i)

            # get chunk of lines for nuclides
            NuclideLines = Lines[i:i+self.NNuclides]
            NuclidesBool = False*self.NNuclides

            # loop through each nuclide line 
            for Line in NuclideLines:
                
                SplitLine = Line.strip().split(" ")
                self.Times[i] = float(SplitLine[0])

                Nuclide = SplitLine[1]

                if Nuclide == "10":
                    if i == 0:
                        self.N10 = np.array(SplitLine[2:],dtype="float64")
                    else:
                        self.N10 = np.vstack((self.N10, np.array(SplitLine[2:],dtype="float64")))
                
                elif Nuclide == "14":
                    if i == 0:
                        self.N14 = np.array(SplitLine[2:],dtype="float64")
                    else:
                        self.N14 = np.vstack((self.N14, np.array(SplitLine[2:],dtype="float64")))
                
                elif Nuclide == "26":
                    if i == 0:
                        self.N26 = np.array(SplitLine[2:],dtype="float64")
                    else:
                        self.N26 = np.vstack((self.N14, np.array(SplitLine[2:],dtype="float64")))
                
                else:
                    sys.exit("Nuclide " + Nuclide + " not recognised!")

        
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

    def PlotCRNConcentrations(self):

        """
        Function to plot final concentrations

        """

        #create blank figure
        fig = plt.figure(1,figsize=(6.6,3.3))
        ax1 = plt.subplot(111)
        
        if "10" in self.Nuclides:
            TempX = np.arange(0,len(self.N10[-1],self.dX))
            ax1.plot(TempX,self.N10[-1],'k-',label="$^{10}$Be")
        
        if "14" in self.Nuclides:
            TempX = np.arange(0,len(self.N14[-1],self.dX))
            ax1.plot(TempX,self.N14[-1],'r-',label="$^{14}$C")

        if "26" in self.Nuclides:
            TempX = np.arange(0,len(self.N26[-1],self.dX))
            ax1.plot(TempX,self.N26[-1],'b-',label="$^{26}$C")

        # tweak the plot
        plt.xlabel("Distance (m)")
        plt.ylabel(r"Concentration (atoms g\textsuperscript{-1})")
        plt.xlim(np.min(TempX[-1]),np.max(TempX[-1]))
        
        # add the legend
        ax1.legend(loc='upper right')

        plt.tight_layout()
        
        plt.savefig('RPM_CRN_plotting_test.png',dpi=300)

    def Animation(self):
        """
        Function to plot an animation of the profile and concentrations evolving through time

        MDH, Feb 2020

        """

