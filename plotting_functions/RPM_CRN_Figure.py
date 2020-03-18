# import plotting tools
from pathlib import Path
import numpy as np
import re
from scipy.stats import mode
from matplotlib import rcParams, cm, gridspec
from matplotlib.lines import Line2D
from cycler import cycler
import matplotlib.pyplot as plt
from RPM_CRN_Plotting_Functions import *

class RPM_CRN_Figure:

    def __init__(self, FigSizeFormat="EPSL", FigWidth_Inches=0., AspectRatio=16./9.):

        """
        This function creates a default matplotlib figure object

        Args:
            FigSizeFormat: the figure size format according to journal for which the figure is intended
                values are geomorphology,ESURF, ESPL, EPSL, JGR, big
                default is ESURF
            
            AspectRatio: The shape of the figure determined by the aspect ratio, default is 16./9.

        Returns:
            matplotlib figure object

        Author: MDH

        """
        self.Figure = None
        self.Axes = None
        
        self.CreateFigure(FigSizeFormat, FigWidth_Inches, AspectRatio)
        
    def CreateFigure(self, FigSizeFormat="EPSL", FigWidth_Inches=0., AspectRatio=16./9.):
    
        """
        This function creates a default matplotlib figure object

        Args:
            FigSizeFormat: the figure size format according to journal for which the figure is intended
                values are geomorphology,ESURF, ESPL, EPSL, JGR, big
                default is ESURF
            
            AspectRatio: The shape of the figure determined by the aspect ratio, default is 16./9.

        Returns:
            matplotlib figure object

        Author: MDH

        """

        # set figure sizes (in inches) based on format
        if FigWidth_Inches > 0:
            FigWidth_Inches = FigWidth_Inches
        elif FigSizeFormat == "geomorphology":
            FigWidth_Inches = 6.25
        elif FigSizeFormat == "big":
            FigWidth_Inches = 16
        elif FigSizeFormat == "small":
            FigWidth_Inches = 3.3
        elif FigSizeFormat == "ESURF":
            FigWidth_Inches = 4.92
        elif FigSizeFormat == "ESPL":
            FigWidth_Inches = 7.08
        elif FigSizeFormat == "EPSL":
            FigWidth_Inches = 7.48
        elif FigSizeFormat == "EPSL_small":
            FigWidth_Inches = 3.74
        elif FigSizeFormat == "JGR":
            FigWidth_Inches = 6.6
        else:
            FigWidth_Inches = 4.92126
            
        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = 10
        rcParams['text.usetex'] = True
        
        # customise the colorcycle for plotting
        rcParams['axes.prop_cycle'] = cycler(color=cm.Dark2.colors)
            
        self.Figure = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio),facecolor=None)

    def PlotProfileAndConcentrationFigure(self, ProfileFile, ConcentrationsFile, Colour=None, Symbol="-", Legend=False, Label=None):

        # if no figure make the default
        if not self.Figure:
            print(self.Figure)
            self.CreateFigure()

        # if axes not created yet add axes as list for subplots and organise labels
        if not self.Axes:
            
            # set up the gridspec
            GridSpec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[2, 1], height_ratios=[1,1])

            # ax0 for profiles, no x axis, y axis on the left
            ax0 = self.Figure.add_subplot(GridSpec[0,0])
            ax0.set_ylabel("Elevation (m)")
            ax0.xaxis.set_visible(False)
            ax0.spines['right'].set_visible(False)
            ax0.spines['top'].set_visible(False)
            ax0.spines['bottom'].set_visible(False)

            # ax1 for concentrations, y axis on the right
            ax1 = self.Figure.add_subplot(GridSpec[1,0])
            ax1.set_yscale("log")
            ax1.set_xlabel("Distance (m)")
            ax1.set_ylabel("Concentration (at g${-1}$)")
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)

            # ax2 axis for time series of retreat rates
            ax2 = self.Figure.add_subplot((GridSpec[0,1]))
            ax2.set_yscale("log")
            ax2.xaxis.set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.yaxis.set_ticks_position('right')
            ax2.yaxis.set_label_position('right')
            ax2.set_ylabel("Cliff Retreat Rate (m yr$^{-1}$)")
            ax2.spines['bottom'].set_visible(False)
            
            # ax3 for time series of maximum concentrations
            ax3 = self.Figure.add_subplot((GridSpec[1,1]))
            ax3.set_yscale("log")
            ax3.set_xlabel("Time (k yrs)")
            ax3.set_ylabel("Max Intertidal Concentration (at g$^{-1}$)")
            ax3.spines['left'].set_visible(False)
            ax3.spines['top'].set_visible(False)
            ax3.yaxis.set_ticks_position('right')
            ax3.yaxis.set_label_position('right')
            

            self.Axes = [ax0, ax1, ax2, ax3]

        # read the profile file
        Times, SeaLevels, Z, X = ReadShoreProfile(ProfileFile)
        LastX = X[-1]

        # find cliff and normalise
        CliffPositions = np.array([mode(EachX[EachX > 1])[0] for EachX in X])
        CliffPositions = np.array([Element for Each in CliffPositions for Element in Each])
        CliffIndices = np.argmin(np.abs(X.T-CliffPositions),axis=0)
        
        LastX -= CliffPositions[-1]
        #self.Axes[0].set_xlim(0, CliffPosition)

        # plot final result on ax0
        Line, = self.Axes[0].plot(LastX, Z, ls=Symbol, color=Colour, label=Label)

        # copy the colour for other plots
        Colour = Line.get_color()
        LineStyles = ['-', '--', ':','-.']
        
        # read the concentrations
        Times2, dX, Concentrations = ReadConcentrationData(ConcentrationsFile)
        
        # populate lines for legend
        LegendLines = []
        LegendLabels = []
        
        for i, key in enumerate(Concentrations.keys()):
            
            N = Concentrations[key][-1]
            XConc = np.arange(0,len(N))*dX
            CliffIndex = np.argmin(np.abs(XConc-CliffPositions[-1]))
            XConc -= XConc[CliffIndex]
            self.Axes[1].plot(XConc[0:CliffIndex], N[0:CliffIndex], color=Colour, ls=LineStyles[i])
            LegendLines.append(Line2D([0], [0], color="grey", ls=LineStyles[i]))
            result = [split for split in re.split('([0-9]+)', key) if split != ""] #lstrip('0123456789')
            Mass = result[0]
            Element = result[1]
            LegendLabels.append("$^{"+Mass+"}$"+Element)
            
            # calculate max concentrations
            MaxN = []
            for Time, CliffPosition, N in zip(Times, CliffPositions, Concentrations[key]):
                XConc = np.arange(0,len(N))*dX
                CliffIndex = np.argmin(np.abs(XConc-CliffPosition))
                MaxN.append(np.max(N[0:CliffIndex]))
            self.Axes[3].plot(Times/1000., MaxN, ls=LineStyles[i], color=Colour)
        
        # calculate cliff retreat rates
        RetreatRates = np.diff(CliffPositions)/np.diff(Times)
        self.Axes[2].plot(Times[1:]/1000,RetreatRates,'-', color=Colour)
        
        # make sure axes line up
        xmin, xmax = self.Axes[0].get_xlim()
        self.Axes[1].set_xlim(xmin, xmax)
        
        # make sure axes line up
        xmin, xmax = self.Axes[2].get_xlim()
        self.Axes[3].set_xlim(xmin, xmax)
        
        # make sure axes line up
        ymin, ymax = self.Axes[1].get_ylim()
        self.Axes[3].set_ylim(ymin, ymax)
        
        # create or update legends
        if Legend:
            self.Axes[0].legend()
            self.Axes[1].legend(LegendLines,LegendLabels)

    def PlotProfileEvolutionFigure(self, ProfileFile, Symbol="-", TimeInterval=1000.):

        """
        """
        # if no figure make the default
        if not self.Figure:
            print(self.Figure)
            self.CreateFigure()

        # if axes not created yet add axes as list for subplots and organise labels
        if not self.Axes:
            
            # ax0 for profiles, no x axis, y axis on the left
            ax0 = self.Figure.add_subplot(111) #GridSpec[0,0])
            ax0.set_ylabel("Elevation (m)")
            #ax0.xaxis.set_visible(False)
            ax0.spines['right'].set_visible(False)
            ax0.spines['top'].set_visible(False)
            #ax0.spines['bottom'].set_visible(False)
            self.Axes = [ax0]

        # read the profile file
        Times, SeaLevels, Z, X = ReadShoreProfile(ProfileFile)
        StartTime = Times[0]
        EndTime = Times[-1]
        Time = StartTime
        OldIndex = -9999

        if StartTime > EndTime:
            TimeInterval *= -1
        
        # set colour map
        ColourMap = cm.bone

        while Time >= 0:
            
            print(Time)
            
            # Find time
            Index = np.argmin(np.abs(Time-Times))
            
            if Index == OldIndex:
                break
            
            OldIndex = Index

            # plot final result on ax0
            Label = str(int(Time)) + " years"
            Colour = ColourMap(Time/np.max([StartTime,EndTime]))
            self.Axes[0].plot(X[Index], Z, ls=Symbol, color=Colour, label=Label)
            Time += TimeInterval
        
        # create or update legends
        self.Axes[0].legend()            

    def SaveFig(self, Outputfilename):
        self.Figure.savefig(Outputfilename)

if __name__ == "__main__":
    Folder = Path(r"C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Results\ModelOutput")
    ProfileFile = Folder / "EnsembleShoreProfile_G1_S_0_T_4_W_0.1_Ws_0.01_R_1_H_2_A_0.1.xz"
    ConcentrationsFile = Folder / "EnsembleConcentrations_G1_S_0_T_4_W_0.1_Ws_0.01_R_1_H_2_A_0.1.xn"
    FigureFile = Folder / "test.png"
    
    MyFigure = RPM_CRN_Figure(FigWidth_Inches=11.)

    MyFigure.PlotProfileAndConcentrationFigure(ProfileFile, ConcentrationsFile, Label="test", Legend=True)
    MyFigure.SaveFig(FigureFile)