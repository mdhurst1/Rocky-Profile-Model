# import plotting tools
from pathlib import Path
import numpy as np
from matplotlib import rcParams
from matplotlib.lines import Line2D
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
            
        self.Figure = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio),facecolor=None)

    def PlotProfileAndConcentrationFigure(self, ProfileFile, ConcentrationsFile, Colour="k", Symbol="-", Legend=False, Label=None):

        # if no figure make the default
        if not self.Figure:
            print(self.Figure)
            self.CreateFigure()

        # if axes not created yet add axes as list for subplots and organise labels
        if not self.Axes:
            
            # ax1 for profiles, no x axis, y axis on the left
            ax1 = self.Figure.add_subplot(211)
            ax1.set_ylabel("Elevation (m)")
            ax1.xaxis.set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)

            # ax2 for concentrations, y axis on the right
            ax2 = self.Figure.add_subplot(212)
            ax2.set_yscale("log")
            ax2.set_xlabel("Distance (m)")
            ax2.set_ylabel("Concentration (at g${-1}$)")
            ax2.yaxis.set_ticks_position('right')
            ax2.yaxis.set_label_position('right')
            ax2.spines['left'].set_visible(False)
            ax2.spines['top'].set_visible(False)

            self.Axes = [ax1, ax2]

        # read the profile file
        Times, Z, X = ReadShoreProfile(ProfileFile)
        
        # plot final result on ax1
        self.Axes[0].plot(X[-1], Z, ls=Symbol, color=Colour, label=Label)

        # read the concentrations
        Times2, dX, Concentrations = ReadConcentrationData(ConcentrationsFile)
        
        # populate lines for legend
        LegendLines = []
        LegendLabels = []
        
        print(Concentrations.keys())

        if "26Al" in Concentrations.keys():
            print("26Al")
            N26Al = Concentrations["26Al"][-1]
            X = np.arange(0,len(N26Al))*dX
            self.Axes[1].plot(X, N26Al, ":", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls=":"))
            LegendLabels.append("$^{26}$Al")
        
        if "14C" in Concentrations.keys():
            print("14C")
            N14C = Concentrations["14C"][-1]
            X = np.arange(0,len(N14C))*dX
            self.Axes[1].plot(X, N14C, "--", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls="--"))
            LegendLabels.append("$^{14}$C")
            
        if "10Be" in Concentrations.keys():
            print("10Be")
            N10Be = Concentrations["10Be"][-1]
            X = np.arange(0,len(N10Be))*dX
            self.Axes[1].plot(X, N10Be, "-", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls="-"))
            LegendLabels.append("$^{10}$Be")

        # make sure axes line up
        xmin,xmax = self.Axes[0].get_xlim()
        self.Axes[1].set_xlim(xmin, xmax)
        
        # create or update legends
        if Legend:
            self.Axes[0].legend()
            self.Axes[1].legend(LegendLines,LegendLabels)

    def SaveFig(self, Outputfilename):
        self.Figure.savefig(Outputfilename)


if __name__ == "__main__":
    Folder = Path(r"C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Results")
    ProfileFile = Folder / "EnsembleShoreProfile_G0_S_0_T_4_W_0.05_Ws_0.01_R_1_H_1_A_0.1.xz"
    ConcentrationsFile = Folder / "EnsembleConcentrations_G0_S_0_T_4_W_0.05_Ws_0.01_R_1_H_1_A_0.1.xn"
    FigureFile = Folder / "test.png"
    
    MyFigure = RPM_CRN_Figure(FigWidth_Inches=11.)
    MyFigure.PlotProfileAndConcentrationFigure(ProfileFile, ConcentrationsFile, Label="test", Legend=True)
    MyFigure.SaveFig(FigureFile)