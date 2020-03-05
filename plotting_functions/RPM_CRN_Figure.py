# import plotting tools
from pathlib import Path
import matplotlib
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt

class RPM_CRN_Figure:

    def __init__(FigSizeFormat="EPSL", FigWidth_Inches=0., AspectRatio=16./9.):

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
        self.Legends = None

        CreateFigure(FigSizeFormat, FigWidth_Inches, AspectRatio)
        
    def CreateFigure(FigSizeFormat="EPSL", FigWidth_Inches=0., AspectRatio=16./9.):
    
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
            
        self.Fig = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio),facecolor=None)

    def PlotProfileAndConcentrationFigure(self, ProfileFile, ConcentrationsFile, Colour="k", Symbol="-", Legend=False, Label=None):

        # if no figure make the default
        if not Figure:
            self.CreateFigure()

        # if axes not created yet add axes as list for subplots and organise labels
        if not Axes:
            
            # ax1 for profiles, no x axis, y axis on the left
            ax1 = Figure.add_subplot(211)
            ax1.set_ylabel("Elevation (m)")
            ax1.xaxis.set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)

            # ax2 for concentrations, y axis on the right
            ax2 = Figure.add_subplot(212)
            ax2.set_yscale("log")
            ax2.set_xlabel("Distance (m)")
            ax2.set_ylabel("Concentration (at g${-1}$)")
            ax2.yaxis.set_ticks_position('right')
            ax2.yaxis.set_label_position('right')

            self.Axes = [ax1, ax2]

        # read the profile file
        Times, Z, X = ReadShoreProfile(ProfileFile)

        # plot final result on ax1
        ax1.plot(X[-1], Z, ls=Symbol, color=Colour, label=Label)

        # read the concentrations
        Times2, dX, Concentrations = ReadConcentrationData(ConcentrationsName)
        
        # populate lines for legend
        LegendLines = []
        LegendLabels = []

        if "10Be" in Concentrations.keys():
            N10Be = Concentrations["10Be"][-1]
            X = np.arange(0,len(N10Be))*dX
            ax2.plot(X, N10Be, "-", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls="-"))
            LegendLabels.append("$^{10}$Be")

        if "14C" in Concentrations.keys():
            N14C = Concentrations["14C"][-1]
            X = np.arange(0,len(N14C))*dX
            ax2.plot(X, N14C, "--", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls="--"))
            LegendLabels.append("$^{14}$C")

        if "26Al" in Concentrations.keys():
            N26Al = Concentrations["26Al"][-1]
            X = np.arange(0,len(N26Al))*dX
            ax2.plot(X, N26Al, ":", color=Colour)
            LegendLines.append(Line2D([0], [0], color="grey", ls=":"))
            LegendLabels.append("$^{26}$Al")

        # create or update legends
        if Legend:
            ax1.legend()
            ax2.legend(LegendLines,LegendLabels)

    def SaveFig(self, Outputfilename):
        Figure.savefig(Outputfilename)


if __name__ == "__main__":
    Folder = Path("C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Results")
    ProfileFile = Folder / "EnsembleConcentrations_G0.1_S_0_T_4_W_0.05_Ws_0.01_R_100_H_1_A_0.1.xz"
    ConcentrationsFile = Folder / "EnsembleConcentrations_G0.1_S_0_T_4_W_0.05_Ws_0.01_R_100_H_1_A_0.1.xn"
    FigureFile = Folder / "test.png"
    
    RPM_CRN_Figure(FigWidth_Inches=11.)
    RPM_CRN_Figure.PlotProfileAndConcentrationFigure(ProfileFile, ConcentrationsFile, Label="test", Legend=True)
    RPM_CRN_Figure.SaveFig(FigureFile)
)