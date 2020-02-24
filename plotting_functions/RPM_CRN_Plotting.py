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


    def ReadShoreProfile(self):

        """
        Function to read shore profile data from file

        MDH, Feb 2020

        """    

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
