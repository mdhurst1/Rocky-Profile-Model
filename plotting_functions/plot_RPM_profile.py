# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:36:34 2020

@author: mdhurst
"""

from pathlib import Path
from RPM_CRN_Figure import *

Folder = Path(r"C:\Users\Martin Hurst\OneDrive - University of Glasgow\Projects\RockCoastCosmo\CoupledModelling\Rocky-Profile-Model\driver_files")
ProfileFile = Folder / "RPM_CRN_ShoreProfile.xz"

FigureFile = Folder / "test.png"
MyFigure = RPM_CRN_Figure(FigWidth_Inches=11.)
MyFigure.PlotProfileEvolutionFigure(ProfileFile)
MyFigure.SaveFig(FigureFile)