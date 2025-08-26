# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:36:34 2020

@author: mdhurst
"""

from pathlib import Path
from RPM_CRN_Figure import *

Folder = Path('/media/14TB_RAID_Array/User_Homes/Martin_Hurst/Rocky-Profile-Model/test')
ProfileFile = Folder / "TestProject_ShoreProfile.xz"

FigureFile = Folder / "test.png"
MyFigure = RPM_CRN_Figure(FigWidth_Inches=11.)
MyFigure.PlotProfileEvolutionFigure(ProfileFile)
MyFigure.SaveFig(FigureFile)