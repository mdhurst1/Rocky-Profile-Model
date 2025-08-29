from RPM_CRN_Run import *


Folder = '/media/14TB_RAID_Array/User_Homes/Martin_Hurst/Rocky-Profile-Model/test/'
Project = "TestProject"

ThisRun = RPM_CRN_Run(Folder,Project)
#ThisRun.PlotProfiles(100)
ThisRun.PlotCRNConcentrations()