from RPM_CRN_Run import *

Folder = "../driver_files/"
Project = "TestProject"
ThisRun = RPM_CRN_Run(Folder,Project)
#ThisRun.PlotProfiles(1000.)
ThisRun.PlotCRNConcentrations()