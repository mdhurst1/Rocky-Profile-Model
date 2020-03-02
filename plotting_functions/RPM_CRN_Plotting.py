from RPM_CRN_Run import *

Folder = "../../RPM_CRN_test/"
Project = "RPM_CRN"
ThisRun = RPM_CRN_Run(Folder,Project)
#ThisRun.PlotProfiles(1000.)
ThisRun.PlotCRNConcentrations()