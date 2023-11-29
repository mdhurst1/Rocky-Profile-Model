#import modules
from RPM_CRN_Analysis_Functions import *
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm

#load an example profile
Folder = "C:/Users/mh322u/OneDrive - University of Glasgow/NZ_2023/Modelling_Holcene_Terraces/Uplift_RPM/"
ProfileFileName = "Uplift_ShoreProfile.xz"
Times, SeaLevels, Z, X = ReadShoreProfile(Folder+ProfileFileName)

# figure properties
# Set up fonts for plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 10

# plot the profile
fig = plt.figure(1,figsize=(16,9))
ax = fig.add_subplot(111)
ax.set_ylabel("Elevation (m)")
ax.set_xlabel("Distance (m)")

StartTime = Times[0]
EndTime = Times[-1]
Time = StartTime
TimeInterval = 1000
OldIndex = -9999

if StartTime > EndTime:
    TimeInterval *= -1

# set colour map
ColourMap = cm.bone

while Time >= EndTime:
    
    # Find time
    Index = np.argmin(np.abs(Time-Times))
    
    if Index == OldIndex:
        break
    
    OldIndex = Index

    # plot final result on ax0
    Label = str(int(Time)) + " years"
    Colour = ColourMap(Time/np.max([StartTime,EndTime]))
    ax.plot(X[Index], Z, ls="-", color=Colour, label=Label)

    Time += TimeInterval

# create or update legends
fig.savefig("Folder+ProfilePlot.png")