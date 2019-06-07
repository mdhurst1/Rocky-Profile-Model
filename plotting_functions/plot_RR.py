# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

@author: mhurst
"""

# import plotting tools and set the back end for running on server
#import matplotlib
#matplotlib.use('Agg')

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=12)
rc('ytick.major',pad=1)
rc('xtick.major',pad=1)
padding = 1

def read_profile(FileName):

    # declare the file and the axis
    ProfileName = FileName+"_ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    #StartTime = float(Lines[1].strip().split(" ")[0])
    #EndTime = float(Lines[-1].strip().split(" ")[0])

    # create a place holder for cliff position
    Times = np.zeros(NoLines-1)
    XPos = np.zeros(NoLines-1)
    RSL = np.zeros(NoLines-1)
    
    print(Lines[0])

    for j in range(1,NoLines):
        
        # extract each line
        Line = Lines[j].split(" ")
        
        # record the time to our array
        Times[j-1] = float(Line[0])

        # record the cliff position to array
        XPos[j-1] = float(Line[2])

        #record the RSL to array
        RSL[j-1] = float(Line[1])

    if Times[0] == -9999:
        Times[0] = 10000.

    Times = Times[0::50]
    XPos = XPos[0::50]
    RSL = RSL[0::50]

    return Times, XPos, RSL

CB_Filename = "../../RPM_JRS/Test_4"
SY_Filename = "../../RPM_JRS/MCMC_1"

CB_Times, CB_XPos, CB_RSL = read_profile(CB_Filename)
SY_Times, SY_XPos, SY_RSL = read_profile(SY_Filename)

print CB_XPos

# calculate retreat rates
CB_Rates = np.diff(CB_XPos)/(CB_Times[1]-CB_Times[0])
print(CB_Times[1], CB_Times[0])
SY_Rates = np.diff(SY_XPos)/(SY_Times[1]-SY_Times[0])

#create blank figure
fig = plt.figure(1,figsize=(6.6,5))


#plt.xlim(np.max(CB_Times),np.min(CB_Times))


ax1 = fig.add_subplot(111)
ax1.set_xlabel("Age(BP)")
ax1.set_ylabel("Retreat Rate (m y$^-1$)")
ax1.plot(CB_Times[1:],CB_Rates,'k-', label='Bideford Retreat Rate')
ax1.plot(SY_Times[1:],SY_Rates,'b-', label='Scalby Retreat Rate')

xmin, xmax = ax1.get_xlim()
ax1.set_xlim(7000,0)

ax1.set_xlim(np.max(CB_Times),np.min(CB_Times))
for label in ax1.xaxis.get_ticklabels():
    label.set_rotation(45)

#plot RSL with different axis

ax2 = ax1.twinx()

ax2.set_ylabel("RSL (m)")
ax2.plot(CB_Times[1:],CB_RSL[1:],'k--', label='Bideford RSL')
ax2.plot(SY_Times[1:],SY_RSL[1:],'b--', label='Scalby RSL')
ax2.set_xlim(7000,0)
#ax1.set_ylim(0.01,0)
#ax1.legend(loc='upper left')
#ax2.legend(loc='upper left')

plt.tight_layout(pad=1)

plt.show()
plt.draw()
fig.savefig('poster_RR_CB_SY.png',dpi=300)
