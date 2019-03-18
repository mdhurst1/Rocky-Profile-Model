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
rc('font',size=8)
rc('ytick.major',pad=1)
rc('xtick.major',pad=1)
padding = 1

def make_plot(FileName,ColourMap):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,3.3))

    # declare the file and the axis
    ProfileName = FileName+"_ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    #StartTime = float(Lines[1].strip().split(" ")[0])
    #EndTime = float(Lines[-1].strip().split(" ")[0])

    # create a place holder for cliff position
    Times = np.zeros(NoLines-1)
    CliffPositionX = np.zeros(NoLines-1)
    
    for j in range(1,NoLines):
        
        # extract each line
        Line = Lines[j].split(" ")
        
        # record the time to our array
        Times[j-1] = float(Line[0])

        # record the cliff position to array
        CliffPositionX[j-1] = float(Line[2])

    Times = Times[1::10]
    CliffPositionX = CliffPositionX[1::10]

    # calculate retreat rates
    Rates = np.diff(CliffPositionX)/(Times[1]-Times[0])
    plt.plot(Times[1:],Rates,'k-')
    plt.xlim(np.max(Times),np.min(Times))

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('plot_RR_4.png',dpi=300)

if __name__ == "__main__":
    FileName = "../../RPM_JRS/Test_SLR_4" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        
