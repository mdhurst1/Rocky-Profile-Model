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
    RSL = np.zeros(NoLines-1)
    
    for j in range(1,NoLines):
        
        # extract each line
        Line = Lines[j].split(" ")
        
        # record the time to our array
        Times[j-1] = float(Line[0])

        # record the cliff position to array
        CliffPositionX[j-1] = float(Line[2])

        #record the RSL to array
        RSL[j-1] = float(Line[1])

    Times = Times[1::50]
    CliffPositionX = CliffPositionX[1::50]
    RSL = RSL[1::50]

    # calculate retreat rates
    Rates = np.diff(CliffPositionX)/(Times[1]-Times[0])
    plt.plot(Times[1:],Rates,'k-')
    plt.xlim(np.max(Times),np.min(Times))


    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Age(BP)")
    ax1.set_ylabel("Retreat Rate (m y$^-1$)")
    ax1.plot(Times[1:],Rates,'k-', label='Retreat Rate')
    ax1.set_xlim(np.max(Times),np.min(Times))
    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)

    #plot RSL with different axis

    ax2 = ax1.twinx()

    ax2.set_ylabel("RSL (m)")
    ax2.plot(Times[1:],RSL[1:],'b-', label='RSL')


    ax1.legend(loc='lower right')
    ax2.legend(loc='upper left')

    plt.tight_layout(pad=1)
    

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('plot_RR_RSL_DB_6.png',dpi=300)


if __name__ == "__main__":
    FileName = "../../RPM_JRS/DB_Test_6" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        
