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

    # Get Header Values
    #HeaderLine = MorphLines[0].strip().split(" ")
    #CliffHeight = float(HeaderLine[0])
    #MinElev = float(HeaderLine[1])
    #dZ = float(HeaderLine[2])
    #NValues = (int)((CliffHeight-MinElev)/dZ+1)
    #Z = np.linspace(CliffHeight, MinElev, NValues)

    # create a place holder for RSL
    Times = np.zeros(NoLines-1)
    RSL = np.zeros(NoLines-1)
    
    for j in range(2,NoLines):
        
        # extract each line
        Line = Lines[j].split(" ")
        
        # record the time to our array
        Times[j-1] = float(Line[0])

        # record the RSL to our array
        RSL[j-1] = float(Line[1])

    
    plt.plot(Times[1:],RSL[1:],'k-',lw=1.5,label='Bideford RSL')
    plt.xlim(np.max(Times),np.min(Times))

    #Plotting several RSL curves
    File2 = "../driver_files/SY_RSL.data"
    Time,RSL = np.loadtxt(File2,unpack=True,skiprows=1,usecols=(0,1),delimiter=" ")
    plt.plot(Time,RSL,'r--',lw=1.5, label='Scalby RSL')

    File3 = "../driver_files/DB_RSL.data"
    Time3,RSL3 = np.loadtxt(File3,unpack=True,skiprows=1,usecols=(0,1),delimiter=" ")
    plt.plot(Time3,RSL3, 'b--',lw=1.5, label='Dunbar RSL')

    File4 = "../driver_files/GM_RSL.data"
    Time4,RSL4 = np.loadtxt(File4,unpack=True,skiprows=1,usecols=(0,1),delimiter=" ")
    plt.plot(Time4,RSL4, 'g--',lw=1.5, label='Glamorgan RSL')

    File5 = "../driver_files/SM_RSL.data"
    Time5,RSL5 = np.loadtxt(File5,unpack=True,skiprows=1,usecols=(0,1),delimiter=" ")
    plt.plot(Time5,RSL5, 'm--',lw=1.5, label='St. Margarets RSL')


    plt.ylabel("RSL (m)")
    plt.xlabel("Age (BP)")
    plt.grid()
    plt.ylim(-40,10)
     
    plt.legend(loc='lower right', numpoints=1)

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('plot_RSL_test_SY.png',dpi=300)

if __name__ == "__main__":
    FileName = "../../RPM_JRS/Test_4" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        
