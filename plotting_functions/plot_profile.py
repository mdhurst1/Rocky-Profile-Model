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

    #load the extracted shore platform profile
    ExProfileName = "../driver_files/Data/CB_profile.txt"
    Xprof, Zprof = np.loadtxt(ExProfileName, unpack=True,skiprows=1,usecols=(0,1))
    XCliffPosition = Xprof[-1]
    Xprof -= XCliffPosition
   
    plt.ylabel("Z (m)")
    plt.xlabel("X (m)")
    plt.grid()
    plt.plot(Xprof[1:],Zprof[1:],'k-')
   
   
    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('plot_profile_CB.png',dpi=300)


    if __name__ == "__main__":
    #FileName = "../../RPM_JRS/Test_4" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        
