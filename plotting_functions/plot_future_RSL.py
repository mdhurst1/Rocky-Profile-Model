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
rc('xtick.major',pad=4)
padding = 1

#create blank figure
plt.figure(1,figsize=(6.6,3.3))

File = "../driver_files/Data/CB_RSL_future.data"
Time,RSL = np.loadtxt(File,unpack=True,skiprows=1,usecols=(0,1),delimiter=" ")
plt.plot(Time,RSL,'k--',lw=1.5, label='Bideford RSL')

plt.ylabel("RSL (m)")
plt.xlabel("Age (BP)")
plt.grid()
plt.ylim(0,0.5)
plt.xlim(1900,2100)
     
plt.legend(loc='upper left', numpoints=1)

fig1 = plt.gcf()
plt.show()
plt.draw()
#fig1.savefig('plot_RSL_SY_CB_talk.png',dpi=300)