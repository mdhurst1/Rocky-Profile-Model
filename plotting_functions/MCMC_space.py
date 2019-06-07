# -*- coding: utf-8 -*-

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
import scipy.stats as st
from scipy.stats import kde

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 1

def make_plot(FileName,ColourMap):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,4))
    area = np.pi*3

    #load the MCMC data and extract x and y
    MCMCFileName = "../driver_files/test4.out"
    Resistance_New, WeatheingRate_New = np.loadtxt(MCMCFileName, unpack=True,skiprows=1,usecols=(1,2), delimiter=" ")
    x = Resistance_New
    y = WeatheingRate_New

    #plot scatter
    #plt.scatter(x, y, s=area, alpha=0.5)

    #plot 2D density plot
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    nbins=300
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
    #plt.show()
 
    # Change colour palette
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.afmhot) #_r to reverse colour gradient
    #plt.colorbar()

    #plotting functions for scatter
    plt.ylabel(" Weathering Rate (K) ")
    plt.xlabel(" Resistance (FR) ")
    plt.xlim(0.025,0.05)
    plt.ylim(0.00046,0.00051)


   
    fig1 = plt.gcf()
    plt.show()
    #plt.draw()
    fig1.savefig('MCMC_SY_test1.png',dpi=300)



if __name__ == "__main__":
     FileName = "../driver_files/test4.out" # /Users/jennyshadrick/RPM_JRS
     ColourMap = cm.gray
     make_plot(FileName,ColourMap)
        