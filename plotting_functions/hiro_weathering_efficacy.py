# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 01:25:08 2017

@author: martin
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#Make weathering shape function based on tidal duration
TideRange=1
dZ = 0.1
NTideValues = (int)(TideRange/dZ);
Elevations = np.arange(0,NTideValues)*dZ
WeatheringEfficacy = np.zeros(NTideValues)

print np.log(0.25)
Mean = 0.25
Std = 1
Mode = np.exp(Mean-Std)
print Mode

if (dZ):
    for i in range(1, NTideValues):
        WeatheringEfficacy[i] = np.exp(-(((np.log(Elevations[i])-0.25*TideRange)**2.)/TideRange))
        
plt.plot(WeatheringEfficacy, Elevations,'k-')
plt.xlabel("Weathering Efficacy")
plt.ylabel("Elevation (m)")
plt.show()