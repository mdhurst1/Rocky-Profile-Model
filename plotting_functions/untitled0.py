#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:50:43 2020

@author: mhurst
"""

import numpy as np
import matplotlib.pyplot as plt

WaveHeight = 3.
Distance = np.arange(0,200,0.1)
WaveHeights= np.zeros(len(Distance))

BreakingPointZInd = 1000
BreakingWaveDist = 28.370597979145998
AttenuationConst = 5.5

EndBreakers = Distance[BreakingPointZInd]+BreakingWaveDist

print(Distance[BreakingPointZInd]+BreakingWaveDist)

for i in range(0,len(Distance)):
     
    if Distance[i] <= Distance[BreakingPointZInd]:
         WaveHeights[i] = WaveHeight
         
    elif Distance[i] <= (Distance[BreakingPointZInd]+BreakingWaveDist):
        WaveHeights[i] = WaveHeight*np.exp(-AttenuationConst*(Distance[i]-Distance[BreakingPointZInd]))
        
    else:
        print(Distance[i])
        WaveHeights[i] = WaveHeight*np.exp(-AttenuationConst*BreakingWaveDist)*np.exp(-AttenuationConst*(Distance[i]-EndBreakers))


plt.plot(Distance,WaveHeights)
plt.plot(EndBreakers,0,'r.')
#Dist = 52.570597979146001
#Xz = 52.6


#print(-AttenuationConst*np.exp(Xz-Dist))