# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 04:40:15 2017

@author: martin
"""

import numpy as np
import math
import matplotlib.pyplot as plt

tide = 1.
Z = np.arange(0,tide,1)
Weathering = np.zeros(len(Z))

print math.ceil(tide/4.)

for i in range(2,math.ceil(tide/4)):
    Weathering[i] = np.exp(-((i-math.ceil(tide/4))**2)/(math.ceil(tide/2)))
for i in range(np.ceil(tide/4),tide):
    Weathering[i] = np.exp(-((i-math.ceil(tide/4))**2)/(tide*math.ceil(tide/10)))
    
plt.plot(Z,Weathering)
plt.show()