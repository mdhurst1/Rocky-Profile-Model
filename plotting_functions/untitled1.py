#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:25:19 2020

@author: mhurst
"""

import numpy as np
import matplotlib.pyplot as plt

Exponents = np.array([-1e-06,-1e-05,-0.0001,-0.001,-0.01,-0.1,-1,-10,-100,-1000])
Slow = np.array([0.999999,0.99999,0.9999,0.999,0.99005,0.904837,0.367879,4.53999e-05,3.72008e-44,0])
Fast = np.array([0.998274,0.998266,0.998179,0.997315,0.988702,0.9057,0.368239,4.54502e-05,-4.31496e+33,1.22171e+28])

plt.figure()
plt.plot(-Exponents,Slow,'o-')
plt.plot(-Exponents,Fast,'o-')
plt.xscale('log')
plt.yscale('log')

