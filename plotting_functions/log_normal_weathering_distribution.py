# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 02:53:09 2017

@author: martin
"""

import numpy as np
import matplotlib.pyplot as plt

TidalRange = 10.
NTidalValues = 11.
XMax = 10
#x = np.arange(0.001,XMax,0.001)
x = np.arange(0,NTidalValues)*10./NTidalValues

Theta = 0
m = 1.1665
#m = np.arange(0.01,1.,0.01)
sigma = 0.5
#m = sigma**2.
#Mode = np.exp(m-sigma**2.)
##np.log(Mode) = m-sigma**2.
##0 = m-sigma**2.
#m2 = sigma**2.
#print m2
#Mode = np.exp(m2-sigma**2.)
#print Mode
#plt.plot(m,Mode)
#plt.show()
P = np.zeros(len(x))
P[1:] = (np.exp(-((np.log(x[1:]-Theta)-m)**2.)/(2*sigma**2.))) / ((x[1:]-Theta)*sigma*np.sqrt(2*np.pi))
XPeak = 2.5
TheoreticalMax = np.exp(-((np.log(XPeak-Theta)-m)**2.)/(2*sigma**2.)) / ((XPeak-Theta)*sigma*np.sqrt(2*np.pi))
P = P/TheoreticalMax
print "theoretical max", TheoreticalMax
Elevation = np.linspace(TidalRange,0,len(x))
#plt.plot(x,P)
print len(P)
print len(Elevation)
plt.plot(P,Elevation)
plt.plot([0,1],[7.5,7.5],'k--')
plt.show()