# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 01:04:47 2017

Code developed by Hironori Matsumoto
 Last update : 28 February 2017

 CREATE TIDAL DURATION FUNCTION
### INPUT : num_tidal_range ... tidal range / gridsize
### INPUT : start_tide ... starting index of tidal range loop
### INPUT : end_tide ... ending index of tidal range loop  
### INPUT : gridsize ... gridsize
### OUTPUT : tr_esf ... tidal duration function

Pythonised by Martin

@author: martin
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


def make_tidal_range(num_tidal_range, start_tide, end_tide, gridsize):

    # INITIAL SETTING
    total = 0;
    num =  len(num_tidal_range)
    tr_max  = np.max(num_tidal_range)+1
    tr_esf  = np.zeros([tr_max,num])
    tr_int = num_tidal_range.astype(int)
    
    ### CHANGE BOTTOM DEPENDING ON GRIDSIZE
    if ( gridsize == 0.1 ): bottom = 1;
    else: bottom  = gridsize / 0.1;
    
    ### ERROR DISPLAY
    for i in range(0,num):
        if (np.ceil(num_tidal_range[i]) != np.floor(num_tidal_range[i])):
            sys.exit('tide/gridsize range must be integer number')
                
    ### DEVELOP
    for j in range(0,num):
        total=0
        tr = np.int(num_tidal_range[j])
        print "tr is ", tr
        
        if (tr < 20):
            if (np.ceil(0.5*tr)!=np.floor(0.5*tr)):
                print "am I here?", tr
                for i in range(0,np.ceil(0.5*tr)+1):
                    tr_esf[i,j] = np.sin((i)*np.pi/np.ceil(tr*0.5))
                    total = total + tr_esf[i,j]*bottom
                
                for i in range(np.floor(0.5*tr),tr):     
                    tr_esf[i,j] = tr_esf[i,j] + np.sin((i-np.floor(tr*0.5))*np.pi/np.ceil(tr*0.5))
                    total = total + tr_esf[i,j]*bottom
                
            else:
                print "am I here?", tr
                for i in range(0,np.int(0.5*tr)):
                    tr_esf[i,j] = np.sin((i)*np.pi/(0.5*tr))
                    total += tr_esf[i,j]*bottom
                    print i, tr_esf[i,j]
                
                for i in range(np.int(0.5*tr),tr_max):
                    tr_esf[i,j] = -np.sin((i)*np.pi/(0.5*tr))
                    total += tr_esf[i,j]*bottom
                    print i, tr_esf[i,j]
            
            tr_esf[:,j] = tr_esf[:,j]/total
            print tr_esf[:,j]
            
        elif ( tr >= 20):
            print "or am I here?", tr
            for i in range(0,np.int(np.ceil(0.55*tr))+1):
                tr_esf[i,j] = np.sin((i)*np.pi/np.ceil(tr*0.55))
                total = total + tr_esf[i,j]*bottom
            
            for i in range(np.int(np.floor(0.45*tr)),tr_max):
                tr_esf[i,j] = tr_esf[i,j] + np.sin((i-(tr-np.floor(tr*0.55)))*np.pi/np.ceil(tr*0.55))
                total = total + tr_esf[i,j]*bottom
            
            tr_esf[:,j] = tr_esf[:,j]/total
    
        print tr_esf[:,0]
    
    #plot results
    for k in range(0,len(tr_int)):
        Weights = tr_esf[:,k]
        Weights = Weights[0:tr_int[k]+1]
        Elevations = np.arange(0.5*tr_int[k],-0.5*tr_int[k]-1,-1)*gridsize
        plt.plot(Weights,Elevations)

    plt.xlabel("Tidal Distribution")
    plt.ylabel("Elevation (m)")
    plt.show()

if __name__ == "__main__":
    tidal_range = np.array([1, 2, 4, 8],dtype=np.int32)
    start_tide = 1
    end_tide = 1
    gridsize = 0.1
    num_tidal_range = (tidal_range/gridsize).astype(int)
    make_tidal_range(num_tidal_range, start_tide, end_tide, gridsize)