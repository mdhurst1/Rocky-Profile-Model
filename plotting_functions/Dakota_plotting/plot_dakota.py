# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=9)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

Folder = "/home/jrs17/Dakota_Results/"
Files = [Folder+'Run_1'] #,Folder+'Run_2',Folder+'Run_3',Folder+'Run_4',Folder+'Run_5',Folder+'Run_6',Folder+'Run_7',Folder+'Run_8',Folder+'Run_9']
NumFiles = len(Files)

ColourMap = cm.viridis

for i, FileName in enumerate(Files):


    #create blank figure
    plt.figure(1,figsize=(8,4))

    #First load the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"_ShoreProfile.xz"
    f = open(ProfileName,'r')
    MorphLines = f.readlines()
    NoLines = len(MorphLines)
    StartTime = float(MorphLines[2].strip().split(" ")[0])
    EndTime = float(MorphLines[-1].strip().split(" ")[0])
    f.close()

    # Get last cliff position
    LastOutput = MorphLines[-1].strip().split(" ") 
    LastProfile = np.array(LastOutput[2:],dtype="float64")
    LCPosition = float(LastOutput[3])   #cliff position 2 for SY, 3 for CB - check this   
    LastProfile = (LCPosition-LastProfile)*-1

    #Need to normalise modelled topography to furthest cliff position?

    #Read morphology
 
    # Get Z Values
    HeaderLine = MorphLines[0].strip().split(" ")
    CliffHeight = float(HeaderLine[0])
    MinElev = float(HeaderLine[1])
    dZ = float(HeaderLine[2])
    NValues = (int)((CliffHeight-MinElev)/dZ+1)
    Z = np.linspace(CliffHeight, MinElev, NValues)

    Z_LastProfile = np.linspace(CliffHeight, MinElev, len(LastProfile))
    
    #Second load CRN concentrations through time
    # declare the file and the axis
    ProfileName = FileName+"Concentrations.xn"
    f = open(ProfileName,'r')
    NLines = f.readlines()
    #NNoLines = len(NLines)
    #NEndTime = float(NLines[-1].strip().split(" ")[0])
    f.close()
    
    # Only plot every 1 000 years
    PlotTime = 8000
    PlotInterval = 100

    
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    ax1.plot(LastProfile,Z_LastProfile,'-',color=ColourMap(float(i)/float(NumFiles)), lw=1.5) 
    
    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
    #    MorphLine = (MorphLines[j].strip().split(" "))
        
    #    Time = float(MorphLine[0]) 
    #    RSL = float(MorphLine[1])
        
        #Read morphology
    #    X = np.array(MorphLine[2:],dtype="float64")
    #    Z = np.linspace(CliffHeight, MinElev, len(X))
                
    #    if (j == 1):
    #        print(Time)
    #        print (np.shape(X))
    #        print (np.shape(Z))
    #        ax1.plot(X,Z,'--',lw=1.5,color=ColourMap((Time)/(StartTime)))   #,label='Modelled Morphology')
    #    if (Time == PlotTime):
    #        #print (Time,StartTime)
    #        ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time)/(StartTime)))             
    #        PlotTime -= PlotInterval

    
        N10Line = (NLines[-1].strip().split(" "))   #reading line -2?   
        N10 = np.array(N10Line[1:],dtype="float64")

        X2 = np.arange(0,len(N10))*0.1

        #normalise so cliff position = 0 
        X2 = (LCPosition-X2)*-1

        mask = [N10!=N10[-1]]
        N10 = N10[mask]
        X2 = X2[mask]


        ax2.plot(X2,N10,'k-',lw=1.5) #, label='Modelled CRN Concentrations')
  
    
#Normalising CRN x positions to final cliff position
File2 = "../driver_files/Data/CB_CRN.data"
X,CRN,Error=np.loadtxt(File2,unpack=True,skiprows=1,usecols=(1,2,3),delimiter=" ")

#NormalisedX = LCPosition - X
NormalisedX = 0 - X

ax2.errorbar(NormalisedX,CRN,fmt='o', yerr=Error,c='grey', label='Measured CRN Concentrations')    
ax2.scatter(NormalisedX,CRN, s=0.1)

#load the extracted shore platform profile (SY - profile read from highest to lowest Z)
ExProfileName = "../driver_files/Data/CB_profile.txt"
Xprof, Zprof = np.loadtxt(ExProfileName, unpack=True,skiprows=1,usecols=(0,1))

#XCliffPosition = Xprof[0]
#Xprof -= XCliffPosition
#Xprof = LCPosition - Xprof
Xprof = 0 - Xprof

    
#axis labels
ax1.set_ylabel("Elevation (m)")
ax2.set_ylabel("Concentration (a g$^-1$)")  #x 10$^3$ 
ax2.set_xlabel("Distance (m)")
#xmin, xmax = ax1.get_xlim()
#ax1.set_xlim(-300,0)
#ax2.set_xlim(-300,0) 
#ax1.set_ylim(-10,5)
#ax2.set_ylim(0,9000)
ax1.plot(Xprof,Zprof,'r-',lw=1.5, label='Extracted Morphology')  #NXprof for scalby, Xprof for CB

ax1.legend(loc='upper left', numpoints=1)
ax2.legend(loc='lower left', numpoints=1)
    
fig1 = plt.gcf()
plt.show()
plt.draw()
#fig1.savefig('Dakota_result_test.png',dpi=300)

        