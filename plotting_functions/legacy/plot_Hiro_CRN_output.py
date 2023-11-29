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
rc('font',size=10)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

def make_plot(FileName,ColourMap):
    
    """ some comments about what the function does """
    
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

    LastProfile = LCPosition-LastProfile
    
    # Get Z Values
    HeaderLine = MorphLines[0].strip().split(" ")
    CliffHeight = float(HeaderLine[0])
    MinElev = float(HeaderLine[1])
    dZ = float(HeaderLine[2])
    NValues = (int)((CliffHeight-MinElev)/dZ+1)
    Z = np.linspace(CliffHeight, MinElev, NValues)
    
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
    
    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
        MorphLine = (MorphLines[j].strip().split(" "))
        
        Time = float(MorphLine[0]) 
        RSL = float(MorphLine[1])
        
        #Read morphology
        X = np.array(MorphLine[2:],dtype="float64")
        Z = np.linspace(CliffHeight, MinElev, len(X))
                
        if (j == 1):
            print(Time)
            print (np.shape(X))
            print (np.shape(Z))
            ax1.plot(X,Z,'--',lw=1.5,color=ColourMap((Time)/(StartTime)), label='Modelled Morphology')
        if (Time == PlotTime):
            #print (Time,StartTime)
            ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time)/(StartTime)))             
            PlotTime -= PlotInterval

    
    
    N10Line = (NLines[-1].strip().split(" "))   #reading line -2?   
    N10 = np.array(N10Line[1:],dtype="float64")
    X2 = np.arange(0,len(N10))*0.1
    mask = [N10!=N10[-1]]
    N10 = N10[mask]
    X2 = X2[mask]

           
    ax2.plot(X2,N10,'k-',lw=1.5, label='Modelled CRN')
  
        
    #tweak the plot
    #ax1.set_xticklabels([])

    #Normalising CRN x positions to final cliff position
    File2 = "../driver_files/Data/CB_CRN.data"
    X,CRN,Error=np.loadtxt(File2,unpack=True,skiprows=1,usecols=(1,2,3),delimiter=" ")
    NormalisedX = LCPosition - X
    ax2.errorbar(NormalisedX,CRN,fmt='o', yerr=Error,c='grey', label="Measured CRN")    
    ax2.scatter(NormalisedX,CRN, s=0.1)

    #File 3 uncorrected CRN data for Dunbar?
    #File3 = "../driver_files/DB_CRN_uncorrected.data"
    #X2,CRN2,Error2=np.loadtxt(File3,unpack=True,skiprows=1,usecols=(1,2,3),delimiter=" ")
    #NormalisedX2 = LCPosition - X2
    #ax2.errorbar(NormalisedX2,CRN2, fmt='o', yerr=Error2, c='k', label='Measured CRN Concentrations (uncorrected)')

    #load the extracted shore platform profile (CB - profile read from lowest to highest Z)
    #ExProfileName = "../driver_files/Data/Swath_Profile_CB.txt"
    #Xprof, Zprof = np.loadtxt(ExProfileName, unpack=True,skiprows=1,usecols=(0,1))
    #XCliffPosition = Xprof[-1]
    #Xprof -= XCliffPosition
    #Xprof += LCPosition

    #load the extracted shore platform profile (SY - profile read from highest to lowest Z)
    ExProfileName = "../driver_files/Data/CB_profile.txt"
    Xprof, Zprof = np.loadtxt(ExProfileName, unpack=True,skiprows=1,usecols=(0,1))
    #XCliffPosition = Xprof[0]
    #Xprof -= XCliffPosition
    Xprof = LCPosition - Xprof
    
    #axis labels
    ax1.set_ylabel("Elevation (m)")
    ax2.set_ylabel("Concentration (a g$^-1$)")  #x 10$^3$ 
    ax2.set_xlabel("Distance (m)")
    #xmin, xmax = ax1.get_xlim()
    ax1.set_xlim(1200,1500)
    ax2.set_xlim(1200,1500) 
    ax1.set_ylim(-10,10)
    ax2.set_ylim(0,20000)

    ax1.plot(Xprof,Zprof,'g-',lw=1.5, label='Extracted Morphology')  #NXprof for scalby, Xprof for CB

   

    ax1.legend(loc='upper left', numpoints=1, fontsize=10)
    ax2.legend(loc='upper left', numpoints=1, fontsize=10)
    

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('best_fit.png',dpi=300)

if __name__ == "__main__":
    FileName = "/Users/jennyshadrick/Dakota_Runs/Test_4" # /Users/jennyshadrick/RPM_JRS
    ColourMap = cm.RdBu
    make_plot(FileName,ColourMap)
        