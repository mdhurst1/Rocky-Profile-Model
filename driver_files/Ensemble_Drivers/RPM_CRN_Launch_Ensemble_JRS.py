# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:23:06 2019

@author: lguest
"""

# Import  modules
import sys, os
import subprocess
from subprocess import call
from subprocess import Popen
from yaml import safe_load

# Set physical constants as used by RPM - these are values read in command line
CRNFlag = 1
Gradient = 1
TidalRanges = 5.82
SubtidalEfficacy = 0.005
WaveHeight = 3

#Wave decay rate (Attenuation constant)
#WaveAttenuationConst = [-1.333]
WaveAttenuationConst = []
WaveAttenuationConst.append(-1.537)
WaveAttenuationConst.append(-1.333)
WaveAttenuationConst.append(-1.152)

#Initialise Resistance       (kg m^2 yr^-1)
Resistances = [1.985]
#Resistances = []
#Resistances.append(1.360)
#Resistances.append(1.985)
#Resistances.append(2.631)
  
#Initialise WeatheringRate as a factor of resistance?
WeatheringRates = [-2.644]
#WeatheringRates = []
#WeatheringRates.append(-3.078)
#WeatheringRates.append(-2.644)
#WeatheringRates.append(-2.226)

#set folder location 
Folder = "/home/jrs17/Main_RPM/Rocky-Profile-Model/driver_files/"

#set up a number to track runs 
Run = 0

#for(int p=0, Np = WaveAttenuationConst.size(); p<Np; p+=2)
for i in range(0, len(WaveAttenuationConst)):
	for j in range(0, len(Resistances)):
		for k in range(0, len(WeatheringRates)):
			Launchstr = Folder+"/Ensemble_Drivers/RPM_CRN_Ensemble.out "+ Folder+" "+"Ensemble "+ str(CRNFlag) +" "+ str(Gradient) +" "+ str(TidalRanges) +" "+ str(SubtidalEfficacy) +" "+ str(WaveHeight) +" "+ str(Resistances[j]) +" "+ str(WaveAttenuationConst[i]) +" "+ str(WeatheringRates[k])
			#Track Run
			Run +=1
			subprocess.call(Launchstr,shell=True)
			print(Launchstr)

