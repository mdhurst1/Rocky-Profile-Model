# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:23:06 2019

@author: lguest
"""

# Import  modules
import sys
import subprocess
from subprocess import call
from yaml import safe_load

# Import RPM model 
#subprocess.Popen(["./RPM_dakota.out"])

# Set physical constants as used by RPM - these are values read in command line

Gradient = 1
TidalRange = 8
SubtidalEfficacy = 0.005
WaveAttenuationConst = 0.01

#########################################
#                                       #
#    Step 1: Use Dakota created         #
#    input files to prepare for         #
#    model run.                         #
#                                       #
#########################################

input_template = "input_template.yml"
inputs = "inputs.yml"
call(["dprepro", sys.argv[1], input_template, inputs])
#call(['rm', input_template])

#########################################
#                                       #
#    Step 2: Run Model                  #
#                                       #
#########################################

# Load parameters from the yaml formatted input.
with open(inputs, "r") as f:
    params = safe_load(f)
    Resistance = params["Resistance"]
    WeatheringRate = params["WeatheringRate"]
  
# set up command to run model 
  
Launchstr = "../../RPM_dakota.out /dakotaQuesoWorking/Rocky-Profile-Model/driver_files/Dakota_Drivers/  "+ sys.argv[2] +" ../driver_files/Data/CB_profile.txt 0  "+ str(Gradient) +" "+ str(TidalRange) +" "+ str(SubtidalEfficacy) +" "+ str(WaveAttenuationConst) +" "+ str(Resistance) +" "+ str(WeatheringRate)
subprocess.Popen(Launchstr,shell=True)
print(Launchstr) 

#########################################
#                                       #
#    Step 3: Write Output in format     #
#    Dakota expects                     #
#                                       #
#########################################

# Each of the metrics listed in the Dakota .in file needs to be written to
# the specified output file given by sys.argv[2]. This is how information is
# sent back to Dakota.

# RPM already produces outputs expected by dakota (argv[2])
# Don't need to write it to expected file 

# Get RMSE/ Likelihood output file from RPM_dakota 
# RMSE/ Likelihood (objective function) currently calculated within RPM model 


  
