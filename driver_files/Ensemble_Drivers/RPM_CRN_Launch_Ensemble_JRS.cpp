/*------------------------------------------------------------------------
	RPM_CRN_Launch_Ensemble.cpp
	
	function to build and launch slurm scripts for running multiple simulations
    of the coupled RPM_CRN model with varying parameters in a one-at-a-time fashion
	
    Martin Hurst, February 2020
	Copyright (C) 2017, Hiro Matsumoto, Martin Hurst
	
	Developer can be contacted
	martin.hurst@glasgow.ac.uk
  
	Martin D. Hurst
	School of Geographical and Earth Sciences
	University of Glasgow
	Glasgow
	Scotland
	G12 8QQ
  
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "--------------------------------------------------------------" << endl;
	cout << "|  Rocky Profile Model (RPM)                                 |" << endl;
	cout << "|  This program launches an ensemble of model runs that      |" << endl;
	cout << "|  simulate the development of shore platforms following     |" << endl;
	cout << "|  model developed by Matsumoto et al. (2016)                |" << endl;
	cout << "|                                                            |" << endl;
	cout << "|  launches qsub jobs across a range of parameter values     |" << endl;
	cout << "|                                                            |" << endl;
	cout << "|  Implemented in C++ by Martin Hurst, University of Glasgow |" << endl;
	cout << "|  for coupling to RockyCoastCRN; model for predicting       |" << endl;
	cout << "|  cosmogenic radionuclide concentrations in shore platforms |" << endl;
	cout << "--------------------------------------------------------------" << endl;
	cout << endl;

	//Test for correct input arguments
	if (nNumberofArgs!=2)
	{
		cout << "Error: This program requires one input: " << endl;
		cout << " * A path to the folder where the model will be run" << endl;
		cout << "-----------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_Driver.out /ProjectFolder/" << endl;
		cout << "-----------------------------------------------------------" << endl;
		exit(EXIT_SUCCESS);
	}

	string Folder = argv[1];

        //declare constant parameters 
        int CRNFlag = 1;
        double Gradient = 1;
        double TidalRanges = 4;
        double SubtidalEfficacy = 0.005;
        double WaveHeight = 3;
	
        //Wave decay rate (Attenuation constant)
	vector<double> WaveAttenuationConst;
	WaveAttenuationConst.push_back(-1.921);
	WaveAttenuationConst.push_back(-1.783);
	WaveAttenuationConst.push_back(-1.587);

        //Initialise Resistance       (kg m^2 yr^-1)
	vector<double> Resistances;
	Resistances.push_back(1.277);
	Resistances.push_back(1.767);
	Resistances.push_back(2.364);  

	//Initialise WeatheringRate as a factor of resistance?
	vector<double> WeatheringRates;
	WeatheringRates.push_back(-4.208);
	WeatheringRates.push_back(-3.130); 
	WeatheringRates.push_back(-2.187);

	// set up a number to track runs
	int Run = 0;

   // for(int p=0, Np = WaveAttenuationConst.size(); p<Np; p+=2)
   // {
        // Track run number
      // ++Run;

        // setup the script
        //string write_sh;
        //char sh_name[128];
        //sprintf(sh_name, "/home/jrs17/Main_RPM/Rocky-Profile-Model/driver_files/Ensemble_Drivers/RPM_CRN_%d.sh", Run);
        
        //write_sh.open(sh_name);
        //write_sh << "#!/bin/bash" << endl;
        //write_sh << "#PBS -wd /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -M martin.hurst@glasgow.ac.uk" << endl;
        //write_sh << "#PBS -m abe" << endl;
        //write_sh << "#PBS -N Run" << Run << endl;
        //write_sh << "#PBS -l cput=12:00:00" << endl;
        //write_sh << "#PBS -l walltime=24:00:00" << endl;
        //write_sh << "#PBS -e /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -o /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        
        //write_sh << "" << endl;

        // set up command to launch the model
      //  string write_sh =  "./RPM_CRN_Ensemble.out " << Folder << " "  << "Ensemble "
      //              << CRNFlag << " "
      //              << Gradient << " "
      //              << TidalRanges << " "
      //              << SubtidalEfficacy << " "
      //              << WaveHeight << " "
      //              << WaveAttenuationConst[p] << " "
      //              << Resistances[1] << " "
      //              << WeatheringRates[1] << endl;
    //}

    


   // for(int n=0, Nn = Resistances.size(); n<Nn; n+=2)
    //{
        // Track run number
      //  ++Run;

        // setup the script
        //string write_sh;
        //char sh_name[128];
        //sprintf(sh_name, "/home/jrs17/Main_RPM/Rocky-Profile-Model/driver_files/Ensemble_Drivers/RPM_CRN_%d.sh", Run);

        //write_sh.open(sh_name);
        //write_sh << "#!/bin/bash" << endl;
        //write_sh << "#PBS -wd /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -M martin.hurst@glasgow.ac.uk" << endl;
        //write_sh << "#PBS -m abe" << endl;
        //write_sh << "#PBS -N Run" << Run << endl;
        //write_sh << "#PBS -l cput=12:00:00" << endl;
        //write_sh << "#PBS -l walltime=24:00:00" << endl;
        //write_sh << "#PBS -e /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -o /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        
        //write_sh << "" << endl;

        // set up command to launch the model
      //  string write_sh =  "./RPM_CRN_Ensemble.out " << Folder << " "  << "Ensemble "
      //              << CRNFlag << " "
       //             << Gradient << " "
      //              << TidalRanges << " "
     //               << SubtidalEfficacy << " "
     //               << WaveHeight << " "
     //               << WaveAttenuationConst[1] << " "
     //               << Resistances[n] << " "
     //               << WeatheringRates[1] << endl;
    //}



    for(int l=0, Nl = WeatheringRates.size(); l<Nl; l+=2)
    {
        // Track run number
        ++Run;

        // setup the script
        //string write_sh;
        //char sh_name[128];
        //sprintf(sh_name, "/home/jrs17/Main_RPM/Rocky-Profile-Model/driver_files/Ensemble_Drivers/RPM_CRN_%d.sh", Run);
        
        //write_sh.open(sh_name);
        //write_sh << "#!/bin/bash" << endl;
        //write_sh << "#PBS -wd /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -M martin.hurst@glasgow.ac.uk" << endl;
        //write_sh << "#PBS -m abe" << endl;
        //write_sh << "#PBS -N Run" << Run << endl;
        //write_sh << "#PBS -l cput=12:00:00" << endl;
        //write_sh << "#PBS -l walltime=24:00:00" << endl;
        //write_sh << "#PBS -e /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        //write_sh << "#PBS -o /home/jrs17/COAST/RPM_CRN_Ensembles/" << endl;
        
        //write_sh << "" << endl;

        // set up command to launch the model
        const char* write_sh = "./RPM_CRN_Ensemble.out " + Folder + " "  + "Ensemble "
                    + CRNFlag + " "
                    + Gradient + " "
                    + TidalRanges + " "
                    + SubtidalEfficacy + " "
                    + WaveHeight + " "
                    + WaveAttenuationConst[1] + " "
                    + Resistances[1] + " "
                    + WeatheringRates[l] + endl;
    }

    for (int job=1; job<=Run; ++job)
    {
        // loop through jobs and launch them all
        char LaunchCommand[256];
        //sprintf(launch,"qsub /home/jrs17/Main_RPM/Rocky-Profile-Model/driver_files/Ensemble_Drivers/RPM_CRN_%i.sh", job);
	sprintf(LaunchCommand,"write_sh"); //, job);
        system("write_sh");
    }

	//a few blank lines to finish
        cout << "write_sh" << endl << endl;
	cout << "All jobs launched" << endl << endl;
}
