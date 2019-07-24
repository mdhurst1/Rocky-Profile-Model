/*------------------------------------------------------------------------

	RPM_CRN_Driver.cpp
	
	Driver file for running the shore platform model of Matsumoto et al. (2016)
	Updated following improvements by Matsumoto et al. (2018)
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model with coupling to Cosmogenic Isotope production by RockCoastCRN/RoBoCoP.

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016a)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Matsumoto, H., Dickson, M.E., and Kench, P.S. (2016b)
	Modelling the Development of Varied Shore Profile Geometry on Rocky Coasts.
	Journal of Coastal Research http://dx.doi.org/10.2112/SI75-120.1

	Matsumoto, H., Dickson, M. E., Kench, P. S., (2018)
	Modelling the relative dominance of wave erosion and weathering processes in shore platform development in micro- to mega-tidal settings
	Earth Surface Processes and Landforms  http://dx.doi.org/10.1002/esp.4422
	
    Hurst, M.D., Rood, D.H., Ellis, M.A., Anderson, R.S., and Dornbusch, U. (2016) Recent acceleration in coastal cliff retreat rates on the south coast of Great Britain. Proceedings of the National Academy of Sciences, http://dx.doi.org/10.1073/PNAS.1613044113

	Hurst, M.D., Rood, D.H., and Ellis, M.A. (2017)
	Controls on the distribution of cosmogenic 10 Be across shore platforms
	Earth Surface Dynamics http://dx.doi.org/10.5194/esurf-5-67-2017

	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	Mark Dickson, University of Auckland
	
	March 2017
	
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
	if (nNumberofArgs!=4)
	{
		cout << "Error: This program requires two inputs: " << endl;
		cout << " * First a path to the folder where the model will be run" << endl;
		cout << " * The name of the project/model run" << endl;
		cout << " * A Flag to run with CRNs (1 = True)" << endl;
		cout << "-----------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_Driver.out /ProjectFolder/ Waipapa" << endl;
		cout << "-----------------------------------------------------------" << endl;
		exit(EXIT_SUCCESS);
	}

	string Folder = argv[1];
	string Project = argv[2];
	int CRNFlag = atoi(argv[3]);
	
	// Changing parameters for exploratory iterarions
    //Initial Gradient
	vector<double> Gradient;
	Gradient.push_back(1.);
	Gradient.push_back(0.1);
	Gradient.push_back(0);   
	// 0 Gradient sets horizontal slope of cliff top: inital array of 1's(rock) vertical cliff

	// Initialise sea level using rate of change
	vector<double> SLR;
	SLR.push_back(0.001);
	SLR.push_back(0.);
	SLR.push_back(-0.001);

	//Initialise Tidal range
	vector<double> TidalRanges;
	TidalRanges.push_back(1.);
	TidalRanges.push_back(4.);
	TidalRanges.push_back(8.);
	
	//Initialise WeatheringRate   (kg m^2 yr^-1)
	vector<double> WeatheringRates;
	WeatheringRates.push_back(0.005);
	WeatheringRates.push_back(0.05);
	WeatheringRates.push_back(0.5);

	// sets relative efficacy of subtidal weathering	
	vector<double> SubtidalEfficacy;
	SubtidalEfficacy.push_back(0.001);
	SubtidalEfficacy.push_back(0.01);
	SubtidalEfficacy.push_back(0.1);

	//Initialise Resistance       (kg m^2 yr^-1)
	vector<double> Resistances;
	Resistances.push_back(10.);
	Resistances.push_back(100.);
	Resistances.push_back(1000.);

	//Wave decay rate (Attenuation constant)
	vector<double> WaveAttenuationConst;
	WaveAttenuationConst.push_back(0.01);
	WaveAttenuationConst.push_back(0.1);
	WaveAttenuationConst.push_back(1.);

	// set up a number to track runs
	int Run = 0;

    cout << "I am here" << endl;

	//loop across parameter space varying one-at-a-time
	for (int i=0, Ni = Gradient.size(); i<Ni; ++i)
	{
	    cout << "Now I am here" << endl;

        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[i] << " "
                    << SLR[1] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[1] << endl;
    
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for (int j=0, Nj = SLR.size(); j<Nj; ++j)
	{
        cout << "and now I am here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[j] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[1] << endl;
    
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for(int k=0, Nk = TidalRanges.size(); k<Nk; ++k)
    {
        cout << "and now... I am here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[1] << " "
                    << TidalRanges[k] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[1] << endl;
    
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for(int l=0, Nl = WeatheringRates.size(); l<Nl; ++l)
    {
        cout << "I have made it to here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[1] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[l] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[1] << endl;
    
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for(int m=0, Nm = SubtidalEfficacy.size(); m<Nm; ++m)
    {
        cout << "and I have made it to here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[1] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[m] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[1] << endl;
    
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for(int n=0, Nn = Resistances.size(); n<Nn; ++n)
    {
        cout << "and now I have made it to here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[1] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[n] << " "
                    << WaveAttenuationConst[1] << endl;
        
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }

    for(int o=0, No = WaveAttenuationConst.size(); o<No; ++o)
    {
        cout << "and now... I have made it to here" << endl;
        
        // Track run number
        ++Run;

        // setup the script
        stringstream SS;
        
        // set up command to launch the model
        SS << "./RPM_CRN_Ensemble.out /export/home/mh322u/RPM_CRN_Ensembles/ Ensemble "
                    << CRNFlag << " "
                    << Gradient[1] << " "
                    << SLR[1] << " "
                    << TidalRanges[1] << " "
                    << WeatheringRates[1] << " "
                    << SubtidalEfficacy[1] << " "
                    << Resistances[1] << " "
                    << WaveAttenuationConst[o] << endl;
        
        string LaunchString = SS.str();
        cout << "Launching command: " << endl;
        cout << "\t" << LaunchString << endl;
        system(LaunchString.c_str());
    }
        
    //a few blank lines to finish
	cout << "All jobs launched" << endl << endl;
}
