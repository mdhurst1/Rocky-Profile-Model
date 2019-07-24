/*------------------------------------------------------------------------

	RPM_CRN_Ensemble_Driver.cpp
	
	Driver file for running the shore platform model of Matsumoto et al. (2016)
	Updated following improvements by Matsumoto et al. (2018)

    This driver file takes model parameters as command line arguments to enable running
    ensembles on a HPC using TORQUE/PBS qsub command
	
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
	
	July 2019
	
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
#include "../RPM.hpp"
#include "../RoBoCoP_CRN/RockyCoastCRN.hpp"
#include "../SeaLevel.hpp"

using namespace std;

template <typename T> string tostr(const T& t)
{ 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}


int main(int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "|  Rocky Profile Model (RPM)                                                     |" << endl;
	cout << "|  This program models the development of shore platforms                        |" << endl;
	cout << "|  following model developed by Matsumoto et al. (2016)                          |" << endl;
	cout << "|                                                                                |" << endl;
	cout << "|  Implemented in C++ by Martin Hurst, University of Glasgow                     |" << endl;
	cout << "|  for coupling to RockyCoastCRN; model for predicting                           |" << endl;
	cout << "|  cosmogenic radionuclide concentrations in shore platforms                     |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << endl;

	//Test for correct input arguments
	if (nNumberofArgs!=11)
	{
		cout << "Error: This program requires 10 (YES TEN, one-zero) command line inputs: " << endl;
		cout << " * First a path to the folder where the model will be run" << endl;
		cout << " * The name of the project/model run" << endl;
		cout << " * A Flag to run with CRNs (1 = True)" << endl;
        cout << " * The initial topographic gradient" << endl;
        cout << " * A rate of sea level rise/fall (positive/negative respectively) (m/yr)" << endl;
        cout << " * The tidal range (m)" << endl;
        cout << " * The Maximum weathering rate (kg/m2/yr)" << endl;
        cout << " * The subtidal weathering efficacy (multiplier)" << endl;
        cout << " * The rock resistance (kg/m2)" << endl;
        cout << " * The wave attenuation constant" << endl;
        cout << " * " << endl;
		cout << "-----------------------------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_Driver.out /ProjectFolder/ Waipapa 1 1. 0.001 4. 0.005 0.01 1000. 10." << endl;
		cout << "-----------------------------------------------------------------------------" << endl;
        cout << endl;
		exit(EXIT_SUCCESS);
	}

    // read parameters from command line arguments
	string Folder = argv[1];
	string Project = argv[2];
	int CRNFlag = atoi(argv[3]);
	double Gradient = atof(argv[4]);
    double SLR = atof(argv[5]);
	double TidalRange = atof(argv[6]);
	double WeatheringRate = atof(argv[7]);
    double SubtidalEfficacy = atof(argv[8]);
    double Resistance = atof(argv[9]);
    double WaveAttenuationConst = atof(argv[10]);

	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;
	double CliffHeight = 15.;
	double MinElevation = -15.;

	//Time control parameters
	//Time runs in yrs bp
	double EndTime = 10000;
	double Time = 0.;
	double TimeInterval = 1;

	//Print Control
	double PrintInterval = 100;
	double PrintTime = Time-PrintInterval;
	
    //setup the output file                  
    string OutputMorphologyFileName = Folder+Project+"ShoreProfile_G"+tostr(Gradient)
                                                +"_S_"+tostr(SLR)
                                                +"_T_"+tostr(TidalRange)
                                                +"_W_"+tostr(WeatheringRate)
                                                +"_Ws_"+tostr(SubtidalEfficacy)
                                                +"_R_"+tostr(Resistance)
                                                +"_A_"+tostr(WaveAttenuationConst)+".xz";
                                        
    string OutputConcentrationFileName = Folder+Project+"Concentrations_G"+tostr(Gradient)
                                                +"_S_"+tostr(SLR)
                                                +"_T_"+tostr(TidalRange)
                                                +"_W_"+tostr(WeatheringRate)
                                                +"_Ws_"+tostr(SubtidalEfficacy)
                                                +"_R_"+tostr(Resistance)
                                                +"_A_"+tostr(WaveAttenuationConst)+".xn";
                                                
	//initialise RPM Model
	RPM PlatformModel = RPM(dZ, dX, Gradient, CliffHeight, MinElevation);
	
	//initialise RockyCoastCRN friend object
	RockyCoastCRN PlatformCRN = RockyCoastCRN();

	if (CRNFlag)
	{
		//Which Nuclides to track 10Be, 14C, 26Al, 36Cl?
		vector<int> Nuclides;
		Nuclides.push_back(10);
		
		//initialise RockyCoastCRN friend object
		PlatformCRN = RockyCoastCRN(PlatformModel, Nuclides);
	}
	
	// Initialise Sea level from datafile
	//string RelativeSeaLevelFile = Folder + Project + "_RSL.tz";
	//SeaLevel RelativeSeaLevel = SeaLevel(RelativeSeaLevelFile);
	
	// initialise sea level using rate of change
	SeaLevel RelativeSeaLevel = SeaLevel(SLR);
	
	// Get initial sea level
	double InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);
	PlatformModel.UpdateSeaLevel(InstantSeaLevel);

	//Initialise Tides
	double TidalPeriod = 12.42;
	PlatformModel.InitialiseTides(TidalRange);
    if (CRNFlag) PlatformCRN.InitialiseTides(TidalRange/2.,TidalPeriod);
		
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 3.;
	double WaveHeight_StD = 0.;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);
	
	//Tectonic Events
	//double UpliftFrequency = 1000.;
	//double UpliftTime = UpliftFrequency;
	//double UpliftMagnitude = 6.5;

	// Wave coefficient constant
	double StandingCoefficient = 0.1;
	double BreakingCoefficient = 10.;
	double BrokenCoefficient = 1.;
	PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient, WaveAttenuationConst);

	//reset the geology
	double CliffFailureDepth = 0.1;
	PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistance, WeatheringRate, SubtidalEfficacy);	

	// print initial condition to file
	double TempTime = -9999;
	PlatformModel.WriteProfile(OutputMorphologyFileName, TempTime);			
	if (CRNFlag) PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, TempTime);

	//Loop through time
	while (Time <= EndTime)
	{
		//Update Sea Level
		InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);
		PlatformModel.UpdateSeaLevel(InstantSeaLevel);

		//Get the wave conditions
		PlatformModel.GetWave();

		//Calculate forces acting on the platform
		PlatformModel.CalculateBackwearing();
		PlatformModel.CalculateDownwearing();

		//Do erosion
		PlatformModel.ErodeBackwearing();
		PlatformModel.ErodeDownwearing();

		//Update the Morphology 
		PlatformModel.UpdateMorphology();	
		
		//Implement Weathering
		PlatformModel.IntertidalWeathering();
		PlatformModel.SubtidalWeathering();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();

		//Check for Mass Failure
		PlatformModel.MassFailure();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();

        //Update the morphology inside RockyCoastCRN
		if (CRNFlag) PlatformCRN.UpdateMorphology(PlatformModel);

		//Update the CRN concentrations
		if (CRNFlag) PlatformCRN.UpdateCRNs();
        	
		//print?
		if (Time >= PrintTime)
		{
			PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
			if (CRNFlag) PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
			PrintTime += PrintInterval;
		}
		
		//update time
		Time += TimeInterval;
	}


}


