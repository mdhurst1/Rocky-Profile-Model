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
#include <map>
//#include <omp.h>
#include <unistd.h>
#include "../RPM.hpp"
#include "../Parameters.hpp"
#include "../RoBoCoP_CRN/RockyCoastCRN.hpp"
#include "../SeaLevel.hpp"

using namespace std;

int main(int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "--------------------------------------------------------------" << endl;
	cout << "|  Rocky Profile Model with Cosmogenic RadioNuclides (RPM_CRN) |" << endl;
	cout << "|  This program models the development of shore platforms      |" << endl;
	cout << "|  following model developed by Matsumoto et al. (2016), and   |" << endl;
	cout << "|  the accumulation of cosmogenic isotopes following Hurst     |" << endl;
	cout << "|  et al. (2016; 2017)                                         |" << endl;
	cout << "|                                                              |" << endl;
	cout << "|  Implemented in C++ by Martin Hurst, University of Glasgow   |" << endl;
	cout << "|  for coupling to RockyCoastCRN; model for predicting         |" << endl;
	cout << "|  cosmogenic radionuclide concentrations in shore platforms   |" << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << endl;

	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "Error: This program requires two inputs: " << endl;
		cout << " * First a path to the folder where the model will be run" << endl;
		cout << " * The name of the input parameter file (must be in the folder where the model will be run)" << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_Driver.out /ProjectFolder/ Waipapa.in" << endl;
		cout << "------------------------------------------------------" << endl;
		exit(EXIT_SUCCESS);
	}

	string Folder = argv[1];
	string TempParamFilename = argv[2];
	string InputParamFilename = Folder+TempParamFilename;
	
	// load parameter parser object
  	Parameters Params(Folder,InputParamFilename);

	//initialisation parameters, these are currently not in parameters object
	double dZ = 0.1;
	double dX = 0.1;
	
	//Time control parameters
	//Time runs in yrs bp
	int Time = Params.StartTime;
	int PrintTime = Time-Params.PrintInterval;
	
	//initialise RPM Model
	RPM PlatformModel = RPM(dZ, dX, Params.InitialGradient, Params.CliffHeight, Params.MinElevation);
	
	//initialise RockyCoastCRN friend object
	RockyCoastCRN PlatformCRN;

	// THIS SHOULD BE IN PARAMETER FILE
	if (Params.CRN_Predictions)
	{
		//Which Nuclides to track 10Be, 14C, 26Al, 36Cl?
		vector<int> Nuclides;
        if (Params.Berylium) Nuclides.push_back(10);
        if (Params.Carbon) Nuclides.push_back(14);
        if (Params.Aluminium) Nuclides.push_back(26);
		
		//initialise RockyCoastCRN friend object
		PlatformCRN = RockyCoastCRN(PlatformModel, Nuclides);
	}
	
	// Initialise Sea level from datafile
	SeaLevel RelativeSeaLevel;
	if (Params.ReadSeaLevelFromFile) RelativeSeaLevel = SeaLevel(Params.SeaLevelFilename);
	else RelativeSeaLevel = SeaLevel(Params.SeaLevelRise);
	
	// Get initial sea level
	float InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);
	PlatformModel.UpdateSeaLevel(InstantSeaLevel);

	//Initialise Tides
	PlatformModel.InitialiseTides(Params.TidalRange);
    if (Params.CRN_Predictions) PlatformCRN.InitialiseTides(Params.TidalRange/2.,Params.TidalPeriod);
		
	//Initialise Waves
	PlatformModel.InitialiseWaves(Params.WaveHeight_Mean, Params.WaveHeight_StD, Params.WavePeriod_Mean, Params.WavePeriod_StD);
	
	//Tectonic Events
	//double UpliftFrequency = 1000.;
	//double UpliftTime = UpliftFrequency;
	//double UpliftMagnitude = 6.5;

	// Wave coefficient constant
	PlatformModel.Set_WaveCoefficients(Params.StandingWaveCoef, Params.BreakingWaveCoef, 
										Params.BrokenWaveCoef, Params.WaveAttenuationConst);

	//reset the geology
	PlatformModel.InitialiseGeology(Params.CliffHeight, Params.CliffFailureDepth, Params.Resistance, 
									Params.WeatheringRate, Params.SubtidalEfficacy);	

	// print initial condition to file
	PlatformModel.WriteProfile(Params.ProfileOutFilename, Params.StartTime);			
	if (Params.CRN_Predictions) PlatformCRN.WriteCRNProfile(Params.ConcentrationsOutFilename, Params.StartTime);
	
	//Loop through time
	while (Time >= Params.EndTime)
	{
		//Do an earthquake?
		//if (Time < UpliftTime)
		//{
			//string ArrayFile1 = "MorphArray1.data";
			//string ArrayFile2 = "MorphArray2.data";

			//PlatformModel.WriteMorphologyArray(ArrayFile1, Time);
			//PlatformModel.TectonicUplift(UpliftMagnitude);
			//UpliftTime -= UpliftFrequency;
			//PlatformModel.WriteMorphologyArray(ArrayFile2, Time);

			//Update the Morphology 
			//PlatformModel.UpdateMorphology();
		//}		
		
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

		//Implement Weathering
		PlatformModel.IntertidalWeathering();
		PlatformModel.SubtidalWeathering();
		
		//Check for Mass Failure
		PlatformModel.MassFailure();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();

        //Update the morphology inside RockyCoastCRN
		if (Params.CRN_Predictions) 
		{
			PlatformCRN.UpdateMorphology(PlatformModel);
			PlatformCRN.UpdateCRNs();
		}
        	
		//print?
		if (Time <= PrintTime)
		{
			PlatformModel.WriteProfile(Params.ProfileOutFilename, Time);
			if (Params.CRN_Predictions) PlatformCRN.WriteCRNProfile(Params.ConcentrationsOutFilename, Time);
			PrintTime -= Params.PrintInterval;
		}
		
		//update time
		Time -= Params.TimeStep;
		
	}
	
	//a few blank lines to finish
	cout << endl << "Done" << endl << endl;

}
