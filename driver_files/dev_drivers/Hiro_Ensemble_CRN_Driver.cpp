/*------------------------------------------------------------------------

	Hiro_Ensemble_CRN_Driver.cpp
	
	Driver file for running the shore platform model of Matsumoto et al. (2016) with Cosmogenic Isotope accumulation (Hurst et al. 2017).
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model with Cosmogenic Isotope production.

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016a)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Matsumoto, H., Dickson, M.E., and Kench, P.S. (2016b)
	Modelling the Development of Varied Shore Profile Geometry on Rocky Coasts.
	Journal of Coastal Research http://dx.doi.org/10.2112/SI75-120.1

	Hurst, M.D., Rood, D.H., Ellis, M.A., Anderson, R.S., and Dornbusch, U. (2016)
	Recent acceleration in coastal cliff retreat rates on the south coast of Great Britain.
	Proceedings of the National Academy of Sciences, http://dx.doi.org/10.1073/PNAS.1613044113

	Hurst, M.D., Rood, D.H., and Ellis, M.A. (2017)
	Controls on the distribution of cosmogenic 10 Be across shore platforms
	Earth Surface Dynamics http://dx.doi.org/10.5194/esurf-5-67-2017

	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	
	March 2017
	
	Copyright (C) 2017, Martin Hurst
	
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

#include <fenv.h>
#include <signal.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <omp.h>
#include <unistd.h>
#include <Hiro.hpp>
#include <RockyCoastCRN.hpp>

using namespace std;

using namespace std;

template <typename T> string tostr(const T& t)
{ 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}

int main()
{
	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;
	double Gradient = 1.;

	//Time control parameters
	double EndTime = 10000.;
	double Time = 0.;
	double TimeInterval = 1;

	//Print Control
	double PrintInterval = 100.;
	double PrintTime = Time;
	string OutputMorphologyFileName = "ShoreProfile.xz";
	string OutputConcentrationFileName = "Concentrations.xn";
	
	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX, Gradient);
	
	//Initialise Tidal Ranges to check
	double TidalPeriod = 12.42;
	vector<double> TidalRanges;
	TidalRanges.push_back(1.);
	TidalRanges.push_back(4.);
	TidalRanges.push_back(8.);
	TidalRanges.push_back(1.);
	TidalRanges.push_back(4.);
	TidalRanges.push_back(8.);
	TidalRanges.push_back(1.);
	TidalRanges.push_back(4.);
	TidalRanges.push_back(8.);
	
	//Initialise WaveHeights to check
	vector<double> WaveHeights;
	WaveHeights.push_back(1.);
	WaveHeights.push_back(1.);
	WaveHeights.push_back(1.);
	WaveHeights.push_back(2.);
	WaveHeights.push_back(2.);
	WaveHeights.push_back(2.);
	WaveHeights.push_back(3.);
	WaveHeights.push_back(3.);
	WaveHeights.push_back(3.);
	
	//Initialise weathering efficacy to check
	vector<double> WeatheringRates;
	WeatheringRates.push_back(0.01);
	WeatheringRates.push_back(0.001);
	WeatheringRates.push_back(0.01);
	WeatheringRates.push_back(0.0001);
	WeatheringRates.push_back(0.0001);
	WeatheringRates.push_back(0.001);
	WeatheringRates.push_back(0.0001);
	WeatheringRates.push_back(0.0001);
	WeatheringRates.push_back(0.001);

	//Initialise resistances to check
	vector<double> Resistances;
	Resistances.push_back(0.1);
	Resistances.push_back(0.1);
	Resistances.push_back(0.1);
	Resistances.push_back(0.001);
	Resistances.push_back(0.001);
	Resistances.push_back(0.01);
	Resistances.push_back(0.01);
	Resistances.push_back(0.1);
	Resistances.push_back(0.001);
	
	//Breaking Wave Coefficients
	double BreakingCoefficient = 1.;
	double BrokenCoefficient = 1.;
	double StandingCoefficient = 0.01;
	
	//Geology
	double CliffHeight = 10.;
	double CliffFailureDepth = 1.;
	
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 1.;
	double WaveHeight_StD = 0;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	
	//Which Nuclides to track 10Be, 14C, 26Al, 36Cl?
	vector<int> Nuclides;
	Nuclides.push_back(10);
	
	//Set sea level rise rate
	double SLR = 0;
	
	//initialise RockyCoastCRN friend object
	RockyCoastCRN BlankPlatformCRN = RockyCoastCRN(PlatformModel, Nuclides);
	RockyCoastCRN PlatformCRN = BlankPlatformCRN;
	
	//Loop across parameter space
	for (int a=0; a<9; ++a)
	{
		//setup the output file
		OutputMorphologyFileName = "ShoreProfile_G"+tostr(Gradient)
												+"_T_"+tostr(TidalRanges[a])
												+"_H_"+tostr(WaveHeights[a])
												+"_W_"+tostr(WeatheringRates[a])
												+"_R_"+tostr(Resistances[a])
												+"_Br_"+tostr(BreakingCoefficient)
												+"_Bo_"+tostr(BrokenCoefficient)+".xz";
		
		OutputConcentrationFileName = "Concentrations_G"+tostr(Gradient)
												+"_T_"+tostr(TidalRanges[a])
												+"_H_"+tostr(WaveHeights[a])
												+"_W_"+tostr(WeatheringRates[a])
												+"_R_"+tostr(Resistances[a])
												+"_Br_"+tostr(BreakingCoefficient)
												+"_Bo_"+tostr(BrokenCoefficient)+".xn";
		
		//Reset the model
		PlatformModel.ResetModel();
		PlatformCRN = BlankPlatformCRN;
		
		//Initialise Tides
		PlatformModel.InitialiseTides(TidalRanges[a]);
		PlatformCRN.InitialiseTides(TidalRanges[a]/2.,TidalPeriod);
		PlatformModel.InitialiseSeaLevel(SLR);
	
		//Setup Morphology
		//Populate geometric metrics
		PlatformModel.UpdateMorphology();
	
		//Initialise Waves
		//Single Wave for now but could use the waveclimate object from COVE!?
		PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);
		
		//reset wave height
		WaveHeight_Mean = WaveHeights[a];
		PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

		//set wave coefficients
		PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient);
		
		//reset the geology
		PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistances[a], WeatheringRates[a]);
						
		//run the model!
		cout << "Running Model with..." << endl;
		cout << "\tInitial Gradient " << setprecision(1) << Gradient << endl;
		cout << "\tTidal Range " << TidalRanges[a] << " m" << endl;
		cout << "\tWave Height " << WaveHeight_Mean << " m" << endl;
		cout << "\tBreaking Coefficient " << setprecision(3) << BreakingCoefficient << endl;
		cout << "\tBroken Coefficient " << BrokenCoefficient << endl;
		cout << "\tMax Weathering Rate " << WeatheringRates[a] << " m/yr" << endl;
		cout << "\tRock Resistance " << Resistances[a] << endl;
		
		Time = 0;
		EndTime = 10000;
		PrintTime = 0;
		PrintInterval = 100;
		TimeInterval = 1;
		
		//Loop through time
		while (Time <= EndTime)
		{
			//Update Sea Level
			PlatformModel.UpdateSeaLevel();

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

			//Update the Morphology 
			PlatformModel.UpdateMorphology();	  
		
			//Check for Mass Failure
			PlatformModel.MassFailure();
		
			//Update the morphology inside RockyCoastCRN
			PlatformCRN.UpdateMorphology(PlatformModel);
		  
			//Update the CRN concentrations
			PlatformCRN.UpdateCRNs();

			//print?
			if (Time >= PrintTime)
			{
				//Create string stream for puting time into filename
	//			stringstream ss;
	//			ss << PrintTime;
	//			string PrintTimeString = ss.str();
	//			ss.clear();
	//			string ConcentrationsFileName ="Concentrations_"+PrintTimeString+".xzn";
				//PlatformCRN.WriteNuclideArray(ConcentrationsFileName, Time, Nuclides[0]);
				PlatformModel.WriteProfile(OutputMorphologyFileName, Time);
				PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
				PrintTime += PrintInterval;
			}

			//update time
			Time += TimeInterval;
		}
	
	}
	//a few blank lines to finish
	cout << endl << endl;
	
}


