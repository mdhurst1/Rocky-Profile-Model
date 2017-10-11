/*------------------------------------------------------------------------

	Hiro_Driver.cpp
	
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

using namespace std;

int main()
{
	//initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;
	double Gradient = 1.;

	//Time control parameters
	double EndTime = 10000;
	double Time = 0.;
	double TimeInterval = 1;

	//Print Control
	double PrintInterval = 100;
	double PrintTime = Time;
	string OutputFileName = "ShoreProfile.xz";
	
	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX, Gradient);
	
	//Initialise Tides
	double TidalRange = 8.;
	PlatformModel.InitialiseTides(TidalRange);
		
	//Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 3.;
	double WaveHeight_StD = 0.;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

	//Sea level rise?
	double SLR = 0;
	PlatformModel.InitialiseSeaLevel(SLR);

	// Wave coefficient constant
	double StandingCoefficient = 0.01;
	double BreakingCoefficient = 1.;
	double BrokenCoefficient = 1.;
	PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient);

	//reset the geology
	double CliffHeight = 10.;
	double CliffFailureDepth = 1.;
	double Resistance = 0.01;
	double WeatheringRate = 0.001;
	PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistance, WeatheringRate);
				

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

		//Update the Morphology 
		PlatformModel.UpdateMorphology();	  
		
		//Implement Weathering
		PlatformModel.IntertidalWeathering();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();

		//Check for Mass Failure
		PlatformModel.MassFailure();
		
		//Update the Morphology 
		PlatformModel.UpdateMorphology();
				
		//print?
		if (Time >= PrintTime)
		{
			PlatformModel.WriteProfile(OutputFileName, Time);
			PrintTime += PrintInterval;
			//cout << endl;
		}
		
		//update time
		Time += TimeInterval;
		
	}
	
	//a few blank lines to finish
	cout << endl << endl;
	
}


