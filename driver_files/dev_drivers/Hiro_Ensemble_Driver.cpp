/*------------------------------------------------------------------------

	Hiro_Ensemble_Driver.cpp
	
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

	//initialise Hiro Model
	Hiro PlatformModel = Hiro(dZ, dX);
	
	//Initialise different starting gradients
	vector<double> Gradients;
	Gradients.push_back(0.);
	Gradients.push_back(1.);
	
	//Initialise Tidal Ranges to check
	vector<double> TidalRanges;
	TidalRanges.push_back(1.);
	TidalRanges.push_back(4.);
	TidalRanges.push_back(8.);
	
	//Initialise WaveHeights to check
	vector<double> WaveHeights;
	WaveHeights.push_back(1.);
	WaveHeights.push_back(2.);
	WaveHeights.push_back(3.);
	
	//Initialise weathering efficacy to check
	vector<double> WeatheringRates;
	WeatheringRates.push_back(0.0001);
	WeatheringRates.push_back(0.001);
	WeatheringRates.push_back(0.01);
	
	//Initialise resistances to check
	vector<double> Resistances;
	Resistances.push_back(0.001);
	Resistances.push_back(0.01);
	Resistances.push_back(0.1);
	
	//Breaking Wave Coefficients
	vector<double> BreakingCoefficients;
	BreakingCoefficients.push_back(10.);
	BreakingCoefficients.push_back(1.);
	
	//Broken Wave Coefficients
	vector<double> BrokenCoefficients;
	BrokenCoefficients.push_back(1.);
	BrokenCoefficients.push_back(0.1);
	
	//Standing wave coefficient constant
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
	
	
	//Loop across parameter space
	for (int a=0, A=Gradients.size(); a<A; ++a)
	{
		for (int b=0, B=TidalRanges.size(); b<B; ++b)
		{
			for (int c=0, C=WaveHeights.size(); c<C; ++c)
			{
				for (int d=0, D=WeatheringRates.size(); d<D; ++d)
				{
					for (int e=0, E=Resistances.size(); e<E; ++e)
					{
						for (int f=0, F=BreakingCoefficients.size(); f<F; ++f)
						{
							for (int g=0, G=BrokenCoefficients.size(); g<G; ++g)
							{
								//setup the output file
								PlatformModel.OutputFileName = "ShoreProfile_G"+tostr(Gradients[a])
																		+"_T_"+tostr(TidalRanges[b])
																		+"_H_"+tostr(WaveHeights[c])
																		+"_W_"+tostr(WeatheringRates[d])
																		+"_R_"+tostr(Resistances[e])
																		+"_Br_"+tostr(BreakingCoefficients[f])
																		+"_Bo_"+tostr(BrokenCoefficients[g])+".xz";
								
								//reset gradient
								PlatformModel.Set_InitialGradient(Gradients[a]);
								
								//reset the model
								PlatformModel.ResetModel();
		
								//Reset tidal range
								PlatformModel.InitialiseTides(TidalRanges[b]);
		
								//reset wave height
								WaveHeight_Mean = WaveHeights[c];
								PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);
		
								//set wave coefficients
								PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficients[f], BrokenCoefficients[g]);
								
								//reset the geology
								PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistances[e], WeatheringRates[d]);
												
								//run the model!
								cout << "Running Model with..." << endl;
								cout << "\tInitial Gradient " << setprecision(1) << Gradients[a] << endl;
								cout << "\tTidal Range " << TidalRanges[b] << " m" << endl;
								cout << "\tWave Height " << WaveHeight_Mean << " m" << endl;
								cout << "\tBreaking Coefficient " << setprecision(3) << BreakingCoefficients[f] << endl;
								cout << "\tBroken Coefficient " << BrokenCoefficients[g] << endl;
								cout << "\tMax Weathering Rate " << WeatheringRates[d] << " m/yr" << endl;
								cout << "\tRock Resistance " << Resistances[e] << endl;
								PlatformModel.EvolveCoast();
								cout << "\nDone\n\n";
							}
						}
					}
				}
			}
		}
	}	
	//a few blank lines to finish
	cout << endl << endl;
	
}


