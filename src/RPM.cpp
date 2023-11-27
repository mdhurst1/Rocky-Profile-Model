/*------------------------------------------------------------------------

	RPM.cpp

	C++ implementation of Hiro Matsumoto's Rocky Platform Model
	with option to couple Cosmogenic Isotope production using RoBoCoP_CRN.

	Last modified 06-10-2017 by Hiro Matsumoto

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016a)
	An exploratory numerical model of rocky shore profile evolution.
	Geomorphology http://doi.org/10.1016/j.geomorph.2016.05.017

	Matsumoto, H., Dickson, M.E., and Kench, P.S. (2016b)
	Modelling the Development of Varied Shore Profile Geometry on Rocky Coasts.
	Journal of Coastal Research http://dx.doi.org/10.2112/SI75-120.1

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

#ifndef RPM_CPP
#define RPM_CPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
//#include <omp.h>
#include "RPM.hpp"
#include "SeaLevel.hpp"
#include "FastExp.hpp"
#include "RockyCoastCRN.hpp"

using namespace std;

void RPM::Initialise()
{
  /* initialise an empty RPM object as default */
  printf("\nRPM.Initialise: Initialised an empty RPM object\n");
}

void RPM::Initialise(double dZ_in, double dX_in)
{
	/* initialise a verictal cliff RPM object */
	printf("\nRPM.Initialise: Initialised a RPM as a vertical cliff\n");

	//PHYSICAL CONSTANTS
	rho_w = 1025.;
	g = 9.81;

	//Define these in an input parameter file?
	SubmarineDecayConst = 0.1;
	StandingWaveConst = 0.01;
	BreakingWaveConst = 10.;
	BrokenWaveConst = 1.;
	//BrokenWaveConst = 0.1;
	BreakingWaveDecay = 0.2;
	BrokenWaveDecay = 0.01;
	WeatheringConst = 0.005;
	RockResistance = 0.5;

	//Cliff control params
	CliffElevation = 10.;
	MaxElevation = CliffElevation;
	MinElevation = -CliffElevation;
	CliffFailureDepth = 1.;

	//Declare spatial stuff
	dZ = dZ_in;
	dX = dX_in;
	NXNodes = 1000;
	NZNodes = (int)(round(2.*CliffElevation/dZ)+1);

	//declare an array of zeros for assignment of Z vector
	Z = vector<double>(NZNodes,0.0);
	for (int i=0; i<NZNodes; ++i) Z[i] = dZ*(0.5*(NZNodes-1)-i);

	//declare arrays of zeros to initalise various other vectors
	vector<double> ZZeros(NZNodes,0);
	vector<double> XZeros(NXNodes,0);
	X = XZeros;
	for (int j=0; j<NXNodes; ++j) X[j] = dX*j;
	Zx = XZeros;
	Xz = ZZeros;

	Bw_Erosion = ZZeros;
	Dw_Erosion = ZZeros;
	Weathering = ZZeros;

	//Indices trackers
	XInd = vector<int>(NZNodes,0);
	ZInd = vector<int>(NXNodes,0);

	//declare an array of ones for assignment of the Morphology Array
	MorphologyArray = vector< vector<int> >(NZNodes,vector<int>(NXNodes,1));
	ResistanceArray = vector< vector<double> >(NZNodes,vector<double>(NXNodes,RockResistance));

	//default time interval
	Time = 0;
	TimeInterval = 1;
	EndTime = 10000;

	//default print conditions to print every timestep
	PrintTime = 0;
	PrintInterval = TimeInterval;

	//Set sea level to zero to begin with, and the ind, this will get updated later
	SeaLevelRise = 0;
	SeaLevel = 0;
	SeaLevelInd = 0;
}


void RPM::Initialise(double dZ_in, double dX_in, double Gradient, double CliffsElevation, double MaximumElevation, double MinimumElevation)
{
	/* initialise a sloping cliff RPM object */
	printf("\nRPM.Initialise: Initialised a RPM as a slope\n");

	//PHYSICAL CONSTANTS
	rho_w = 1025.;
	g = 9.81;

	//Define these in an input parameter file?
	SubmarineDecayConst = 0.1;
	StandingWaveConst = 0.01;
	BreakingWaveConst = 10.;
	BrokenWaveConst = 1.;
	//BrokenWaveConst = 0.1;
	BreakingWaveDecay = 0.1;
	BrokenWaveDecay = 0.01;
	WeatheringConst = 0.005;
	RockResistance = 0.5;

	//Cliff control params
	CliffFailureDepth = 1.;
	CliffElevation = CliffsElevation;
	MaxElevation = MaximumElevation;
	MinElevation = MinimumElevation;

	if (MaxElevation < CliffElevation)
	{
		 printf("\tWarning: MaxElevation is less than the CliffElevation\n");
	}

	//Declare spatial stuff
	dZ = dZ_in;
	dX = dX_in;
	InitialGradient = Gradient;
	
	// figure out initial domain size
	NZNodes = (int)(round((MaxElevation-MinElevation)/dZ)+1);
	if (InitialGradient > 0) NXNodes = NZNodes/InitialGradient;
	else NXNodes = 1000.;
	
	//declare an array of zeros for assignment of Z vector
	vector<double> TempZ(NZNodes,0);
	Z = TempZ;
	for (int i=0; i<NZNodes; ++i)   Z[i] = MaxElevation-i*dZ;

	//declare arrays of zeros to initalise various other vectors
	vector<double> ZZeros(NZNodes,0);
	vector<double> XZeros(NXNodes,0);
	X = XZeros;
	for (int j=0; j<NXNodes; ++j)   X[j] = dX*j;

    Zx = XZeros;
	Xz = ZZeros;

	Bw_Erosion = ZZeros;
	Dw_Erosion = ZZeros;
	Weathering = ZZeros;

	//Indices trackers
	XInd = vector<int>(NZNodes,0);
	ZInd = vector<int>(NXNodes,0);

	//declare an array of ones for assignment of the Morphology Array
	MorphologyArray = vector< vector<int> >(NZNodes,vector<int>(NXNodes,1));
	ResistanceArray = vector< vector<double> >(NZNodes,vector<double>(NXNodes,RockResistance));

	if (InitialGradient != 0)
	{
		for (int i=0; i<NZNodes; ++i)
		{
			Xz[i] = (Z[i]-Z[NZNodes-1])/InitialGradient;
			for (int j=0;j<NXNodes;++j)
			{
				// need an Epsilon on this check?
				if (j == NXNodes-1)
				{
					Zx[j-1] = Z[i];
				}
				else if (X[j] < Xz[i]+0.00001)
				{
					MorphologyArray[i][j]=0;
					ResistanceArray[i][j]=0;
				}
				else 
				{
					Zx[j-1] = Z[i];
					break;
				}
			}
		}
	}

	//default time interval
	Time = 0;
	TimeInterval = 1;
	EndTime = 10000;

	//default print conditions to print every timestep
	PrintTime = 0;
	PrintInterval = TimeInterval;

	//Set sea level to zero to begin with, and the ind, this will get updated later
	SeaLevelRise = 0;
	SeaLevel = -99;
	SeaLevelInd = 0;

}


void RPM::InitialiseTides(double TideRange)
{
	//setup tides

	//declare temporary variables
	double Total = 0;
	double Max = 0;

	TidalRange = TideRange;
	UpdateMorphology();

	// Make erosion shape function based on tidal duration
	NTideValues = (int)(TidalRange/dZ)+1;
	vector<double> EmptyTideVec(NTideValues,TidalRange/2.);
	ErosionShapeFunction = EmptyTideVec;

	// Loop over tidal range and assign weights
	// If narrow tidal range just use two separate sin functions
	// If wide use two overlapping tidal functions
	if (NTideValues < 20)
	{
		for (int i=0; i<NTideValues; ++i)
		{
			ErosionShapeFunction[i] = sin(i*dZ*M_PI/(0.5*TidalRange));
			if (i == (NTideValues-1)/2.) ErosionShapeFunction[i] += ErosionShapeFunction[i-1];
			if (i*dZ>0.5*TidalRange) ErosionShapeFunction[i] *= -1;

			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
	}
	else
	{
		for (int i=0; i<0.55*NTideValues; ++i)
		{
			ErosionShapeFunction[i] += sin(i*dZ*M_PI/(0.55*TidalRange));
			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
		for (int i=(int)(round(0.45*NTideValues)); i<NTideValues; ++i)
		{
			ErosionShapeFunction[i] += sin((i*dZ-0.45*TidalRange)*M_PI/(0.55*TidalRange));
			Total += ErosionShapeFunction[i];
			if (ErosionShapeFunction[i] > Max) Max = ErosionShapeFunction[i];
		}
	}

	//Normalise values to total
	//Or should I be normalising to Max Value!?
	for (int i=0; i<NTideValues; ++i) ErosionShapeFunction[i] /= Total;
	//for (int i=0; i<NTideValues; ++i) ErosionShapeFunction[i] /= Max;

	//Initialise weathering shape function
	InitialiseWeathering();

	//Initialize BreakingWaveDist and BreakingWaveConst_new
	BreakingWaveDist_New = vector<double>(NTideValues,0.0);
    BreakingWaveConst_New = vector<double>(NTideValues,0.0);
}


void RPM::InitialiseGeology(double CliffElevationNew, double CliffFailureDepthNew, double RockResistanceNew, double WeatheringConstNew, double SubtidalEfficacy)
{
	/* Function to set the cliff height, failure depth, rock resistance and
		weathering rate constant */

	CliffElevation = CliffElevationNew;
	CliffFailureDepth = CliffFailureDepthNew;
	RockResistance = RockResistanceNew;
	WeatheringConst = WeatheringConstNew;
	MinWeatheringEfficacy = SubtidalEfficacy;

	//Loop across the resistance array and reset all values
	for (int i=0;i<NZNodes; ++i)
	{
		for (int j=0; j<NXNodes; ++j)
		{
			ResistanceArray[i][j] = RockResistance;
		}
	}

	InitialiseWeathering();
}

void RPM::InitialiseWeathering()
{
	/* Weathering Efficacy Function guided by Trenhaile and Kanayay (2005)
	using a log-normal distribution across the tidal range */

	// declare control parameters for distribution
	//double sigma = 0.5;sv
	//double Theta = 0;
	MaxWeatheringEfficacy = 1.;
	
	/* need a check here the MinWeatheringRate has been set. InitialiseGeology must have been run before InitialiseWeathering/InitialiseTides */

	//MinWeatheringEfficacy = SubtidalEfficacy;

	// This m value is tailored to cause a distribution peak at 3/4 of the tidal range
	// as in Matsumoto et al. (2016)
	//double m = 1.1665;

	// Make weathering shape function based on tidal duration
	NTideValues = (int)(TidalRange/dZ)+1;
	vector<double> EmptyTideVec(NTideValues,0);
	WeatheringEfficacy = EmptyTideVec;

	//Create a vector ranging from 0 to 10 with NTideValues
	vector<double> LogNormalDistX(NTideValues,0);
	
	for (int i=0; i<NTideValues; ++i) LogNormalDistX[i] = 10.*i/(NTideValues-1);

	//Create weathering efficacy shape function
	for (int i=1; i<NTideValues ;++i)
	{
		if (i<((NTideValues-1)/4)) WeatheringEfficacy[i] = fastexp(-(pow(i-(NTideValues/4.),2.)/(NTideValues/2.)));
        else WeatheringEfficacy[i] = (MaxWeatheringEfficacy-MinWeatheringEfficacy)*fastexp(-(pow(i-(NTideValues/4.),2.))/(NTideValues*NTideValues/10.))+MinWeatheringEfficacy;
	}
}

void RPM::InitialiseWaves(double WaveHeight_Mean, double WaveHeight_StD, double WavePeriod_Mean, double WavePeriod_StD)
{
  /* intialise waves as a single wave */

  MeanWavePeriod = WavePeriod_Mean;
  StdWavePeriod = WavePeriod_StD;
  MeanWaveHeight = WaveHeight_Mean;
  StdWaveHeight = WaveHeight_StD;
}

void RPM::InitialiseWavePressure_Rectangle(double WaveHeight)
{
    /* initialise wave pressure */
    int NWPValues;
    NWPValues = (int)((WaveHeight/dZ)+1.0);
    vector<double> EmptyTideVec(NWPValues,dZ/WaveHeight);
	//vector<double> EmptyTideVec(NWPValues,0.1);
	WavePressure = EmptyTideVec;

}

void RPM::InitialiseWavePressure_Triangle(double WaveHeight)
{
    /* initialise wave pressure */
    int NWPValues,ccc;
    double sum=0.;
    //double sum1=0.;
    NWPValues = (int)((WaveHeight/dZ)+1.0);
    vector<double> EmptyTideVec1(NWPValues,dZ/WaveHeight);
	StandingWavePressure = EmptyTideVec1;

    vector<double> EmptyTideVec2(NWPValues,0.);
	BreakingWavePressure = EmptyTideVec2;
	BrokenWavePressure = EmptyTideVec2;


	if ( NWPValues%2 == 0)
    {
        for (int i=0; i<NWPValues/2; ++i){
            BreakingWavePressure[i] = 2.*i/NWPValues;
            BrokenWavePressure[i] = 2.*i/NWPValues;
            sum = sum + BreakingWavePressure[i];
        }
        ccc=0;
        for (int i=NWPValues/2; i<NWPValues; ++i){
            BreakingWavePressure[i] = BreakingWavePressure[NWPValues/2-1-ccc];
            BrokenWavePressure[i] = BrokenWavePressure[NWPValues/2-1-ccc];
            sum = sum + BreakingWavePressure[i];
            ccc=ccc+1;
        }
    }
    else
    {
        ccc=0;
		for (int i=0; i<NWPValues; ++i)
		{
            if (i<NWPValues/2)
			{
				BreakingWavePressure[i] = 2.*i/NWPValues;
				BrokenWavePressure[i] = 2.*i/NWPValues;
				sum = sum + BreakingWavePressure[i];
			}
			else
			{
				BreakingWavePressure[i] = BreakingWavePressure[NWPValues/2-ccc];
				BrokenWavePressure[i] = BrokenWavePressure[NWPValues/2-ccc];
				sum = sum + BreakingWavePressure[i];
				ccc=ccc+1;
			}
        }
    }
    for (int i=0; i<NWPValues; ++i){
        BreakingWavePressure[i] = BreakingWavePressure[i]/sum;
        BrokenWavePressure[i] = BrokenWavePressure[i]/sum;
    }
}


void RPM::GetWave()
{
	//declare temp variables
	//double OffshoreWaveHeight;
	//double rand1, rand2;

	//Get two random numbers and generate wave data
	//rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	//WavePeriod = MeanWavePeriod + StdWavePeriod*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	//rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	//OffshoreWaveHeight = MeanWaveHeight + StdWaveHeight*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));

	//Breaking Wave Height calculated following Komar and Gaughan (1972)
	//BreakingWaveHeight = 0.39*pow(g,0.2)*pow(WavePeriod,0.4)*pow(OffshoreWaveHeight,2.4);
	WaveHeight = MeanWaveHeight;
	BreakingWaveHeight = MeanWaveHeight;
	
	//Water depth of breaking wave
	BreakingWaveWaterDepth = BreakingWaveHeight/0.78;

	//Get Wave Pressure Distribution as Uniform
	PressureDistMaxInd = (int)(0.5*BreakingWaveHeight/dZ);
	PressureDistMinInd = (int)(-0.5*BreakingWaveHeight/dZ);
	
	// test the initialise wave pressure functions
	InitialiseWavePressure_Rectangle(WaveHeight);
	InitialiseWavePressure_Triangle(WaveHeight);
}

void RPM::InitialiseSeaLevel(double SLR)
{
	SeaLevelRise = SLR;
	// dont let sea level rise be zero, just make it tinee!
	if (SeaLevelRise == 0) SeaLevelRise = 0.0000000001;
	SLR_sum = 0.;
	SeaLevel = 0;
}

void RPM::UpdateSeaLevel()
{

	/*Update sea level based on a constant sea level rise rate*/

	//In case sea level rate is less than cm/year
	if ( abs(SeaLevelRise) < 0.1 )
	{
		SLR_sum = SLR_sum + SeaLevelRise;
		if( SLR_sum >= 0.1)
		{
			//SeaLevel += SeaLevelRise*TimeInterval;
			SeaLevel = SeaLevel + (round(SLR_sum*10))*0.1;
			SLR_sum = 0.;
		}
	}
	//In case sea level rate is less than cm/year
	else
	{
		//SeaLevel += SeaLevelRise*TimeInterval;
		SeaLevel = SeaLevel + (round(SeaLevelRise*10))*0.1;
	}

}

// UpdateSeaLevel_v1 written by Hiro
// Sea level is shifted every 10cm
void RPM::UpdateSeaLevel_v1(double InputSeaLevel)
{
    if ( InputSeaLevel > 0 ){
	/*Update sea level based on a constant sea level rise rate*/
	//In case sea level rate is less than cm/year
	if ( abs(SeaLevelRise) < 0.1 ){
		SLR_sum = SLR_sum + InputSeaLevel;
		if( SLR_sum >= 0.1){
			//SeaLevel += SeaLevelRise*TimeInterval;
			SeaLevel += (round(SLR_sum*10))*0.1;
			SLR_sum = 0.;
		}
	}
	//In case sea level rate is less than cm/year
	else{
		//SeaLevel += SeaLevelRise*TimeInterval;
		SeaLevel += (round(InputSeaLevel*10))*0.1;
	}
    }
}


void RPM::UpdateSeaLevel(double InputSeaLevel)
{
	/*Update sea level based on a new mean sea level InputSeaLevel*/
	//First catch large difference in SeaLevel
	if (fabs(InputSeaLevel-SeaLevel) > 1)
	{
		for (int i=0; i<NZNodes; ++i)
		{
			if (InputSeaLevel > Z[i])
			{
				SeaLevel = Z[i-1];
				SeaLevelInd = i-1;
				break;
			}
		}
	}
	else
	{
		for (int i=SeaLevelInd-11; i<SeaLevelInd+11; ++i)
		{
			if (Z[i] < InputSeaLevel)
			{
				SeaLevel = Z[i-1];
				SeaLevelInd = i-1;
				break;
			}
		}
	}
}

// Only bigger than 10 cm.
void RPM::TectonicUplift(double UpliftAmplitude)
{
    printf("\nUplift Event\n");
	printf("Elevations changed by %1.1f\n",UpliftAmplitude);

	double Uplifted=0; 

	Uplifted=0;
	while (Uplifted < UpliftAmplitude+0.00001)
	{
		MorphologyArray.push_back(MorphologyArray[NZNodes-1]);
		MorphologyArray.erase(MorphologyArray.begin());
		Uplifted += dZ;
	}
	UpdateMorphology();
}

void RPM::CalculateBackwearing()
{
	//Declare temporary variables
	double WaveForce, SurfZoneBottomZ; //, SurfZoneBottomX;
	double dZ = 1.; 
	double dX = 1.;
	double H = sqrt(2.);
	int WaveType;
	BreakingPointZInd = 0;

	//Reset backwear vector
	vector<double> ZZeros(NZNodes,0);
	Bw_Erosion = ZZeros;

	//Loop across all intertidal elevations
	for (int i=MaxTideZInd+1; i<MinTideZInd; ++i)
	{
		//Estimate horizontal breaking point
		//Elevation of breaker point
		SurfZoneBottomZ = Z[i]-BreakingWaveWaterDepth;
		for (int ii=i; ii<NZNodes; ++ii)
		{
			if (Z[ii] < SurfZoneBottomZ)
			{
				BreakingPointZInd = ii;
				break;
			}
		}

		//If waves are breaking at the seaward edge, find the top of the seaward edge
		for (int ii=BreakingPointZInd; ii>i; --ii)
		{
			if (Xz[ii] == Xz[BreakingPointZInd]) BreakingPointZInd = ii;
			else break;
		}

		//Find the location where the broken wave starts
		//BreakingPointX = Xz[BreakingPointZInd];
		
		//Determine Surfzone Gradient
		//Get surf zone mean platform gradient
		//This is critical!
		if ((Xz[MaxXZInd] != Xz[BreakingPointZInd]))
		{
			//Find gradient of local slope
        	for (int ii=BreakingPointZInd+1; ii<NZNodes; ++ii)
        	{
				if ( Xz[BreakingPointZInd]-Xz[ii] > WaveHeight )
            	{
					dX = Xz[BreakingPointZInd]-Xz[ii];
					dZ = Z[BreakingPointZInd]-Z[ii];
					H = sqrt(dX*dX+dZ*dZ);
					SurfZoneGradient = abs(dZ/dX);
					break;
				}
			}
		}
		else 
		{
			SurfZoneGradient = 1.;
			dZ = 1.; 
			dX = 1.;
			H = sqrt(2.);
		}

		//Limit SurfZoneGradient to 45 degrees!
		//if (SurfZoneGradient > 1.) SurfZoneGradient = 1.;

		// Get breaking wave conditions and constants
		BreakingWaveDist = 10.*WaveHeight*(dX/H);
		BreakingWaveConst_New[i-MaxTideZInd] = BreakingWaveConst*(dZ/H);
		if (BreakingWaveConst_New[i-MaxTideZInd] < BrokenWaveConst) BreakingWaveConst_New[i-MaxTideZInd] = BrokenWaveConst;

		//Set Wave Type
		if (Xz[i] == 0) WaveType = 1;
		else if (Xz[i]-Xz[BreakingPointZInd]<=0) WaveType = 1;
		else if ((Xz[i]-Xz[BreakingPointZInd])<BreakingWaveDist) WaveType = 2;
		else WaveType = 3;

		//Determine Backwear erosion
		//Standing wave
		if (WaveType == 1)
		{
			//Loop across pressure distribution function, currently a const
			for (int ii=i+PressureDistMinInd; ii<=i+PressureDistMaxInd; ++ii)
			{
				// Calculate wave force and update backwear at each elevation
				WaveForce = StandingWaveConst*g*rho_w*WaveHeight*ErosionShapeFunction[i-MaxTideZInd];
				Bw_Erosion[ii] += WaveForce;
			}
		}
		//Breaking wave
		//For a breaking wave, first deal with backwear for the standing wave part,
		// then the breaking part
		else if (WaveType == 2 || WaveType == 3)
		{
			//Loop across pressure distribution function
			//This may have some problems!
			for (int ii=i+PressureDistMinInd; ii<=i+PressureDistMaxInd; ++ii)
			{
				//need to add for condition where changes to broken wave above water level in pressure distribution function
				if (Xz[ii] < Xz[BreakingPointZInd])
				{
					WaveForce = StandingWaveConst*g*rho_w*WaveHeight*ErosionShapeFunction[i-MaxTideZInd];

				}
				else if (Xz[ii] <= (Xz[BreakingPointZInd]+BreakingWaveDist))
				{
					BreakingWaveHeight = WaveHeight*fastexp(-WaveAttenuConst*(Xz[ii]-Xz[BreakingPointZInd]));
					WaveForce = BreakingWaveConst_New[i-MaxTideZInd]*g*rho_w*BreakingWaveHeight*ErosionShapeFunction[i-MaxTideZInd];
				}
				else
				{
					BrokenWaveHeight = WaveHeight*fastexp(-WaveAttenuConst*BreakingWaveDist)*fastexp(-WaveAttenuConst*(Xz[ii]-(Xz[BreakingPointZInd]+BreakingWaveDist)));
					WaveForce = BrokenWaveConst*g*rho_w*BrokenWaveHeight*ErosionShapeFunction[i-MaxTideZInd];
				}
				Bw_Erosion[ii] += WaveForce;
			}
		}
	}
}

void RPM::CalculateDownwearing()
{
	//Declare temporary variables
	double WaveForce = 0;
	double WaterDepth = 0;

	//Reset downwear vector
	vector<double> ZZeros(NZNodes,0);
	Dw_Erosion = ZZeros;

	//Loop across the tidal range water levels
	for (int i=MaxTideZInd; i<=MinTideZInd; ++i)
	{
		//Get wave function needs to calculate a bunch of stuff?
		GetWave();

		//Standing Waves
		if (Xz[i] < Xz[BreakingPointZInd])
		{
			WaveForce = StandingWaveConst*g*rho_w*WaveHeight*ErosionShapeFunction[i-MaxTideZInd];
			DepthDecay = -log(SubmarineDecayConst)/WaveHeight;
		}
		//Breaking Waves
		else if (Xz[i]<(Xz[BreakingPointZInd]+BreakingWaveDist))
		{
			WaveForce = BreakingWaveConst_New[i-MaxTideZInd]*g*rho_w*WaveHeight*ErosionShapeFunction[i-MaxTideZInd]*fastexp(-WaveAttenuConst*(Xz[i]-Xz[BreakingPointZInd]));
			DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*fastexp(-BreakingWaveDecay*(Xz[i]-Xz[BreakingPointZInd])));
		}
		//Broken Waves
		else
		{
			WaveForce = BrokenWaveConst*g*rho_w*WaveHeight*fastexp(-WaveAttenuConst*BreakingWaveDist)*ErosionShapeFunction[i-MaxTideZInd]*fastexp(-WaveAttenuConst*(Xz[i]-(Xz[BreakingPointZInd]+BreakingWaveDist)));
			DepthDecay = -log(SubmarineDecayConst)/(WaveHeight*fastexp(-BrokenWaveDecay*(Xz[i]-(Xz[BreakingPointZInd]+BreakingWaveDist))));
		}
		//Loop from water level down and determine force
		for (int ii=i; ii<i+(3./dZ); ++ii)
		{
			if (ii > NZNodes-1) break; 
			WaterDepth = Z[i]-Z[ii];
			Dw_Erosion[ii] += WaveForce*fastexp(-DepthDecay*WaterDepth);
		}
	}
}

void RPM::SupratidalWeathering()
{
	//add this later
}

void RPM::IntertidalWeathering()
{
	//Declare temporay variables
	double RemainingResistance, WeatheringForce;

	//Reset weathering vector
	vector<double> ZZeros(NZNodes,0);
	Weathering = ZZeros;

	//Loop across the tidal range
	for (int i=MaxTideZInd+1; i<MinTideZInd; ++i)
	{
		//Calculate Weathering
		WeatheringForce = WeatheringConst*WeatheringEfficacy[i-MaxTideZInd];

		//How are we going to get j? i.e. x-position in the array?
		//Need a loop in here moving from bottom to top of tidal range in x-position
		for (int j=0; j<=MaxXXInd; ++j)
		{
			//Check we're at a a surface cell
			if ((MorphologyArray[i][j] == 1) && (j == 0))
			{
				RemainingResistance = ResistanceArray[i][j];
				ResistanceArray[i][j] -= WeatheringForce;

				//If resistance is less than zero then cell is lost to weathering erosion
				//excess weathering force is applied to the block behind
				if (ResistanceArray[i][j] < 0)
				{
					MorphologyArray[i][j] = 0;
					ResistanceArray[i][j] = 0;
					WeatheringForce -= RemainingResistance;
				}
			}
			else if ((MorphologyArray[i][j] == 1) && ((MorphologyArray[i-1][j] == 0) || (MorphologyArray[i][j-1] == 0)))
			{
				RemainingResistance = ResistanceArray[i][j];
				ResistanceArray[i][j] -= WeatheringForce;

				//If resistance is less than zero then cell is lost to weathering erosion
				//excess weathering force is applied to the block behind
				if (ResistanceArray[i][j] < 0)
				{
					MorphologyArray[i][j] = 0;
					ResistanceArray[i][j] = 0;
					WeatheringForce -= RemainingResistance;
				}
			}
		}
	}
}

void RPM::SubtidalWeathering()
{
	//Declare temporay variables
	double RemainingResistance, WeatheringForce;

	//Reset weathering vector
	vector<double> ZZeros(NZNodes,0);
	Weathering = ZZeros;

	//Loop across the subtidal
	for (int i=MinTideZInd; i<NZNodes; ++i)
	{
		//Calculate Weathering
		WeatheringForce = WeatheringConst*MinWeatheringEfficacy;

		//How are we going to get j? i.e. x-position in the array?
		//Need a loop in here moving from bottom to top of tidal range in x-position
		for (int j=0; j<=MaxXXInd; ++j)
		{
			//Check we're at a a surface cell
			if ((MorphologyArray[i][j] == 1) && (j == 0))
			{
				RemainingResistance = ResistanceArray[i][j];
				ResistanceArray[i][j] -= WeatheringForce;

				//If resistance is less than zero then cell is lost to weathering erosion
				//excess weathering force is applied to the block behind
				if (ResistanceArray[i][j] < 0)
				{
					MorphologyArray[i][j] = 0;
					ResistanceArray[i][j] = 0;
					WeatheringForce -= RemainingResistance;
				}
			}
			else if ((MorphologyArray[i][j] == 1) && ((MorphologyArray[i-1][j] == 0) || (MorphologyArray[i][j-1] == 0)))
			{
				RemainingResistance = ResistanceArray[i][j];
				ResistanceArray[i][j] -= WeatheringForce;

				//If resistance is less than zero then cell is lost to weathering erosion
				//excess weathering force is applied to the block behind
				if (ResistanceArray[i][j] < 0)
				{
					MorphologyArray[i][j] = 0;
					ResistanceArray[i][j] = 0;
					WeatheringForce -= RemainingResistance;
				}
			}
		}
	}
}

void RPM::ErodeBackwearing()
{
	//Loop over all wet cells
	int j=0;
	double RemainingForce;

	for (int i=MinTideZInd+PressureDistMaxInd; i>=MaxTideZInd+PressureDistMinInd; --i)
	{
		//Find j ind somehow
		//loop across the active shoreface
		j=0;
		while (j < MaxXXInd)
		{
			if (MorphologyArray[i][j] == 1) break;
			++j;
		}

		//Check Backwear Force vs Resistance
		RemainingForce = Bw_Erosion[i];

		while (RemainingForce > ResistanceArray[i][j])
		{
			//Update remaining force
			RemainingForce -= ResistanceArray[i][j];

			//For now assume that only one block can be removed at a time
			//Hiro has code that allows multiple blocks to be removed
			//This will be critical for wave dominated conditions?
			MorphologyArray[i][j] = 0;
			ResistanceArray[i][j] = 0;

			//iterate landward
			++j;
			if (j > MaxXXInd) MaxXXInd = j;

			//May also need soemthing to move ix_max landward by 1
			//If we run out of land dynamically add more columns
			if (j == NXNodes-2)
			{
				//Grow them
				X.push_back(NXNodes*dX);
				Zx.push_back(Zx[NXNodes-1]);
				ZInd.push_back(ZInd[NXNodes-1]);
				for (int i=0; i<NZNodes; ++i)
				{
					MorphologyArray[i].push_back(MorphologyArray[i][NXNodes-1]);
					ResistanceArray[i].push_back(ResistanceArray[i][NXNodes-1]);
				}
				++NXNodes;
			}
		}
	}
}

void RPM::ErodeDownwearing()
{
	//is there any reason why this needs to be a separate function?
	//Loop over all cells that get wet
	for (int j=0; j<MaxTideXInd; ++j)
	{
		// loop over the tidal range?
		for (int i=MaxTideZInd; i<=MinTideZInd; ++i)
		{
			// Check Downwear Force vs Resistance
			if (Dw_Erosion[j] > ResistanceArray[i][j])
			{
				//For now assume that only one block can be removed at a time
				//Hiro has code that allows multiple blocks to be removed
				//I doubt this happens very often with downwear
				MorphologyArray[i][j] = 0;
				ResistanceArray[i][j] = 0;

				//Hiro then has some code to count the number of blocks removed
				//but not clear why this is needed
			}
		}
	}
}

void RPM::DestroyOffshore()
{
	// Function to turn on flag to destroy offshore cells in arrays and vectors
	// MDH, Jan 2020

	// find lowest point to maintain in simulation.
	double MaxWaterDepth = 0.5*TidalRange+3.*MeanWaveHeight;

	// find new min elevation in vertical
	for (int i=0; i<NZNodes; ++i)
	{
		if (Z[i] < SeaLevel-MaxWaterDepth)
		{
			OffshoreZInd = i-1;
			break;
		}
	}

	// find in horizontal
	for (int j=0; j<NXNodes; ++j)
	{
		if (MorphologyArray[OffshoreZInd][j] == 1)
		{
			OffshoreXInd = j;
			break;
		}
	}
	
	//Destroy the left hand side of X direction vectors and arrays dynamically
	//Destroy them first in horizontal
	X.erase(X.begin(),X.begin()+OffshoreXInd);
	  		
	for (int i=0; i<NZNodes; ++i)
	{
		MorphologyArray[i].erase(MorphologyArray[i].begin(),MorphologyArray[i].begin()+OffshoreXInd);
		ResistanceArray[i].erase(ResistanceArray[i].begin(),ResistanceArray[i].begin()+OffshoreXInd);
	}
	
	// update size value
	NXNodes -= OffshoreXInd;
	
	// then destroy in the vertical
	MorphologyArray.erase(MorphologyArray.begin()+OffshoreZInd,MorphologyArray.end());
	ResistanceArray.erase(ResistanceArray.begin()+OffshoreZInd,ResistanceArray.end());
	Z.erase(Z.begin()+OffshoreZInd,Z.end());
	
	NZNodes -= (NZNodes-OffshoreZInd);
	
	//declare arrays of zeros to initalise various other vectors
	vector<double> ZZeros(NZNodes,0);
	Bw_Erosion = ZZeros;
	Dw_Erosion = ZZeros;
	Weathering = ZZeros;
	
	//Indices trackers
	XInd = vector<int>(NZNodes,0);
	ZInd = vector<int>(NXNodes,0);
	
	// rebuild other vectors as appropriate
	UpdateMorphology();
	
}

void RPM::UpdateMorphology()
{
	/*
	Retrieves the morphology of the shore platform and cliff from the resistance and
	morphology arrays. Vectors and arrays are grown dynamically as the coast retreats.
	Vectors and arrays can optionally be shrunk dynamically at the offshore boundary. 
	This is useful for model simulations with continuous RSLR where the offshore domain
	becomes inactive.

	MDH

	*/

	// Find Sea Level in vertical
	// Only need to do this once if sea level isnt changing
	// But need to do it if tectonics are in play
	if ((SeaLevelInd == 0) || (SeaLevelRise != 0))
	{
		for (int i=0; i<NZNodes; ++i)
		{
			if (Z[i] < SeaLevel)
			{
				SeaLevelInd = i-1;
				break;
			}
		}		
	}

	//Loop across all intertidal elevations
	MinTideZInd = (int)round(SeaLevelInd+0.5*TidalRange/dZ);
	MaxTideZInd = (int)round(SeaLevelInd-0.5*TidalRange/dZ);

	//Determine indices in X-direction
	bool LowTideFlag = false;
	bool HighTideFlag = false;
	
	for (int j=0; j<NXNodes; ++j)
	{
		if ((MorphologyArray[MaxTideZInd][j] == 1) && (HighTideFlag == false))
		{
			MaxTideXInd = j;
			HighTideFlag = true;
		}
		if ((MorphologyArray[MinTideZInd][j] == 1) && (LowTideFlag == false))
		{
			 MinTideXInd = j;
			 LowTideFlag = true;
		}
		if ((HighTideFlag == 1) && (LowTideFlag == 1)) break;
	}

	//Populate vector of X values in Z
	MaxXXInd = 0;
	MaxXZInd = 0;
	for (int i=0; i<NZNodes; ++i)
	{
		for (int j=XInd[i]; j<NXNodes; ++j)
		{
			if (MorphologyArray[i][j] == 1)
			{
				Xz[i] = X[j];
				XInd[i] = j;

				if ((i>MaxTideZInd) && (i < MinTideZInd) && (XInd[i] >= MaxXXInd))
				{
					MaxXXInd = XInd[i];
					MaxXZInd = i;
				}
				break;
			}
		}
	}

	//Grow the X direction arrays dynamically shoreward as required
	if (MaxXXInd > NXNodes-10)
	{
		//Grow them
		X.push_back(NXNodes*dX);
		Zx.push_back(Zx[NXNodes-1]);
		ZInd.push_back(ZInd[NXNodes-1]);
		
		for (int i=0; i<NZNodes; ++i)
		{
			MorphologyArray[i].push_back(MorphologyArray[i][NXNodes-1]);
			ResistanceArray[i].push_back(ResistanceArray[i][NXNodes-1]);
		}
		++NXNodes;
	}

	//Populate vector of Z values in X
	for (int j=0; j<NXNodes; ++j)
	{
		for (int i=ZInd[j]; i<NZNodes; ++i)
		{
			if (MorphologyArray[i][j] == 1)
			{
				Zx[j] = Z[i];
				ZInd[j] = i;
				break;
			}
		}
	}
}

void RPM::MassFailure()
{
	//simple implementation for now, talk to Hiro about this
	//Cliff position taken from Highest elevation. <- this could be the problem

	//Find X position of notch and cliff
	double XMax = 0;
	int XMaxZInd = 0;

	// find most landward point in intertidal zone
	for (int i=MinTideZInd+PressureDistMaxInd; i>=MaxTideZInd+PressureDistMinInd; --i)
	{
		if (Xz[i] > XMax)
		{
			XMax = Xz[i];
			XMaxZInd = i;
		}
	}

	// find most seaward point above this notch
	double XMin = XMax;
	for (int i=XMaxZInd; i>0; --i)
	{
		if (Xz[i] < XMin)
		{
			XMin = Xz[i];
		}
		else if (Xz[i] > XMax) break;
	}


	//if big enough then delete
	if ((XMax-XMin) > CliffFailureDepth)
	{
		for (int i=0; i<=XMaxZInd; ++i)
		{
			for (int j=0; j<XMax/dX; ++j)
			{
				MorphologyArray[i][j] = 0;
				ResistanceArray[i][j] = 0;
			}
		}
	}

	// check rest of the intertidal for overhangs and delete to make steps as appropriate
	XMax = 0;
	for (int i=MinTideZInd+PressureDistMaxInd; i>XMaxZInd; --i)
	{
		if (Xz[i] > XMax)
		{
			XMax = Xz[i];
		}
		else if (Xz[i] < XMax)
		{
			Xz[i] = XMax;
			for (int j=0; j<XMax/dX; ++j)
			{
				MorphologyArray[i][j] = 0;
				ResistanceArray[i][j] = 0;
			}
		}
	}
}

void RPM::RunModel(Parameters Params, RockyCoastCRN PlatformCRN)
{
	/*
		Function to evolve the coastal profile through time following
		Matsumoto et al. 2016

		This is the main model loop. Can be executed here or from a driver file

		Martin Hurst 30/3/2017
	*/

 	//Loop through time
	//Loop through time
	while (Time >= Params.EndTime)
	{
		//Do an earthquake?
		if (Params.Earthquakes && Time < UpliftTime)
		{
			TectonicUplift(Params.UpliftMagnitude);
			UpliftTime -= Params.UpliftFrequency;
			
			//Update the Morphology 
			UpdateMorphology();
		}		
		
		//Update Sea Level
		InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);
		UpdateSeaLevel(InstantSeaLevel);

		//Get the wave conditions
		GetWave();

		//Calculate forces acting on the platform
		CalculateBackwearing();
		CalculateDownwearing();

		//Do erosion
		ErodeBackwearing();
		ErodeDownwearing();

		//Implement Weathering
		IntertidalWeathering();
		SubtidalWeathering();
		
		//Check for Mass Failure
		MassFailure();
		
		//Update the Morphology 
		UpdateMorphology();

        //Update the morphology inside RockyCoastCRN
		if (Params.CRN_Predictions) 
		{
			PlatformCRN.UpdateMorphology(PlatformModel);
			PlatformCRN.UpdateCRNs();
		}
        	
		//print?
		if (Time <= PrintTime)
		{
			WriteProfile(Params.ProfileOutFilename, Time);
			if (Params.CRN_Predictions) PlatformCRN.WriteCRNProfile(Params.ConcentrationsOutFilename, Time);
			PrintTime -= Params.PrintInterval;
		}
		
		//update time
		Time -= Params.TimeStep;
	}
}

void RPM::WriteProfile(string OutputFileName, double Time, bool Print2Screen)
{
  /* Writes a RPM object X coordinates to file, each value spans dZ in elevation
		File format is

		StartZ dZ
		Time | SeaLevel | X[0] | X[1] | X[2] =====> X[NoNodes] */


	//Print to screen ?
	if (Print2Screen)
	{
		cout.flush();
		cout << "RPM: Writing output at Time " << setprecision(2) << fixed << Time << " years\r";
	}

	//test if output file already exists
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteCoastFile;
	if (FileExists == 0)
	{
		WriteCoastFile.open(OutputFileName.c_str());
		if (WriteCoastFile.is_open()) WriteCoastFile << MaxElevation << " " << MinElevation << " " << dZ << endl;
	}
	WriteCoastFile.close();

	//open output filestream again to  coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write X
		WriteCoastFile << setprecision(4) << Time << " " << setprecision(4) << SeaLevel;
		for (int i=0; i<NZNodes; ++i) WriteCoastFile << setprecision(10) << " " << Xz[i];
		WriteCoastFile << endl;
	}
	else
	{
		//report errors
		cout << "RPM.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
	WriteCoastFile.close();
}

void  RPM::WriteResistanceArray(string OutputFileName, double Time)
{
  /* Writes a RPM object Resistance matrix coordinates to file
		File format is

		Time
			X[0][0]     |    X[0][1]    |   X[0][2]     =====>    X[0][NXNodes]
			X[1][0]     |    X[1][1]    |   X[1][2]     =====>    X[0][NXNodes]
			X[2][0]     |    X[2][1]    |   X[2][2]     =====>    X[0][NXNodes]
		      ||               ||             ||                      ||
		      \/               \/             \/                      \/
		X[NZNodes][0]  | X[NZNodes][1] | X[NZNodes][2] =====> X[NZNodes][NXNodes] */

  	//open the output filestream and write headers
	ofstream WriteFile;
	WriteFile.open(OutputFileName.c_str());
	WriteFile << Time << " " << dZ << " " << dX << endl;

	//Check if file exists if not open a new one and write headers
	if (WriteFile.is_open())
	{
		//write Resistance
		for (int i=0; i<NZNodes; ++i)
		{
			for (int j=0;j<NXNodes; ++j)
			{
				WriteFile << setprecision(5) << ResistanceArray[i][j] << " ";
			}
			WriteFile << endl;
		}
	}
	else
	{
		//report errors
		cout << "RPM.WriteResistance: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}

void  RPM::WriteMorphologyArray(string OutputFileName, double Time)
{
  /* Writes a RPM object Resistance matrix coordinates to file
		File format is

		Time
			X[0][0]     |    X[0][1]    |   X[0][2]     =====>    X[0][NXNodes]
			X[1][0]     |    X[1][1]    |   X[1][2]     =====>    X[0][NXNodes]
			X[2][0]     |    X[2][1]    |   X[2][2]     =====>    X[0][NXNodes]
		      ||               ||             ||                      ||
		      \/               \/             \/                      \/
		X[NZNodes][0]  | X[NZNodes][1] | X[NZNodes][2] =====> X[NZNodes][NXNodes] */

  	//open the output filestream and write headers
	ofstream WriteFile;
	WriteFile.open(OutputFileName.c_str());
	WriteFile << Time << " " << dZ << " " << dX << endl;

	//Check if file exists if not open a new one and write headers
	if (WriteFile.is_open())
	{
		//write Resistance
		for (int i=0; i<NZNodes; ++i)
		{
			for (int j=0;j<NXNodes; ++j)
			{
				WriteFile << setprecision(5) << MorphologyArray[i][j] << " ";
			}
			WriteFile << endl;
		}
	}
	else
	{
		//report errors
		cout << "RPM.WriteMorphology: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}


// void  RPM::WriteLocalangleArray(string OutputFileName, double Time)
// {
//   /* Writes a RPM object Resistance matrix coordinates to file
// 		File format is

// 		Time
// 			X[0][0]     |    X[0][1]    |   X[0][2]     =====>    X[0][NXNodes]
// 			X[1][0]     |    X[1][1]    |   X[1][2]     =====>    X[0][NXNodes]
// 			X[2][0]     |    X[2][1]    |   X[2][2]     =====>    X[0][NXNodes]
// 		      ||               ||             ||                      ||
// 		      \/               \/             \/                      \/
// 		X[NZNodes][0]  | X[NZNodes][1] | X[NZNodes][2] =====> X[NZNodes][NXNodes] */

//   	//open the output filestream and write headers
// 	ofstream WriteFile;
// 	WriteFile.open(OutputFileName.c_str());
// 	WriteFile << Time << " " << dZ << " " << dX << endl;

// 	//Check if file exists if not open a new one and write headers
// 	if (WriteFile.is_open())
// 	{
// 		//write Resistance
// 		for (int i=0; i<NZNodes; ++i)
// 		{
// 			for (int j=0;j<NXNodes; ++j)
// 			{
// 				WriteFile << setprecision(5) << LocalangleArray[i][j] << " ";
// 			}
// 			WriteFile << endl;
// 		}
// 	}
// 	else
// 	{
// 		//report errors
// 		cout << "RPM.LocalAngle: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
// 		exit(EXIT_FAILURE);
// 	}
// }



#endif
