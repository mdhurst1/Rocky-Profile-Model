/*------------------------------------------------------------------------

	SeaLevel.cpp
	
	Sea Level Object
	
	Martin D. Hurst, University of Glasgow, April 2017

	Copyright (C) 2017, Martin Hurst
	
	Developer can be contacted:
	martin.hurst@glasgow.ac.uk

	Martin D. Hurst
	
  
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

#ifndef SeaLevel_CPP
#define SeaLevel_CPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "SeaLevel.hpp"






using namespace std;

void SeaLevel::Initialise()
{
	// Initialising SeaLevel using pretend Milankovitch Cycles
	// Parameterisation of this is a complete fabrication to get
	// the right sort of pattern
	
	// Set the amplitude and frequency of each cyclicity
	double Obliquity_Amp = 5.;
	double Obliquity_Period = 41000.;
	double Obliquity_Offset = (double)rand()/RAND_MAX;
	
	double Precession_Amp = 8.;
	double Precession_Period = 26000.;
	double Precession_Offset = (double)rand()/RAND_MAX;
	
	double Eccentricity_Amp = 90.;
	double Eccentricity_Period = 100000.;
	double Eccentricity_Offset = 10000.; //set to roughly match end of Holocene
	
	
	//Set up Times vector
	// Calculate the relative contributions
	vector<double> Empty(710000,0.);
	Times = Empty;
	MeanSeaLevels = Empty;
	for (int t=0, T=Times.size(); t<T; ++t) 
	{
		Times[t] = -700000.+t;
		// Obliquity signal
		MeanSeaLevels[t] += Obliquity_Amp*cos(2.*M_PI*(Times[t]+Obliquity_Offset*Obliquity_Period)/Obliquity_Period);
		// Precession signal
		MeanSeaLevels[t] += Precession_Amp*cos(2.*M_PI*(Times[t]+Precession_Offset*Precession_Period)/Precession_Period);
		// Eccentricity signal
		MeanSeaLevels[t] += (Eccentricity_Amp/M_PI)*atan(1./tan(M_PI*(Times[t]+Eccentricity_Offset)/Eccentricity_Period));
	}
	// break here and check
	cout << "SeaLevel.Initialise: Sea Level History Created." << endl;
}

void SeaLevel::Initialise(string SeaLevelDataFile)
{
	// Initialising Sea Level from a datafile
	// Use this for SeaLevel output from other models or other RSL records
	double indata;
	char Dummy[32];
	
	//read in RSL data from Bradley et al. (2011) model
	ifstream SeaLevelIn(SeaLevelDataFile.c_str());
	if (!SeaLevelIn)
	{
	  printf("SeaLevel::%s: line %d Sea Level data file \"%s\" doesn't exist\n\n", __func__, __LINE__, SeaLevelDataFile.c_str());
	  printf("Setting Realtive Sea Level to zero");
	}
	else
	{
		SeaLevelIn >> Dummy;
		SeaLevelIn >> Dummy;
		while (!SeaLevelIn.eof())
		{
		  SeaLevelIn >> indata;
		  Times.push_back(indata);
		  SeaLevelIn >> indata;
		  MeanSeaLevels.push_back(indata);
		}
		SeaLevelIn.close();
	}

	cout << "Sea level read in: " << endl;
	
	for (int i=0, N=MeanSeaLevels.size(); i<N; ++i)
	{
		cout << Times[i] << ", " << MeanSeaLevels[i] << endl;
	}
}

void SeaLevel::Initialise(double SLR)
{
	// Initialising Sea Level at zero and set sea level rise rate
	MeanSeaLevel = 0;
	SeaLevelRise = SLR;

	// initialise an empty vector
	int MaxTime = 10000;
	vector<double> Empty(MaxTime,0);
	Times = Empty;
	MeanSeaLevels = Empty;
	
	for (int t=0, T=Times.size(); t<T; ++t)
	{
		// Calculate the time 
		Times[t] = MaxTime-t;

		// sea level as a function of sea level rise rate
		MeanSeaLevels[t] = MeanSeaLevel + SeaLevelRise*(MaxTime-Times[t]);
	}
}

double SeaLevel::get_SeaLevel(double Time)
{
	//interpolate to get value
	double Factor;
	int TimeCondition = 0;
	int ind = 0;
	while (TimeCondition == 0)
	{
		if (Time == Times[ind])
		{
			TimeCondition = 1;
			++ind;
		}
		else if (Time < Times[ind]) ++ind;
		else TimeCondition = 1;
	}
	Factor = (Time-Times[ind-1])/(Times[ind]-Times[ind-1]);
	MeanSeaLevel = MeanSeaLevels[ind-1]+Factor*(MeanSeaLevels[ind]-MeanSeaLevels[ind-1]);

	return MeanSeaLevel;
}
#endif

