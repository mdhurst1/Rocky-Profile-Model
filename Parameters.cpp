/*------------------------------------------------------------------------

	Parameters.cpp

	C++ object for managing parameters and parameter files for RPM 
    and RockCoastCRN coupled model 

	Martin D. Hurst, University of Glasgow
	
	Feb 2020

	Copyright (C) 2020, Martin Hurst

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

#ifndef Params_CPP
#define Params_CPP

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
#include "RPM.hpp"
#include "SeaLevel.hpp"
#include "FastExp.hpp"

using namespace std;

void Parameters::Initialise()
{
    /* initialise an empty RPM object as default */
    printf("\nParameters.Initialise: Initialised an empty Parameters object\n");
    printf("\tAll Parameters will be set to default values\n")

    SetDefaultValues();
}

void Parameters::Initialise(string ParameterFilename)
{
    printf("\nParameters.Initialise: Initialised parameters from %s\n", ParameterFilename);

    // initialise with default values
    SetDefaultValues();

    // update to include values parsed form the input file

}

void Parameters::SetDefaultValues()
{
    // Cosmogenic Isotopes
	bool_DefaultParams["CRN_predictions"] = true;
	bool_DefaultParams["Berylium"] = true;
	bool_DefaultParams["Carbon"] = true;
	bool_DefaultParams["Aluminium"] = true;
	
	// Hydrodynamics
	float_DefaultParams["SeaLevelRise"] = 0.001;
	float_DefaultParams["TidalRange"] = 2.;
	float_DefaultParams["WaveHeight"] = 2.;
	float_DefaultParams["StandingWaveCoef"] = 0.1;
	float_DefaultParams["BreakingWaveCoef"] = 10.;
	float_DefaultParams["BrokenWaveCoef"] = 0.01;

	// geology
	float_DefaultParams["InitialGradient"] = 1.;
	float_DefaultParams["CliffHeight"] = 15.;
	float_DefaultParmas["MinElevation"] = -15.;
	float_DefaultParams["Resistance"] = 0.002;
	float_DefaultParams["WeatheringRate"] = 0.0001;
	float_DefaultParams["SubtidalEfficacy"] = 0.1;
	float_DefaultParams["CliffFailureDepth"] = 0.1;

	// time control
	int_DefaultParmas["StartTime"] = -8000;
	int_DefaultParmas["EndTime"] = 0;
	int_DefaultParams["TimeStep"] = 1;
	int_DefaultParmas["PrintInterval"] = 100;

	// output files
	string_DefaultParams["OutputProfileFilename"] = "RPM_ShoreProfile.xz";
	string_DefaultParams["OutputConcentrationFilename"] = "RPM_Concentrations.xn";
}