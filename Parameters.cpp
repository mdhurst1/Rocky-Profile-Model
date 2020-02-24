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
#include <sstream>
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
#include "Parameters.hpp"

using namespace std;

void Parameters::Initialise()
{
    /* initialise an empty RPM object as default */
    printf("\nParameters.Initialise: Initialised an empty Parameters object\n");
    printf("\tAll Parameters will be set to default values\n");

    SetDefaultValues();
}

void Parameters::Initialise(string Folder, string ParameterFilename)
{
    printf("\nParameters.Initialise: Initialised parameters from:\n\tfolder: %s\n\tfile%s\n", 
            Folder.c_str(), ParameterFilename.c_str());

    // set filename
    Folder = Folder;
    Filename = ParameterFilename;

    // initialise with default values
    SetDefaultValues();

    // update to include values parsed form the input file
    ParseValuesFromFile();
}

void Parameters::SetDefaultValues()
{
    // Cosmogenic Isotopes
	bool_Params["CRN_predictions"] = true;
	bool_Params["Berylium"] = true;
	bool_Params["Carbon"] = true;
	bool_Params["Aluminium"] = true;
	
	// Hydrodynamics
	float_Params["SeaLevelRise"] = 0.001;
	float_Params["TidalRange"] = 2.;
	float_Params["WaveHeight"] = 2.;
	float_Params["StandingWaveCoef"] = 0.1;
	float_Params["BreakingWaveCoef"] = 10.;
	float_Params["BrokenWaveCoef"] = 0.01;

	// geology
	float_Params["InitialGradient"] = 1.;
	float_Params["CliffHeight"] = 15.;
	float_Params["MinElevation"] = -15.;
	float_Params["Resistance"] = 0.002;
	float_Params["WeatheringRate"] = 0.0001;
	float_Params["SubtidalEfficacy"] = 0.1;
	float_Params["CliffFailureDepth"] = 0.1;

	// time control
	int_Params["StartTime"] = -8000;
	int_Params["EndTime"] = 0;
	int_Params["TimeStep"] = 1;
	int_Params["PrintInterval"] = 100;

	// output files
	string_Params["OutputProfileFilename"] = "RPM_ShoreProfile.xz";
	string_Params["OutputConcentrationFilename"] = "RPM_Concentrations.xn";
}

void Parameters::ParseValuesFromFile()
{
    // temp variables to hold values read from file
    string Line, Character, ParameterName, Value;
    string Parameter, value, lower, lower_val;
    string bc;
    int ValuePosition;
    bool GotParameter = false;
    
    // file stream to read contents
    ifstream infile;
    infile.open(Filename.c_str());

    // now ingest parameters
    while (infile.good())
    {
        // read the line
        getline(infile,Line);

        for (unsigned int i = 0; i<Line.length(); i++)
        {
            // ignore if commented out
            Character = Line[i];
            if (Character == "#") break;
            
            // find whitespace characters and colon
            else if (Character == ":" || Character == "\t" || Character == " " || Character == ",")
            {
                // split to get parameter name and value and overwrite default
                if (!GotParameter)
                {
                    Parameter = Line.substr(0,i);
                    ValuePosition = i+1;
                    GotParameter = true;
                }
                else
                {
                    ValuePosition = i+1;
                }
            }

            else if (Character == "\n")
            {
                Value = Line.substr(ValuePosition,i);
                break;
            }
        }

        // check if parameter can be found as a key in parameter maps and if so update with value 
        if (bool_Params.find(Parameter) != bool_Params.end())
        {
            bool boolValue;
            istringstream is(Value);
            is >> std::boolalpha >> boolValue;
            bool_Params[Parameter];
        } 
        else if (float_Params.find(Parameter) != float_Params.end()) float_Params[Parameter] = stof(Value);
        else if (int_Params.find(Parameter) != int_Params.end()) int_Params[Parameter] = stoi(Value);
        else if (string_Params.find(Parameter) != string_Params.end()) string_Params[Parameter] = Value;
        else printf("Parameter %s not found, ignoring value\n", Parameter.c_str());
    }

    // check if reading sea level from file
    if (SeaLevelFile != "NULL") ReadSeaLevelFromFile = true;
}

void Parameters::WriteToFile()
{
    ofstream ParamsOut;
    ParamsOut.open(ParameterOutFilename.c_str());

    ParamsOut << "# Here are the parameters used by RPM_CRN:" << endl;
    ParamsOut << "# The folder and file names are: " << endl;
    ParamsOut << "Folder: " << Folder << endl;
    ParamsOut << "ParameterFile: " << Filename << endl;
    ParamsOut << "ProfileOutFilename: " << ProfileOutFilename << endl;
    ParamsOut << "ConcentrationsOutFilename: " << ConcentrationsOutFilename << endl;
    ParamsOut << "# _________________________________________"  << endl;
    ParamsOut << "# Here are the parameters used in your run" << endl;
    
    for (auto it: bool_Params)
    {
        ParamsOut << it.first << ": " << it.second << endl;
    }

    /*
    temp comment out until tested above
    
    for (unsigned int i=0; i<int_Params.size(); ++i)
    {
        ParamsOut << int_Params[i] << endl;
    }

    for (unsigned int i=0; i<float_Params.size(); ++i)
    {
        ParamsOut << float_Params[i] << endl;
    }
    */
}

void Parameters::PrintToScreen()
{
    cout << "# Here are the parameters used by RPM_CRN:" << endl;
    cout << "# The folder and file names are: " << endl;
    cout << "Folder: " << Folder << endl;
    cout << "ParameterFile: " << Filename << endl;
    cout << "ProfileOutFilename: " << ProfileOutFilename << endl;
    cout << "ConcentrationsOutFilename: " << ConcentrationsOutFilename << endl;
    cout << "# _________________________________________"  << endl;
    cout << "# Here are the parameters used in your run" << endl;
    
    for (auto it: bool_Params)
    {
        cout << it.first << ": " << it.second << endl;
    }

    /*
    temp comment

    for (int i=0; i<int(bool_Params); ++i)
    {
        cout << int_Params[i] << endl;
    }

    for (int i=0; i<int(bool_Params); ++i)
    {
        cout << float_Params[i] << endl;
    }
    */
}

#endif