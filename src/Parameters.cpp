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
    printf("\nParameters.Initialise: Initialised parameters from:\n\tfolder: %s\n\tfile: %s\n\n", 
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
	bool_Params["CRN_Predictions"] = true;
	bool_Params["Berylium"] = true;
	bool_Params["Carbon"] = true;
	bool_Params["Aluminium"] = true;
	
	// Hydrodynamics
	float_Params["SeaLevelRise"] = 0.;
	float_Params["TidalRange"] = 2.;
    float_Params["TidalPeriod"] = 12.42;
	float_Params["WaveHeight_Mean"] = 2.;
    float_Params["WaveHeight_StD"] = 0.;
    float_Params["WavePeriod_Mean"] = 6.;
    float_Params["WavePeriod_StD"] = 0.;
	float_Params["StandingWaveCoef"] = 0.1;
	float_Params["BreakingWaveCoef"] = 10.;
	float_Params["BrokenWaveCoef"] = 0.01;
    float_Params["WaveAttenuationConst"] = 0.01;

	// geology
	float_Params["InitialGradient"] = 1.;
	float_Params["CliffElevation"] = 15.;
    float_Params["MaxElevation"] = 15.;
	float_Params["MinElevation"] = -15.;
	float_Params["Resistance"] = 20.;
	float_Params["WeatheringRate"] = 0.1;
	float_Params["SubtidalEfficacy"] = 0.1;
	float_Params["CliffFailureDepth"] = 0.1;

	// time control
	float_Params["StartTime"] = 8000.;
	float_Params["EndTime"] = 0.;
	float_Params["TimeStep"] = 1.;
	float_Params["PrintInterval"] = 100.;

	// output files
    string_Params["Folder"] = Folder;
    string_Params["ProjectName"] = "RPM_CRN";
    string_Params["SeaLevelFilename"] = "NULL";
	string_Params["ProfileOutFilename"] = string_Params["ProjectName"] + "_ShoreProfile.xz";
	string_Params["ConcentrationsOutFilename"] = string_Params["ProjectName"] + "_Concentrations.xn";
}

void Parameters::ParseValuesFromFile()
{
    // temp variables to hold values read from file
    string Line, Character, ParameterName, Value;
    string Parameter, value, lower, lower_val;
    string bc;
    int ValuePosition = 0;
    bool GotParameter, GotValue;
    
    // file stream to read contents
    ifstream infile;
    infile.open(Filename.c_str());

    // now ingest parameters
    while (infile.good())
    {
        // read the line
        getline(infile,Line);
        GotParameter = false;
        GotValue = false;

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
                else if (!GotValue)
                {
                    ValuePosition = i+1;
                }
                else
                {
                    break;
                }
            }

            else if (Character == "\r" || Character == "\n")
            {
                if (GotParameter) Value = Line.substr(ValuePosition, i-ValuePosition);
                break;
            }
            else if (GotParameter && !GotValue) GotValue = true;
        }

        if (!GotParameter) continue;

        // check if parameter can be found as a key in parameter maps and if so update with value 
        if (bool_Params.find(Parameter) != bool_Params.end())
        {
            bool boolValue;
            istringstream is(Value);
            is >> std::boolalpha >> boolValue;
            bool_Params[Parameter] = boolValue;
        } 
        else if (float_Params.find(Parameter) != float_Params.end()) 
        {
            float_Params[Parameter] = stof(Value);
        }
        // currently no integer parameters
        //else if (int_Params.find(Parameter) != int_Params.end()) 
        //{
        //    int_Params[Parameter] = stoi(Value);
        //}
        else if (string_Params.find(Parameter) != string_Params.end()) 
        {
            string_Params[Parameter] = Value;
        }
        else printf("Parameter %s not found, ignoring value\n", Parameter.c_str());
    }

    // assign all values from variable maps to variables
    
    // Cosmogenic Isotopes
	CRN_Predictions = bool_Params["CRN_Predictions"];
    Berylium = bool_Params["Berylium"];
    Carbon = bool_Params["Carbon"];
    Aluminium = bool_Params["Aluminium"];
    ReadSeaLevelFromFile = bool_Params["ReadSeaLevelFromFile"];
		
    // Hydrodynamics
    SeaLevelRise = float_Params["SeaLevelRise"];
    TidalRange = float_Params["TidalRange"];
    TidalPeriod = float_Params["TidalPeriod"];
    WaveHeight_Mean = float_Params["WaveHeight_Mean"];
    WaveHeight_StD = float_Params["WaveHeight_StD"];
    WavePeriod_Mean = float_Params["WavePeriod_Mean"];
    WavePeriod_StD = float_Params["WavePeriod_StD"];
    StandingWaveCoef = float_Params["StandingWaveCoef"];
    BreakingWaveCoef = float_Params["BreakingWaveCoef"];
    BrokenWaveCoef = float_Params["BrokenWaveCoef"];
    WaveAttenuationConst = float_Params["WaveAttenuationConst"];

    // geology
    InitialGradient = float_Params["InitialGradient"];
    CliffElevation = float_Params["CliffElevation"];
    MaxElevation = float_Params["MaxElevation"];
    MinElevation = float_Params["MinElevation"];
    Resistance = float_Params["Resistance"];
    WeatheringRate = float_Params["WeatheringRate"];
    SubtidalEfficacy = float_Params["SubtidalEfficacy"];
    CliffFailureDepth = float_Params["CliffFailureDepth"];

    // time control
    StartTime = float_Params["StartTime"];
    EndTime = float_Params["EndTime"];
    TimeStep = float_Params["TimeStep"];
    PrintInterval = float_Params["PrintInterval"];

    // output files
    Folder = string_Params["Folder"];
    Filename = string_Params["Filename"];
    ProjectName = string_Params["ProjectName"];
    SeaLevelFilename = string_Params["SeaLevelFilename"];

    ProfileOutFilename = ProjectName + "_ShoreProfile.xz";
	ConcentrationsOutFilename = ProjectName + "_Concentrations.xn";
    ParameterOutFilename = ProjectName + ".params";
    
    // check if reading sea level from file
    if (SeaLevelFilename != "NULL" && SeaLevelFilename != "") ReadSeaLevelFromFile = true;

    WriteToFile();
    PrintToScreen();
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
    
    for (auto it: string_Params) ParamsOut << "\t" << it.first << ": " << it.second << endl;

    for (auto it: bool_Params) ParamsOut << "\t" << it.first << ": " << it.second << endl;
    
    //for (auto it: int_Params) ParamsOut << "\t" << it.first << ": " << it.second << endl;
    
    for (auto it: float_Params) ParamsOut << "\t" << it.first << ": " << it.second << endl;
}

void Parameters::PrintToScreen()
{
    cout << "Parameters: Here are the parameters used by RPM_CRN:" << endl;
    
    for (auto it: string_Params) cout << "\t" << it.first << ": " << it.second << endl;

    for (auto it: bool_Params) cout << "\t" << it.first << ": " << it.second << endl;
    
    //for (auto it: int_Params) cout << "\t" << it.first << ": " << it.second << endl;
    
    for (auto it: float_Params) cout << "\t" << it.first << ": " << it.second << endl;

}

#endif