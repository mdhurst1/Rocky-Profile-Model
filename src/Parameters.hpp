/*------------------------------------------------------------------------

	Parameters.hpp

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

#include <map>

#ifndef Params_HPP
#define Params_HPP

class Parameters
{
	private:

		// maps for setting default parameters
 		// map<string,int> int_Params; there currently aren't any
 		map<string,float> float_Params;
  		map<string,bool> bool_Params;
  		map<string,string> string_Params;

		void Initialise();	
		void Initialise(string Folder, string ParameterFilename);
		void SetDefaultValues();
		void ParseValuesFromFile();
		void WriteToFile();
		void PrintToScreen();
	
	public:

		// resolution params
		float dX, dZ;

		// Cosmogenic Isotopes
		bool CRN_Predictions, Berylium, Carbon, Aluminium, ReadSeaLevelFromFile, Earthquakes;
		vector<int> Nuclides;

		// Hydrodynamics
		float SeaLevelRise, TidalRange, TidalPeriod, WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD, StandingWaveCoef, BreakingWaveCoef, BrokenWaveCoef, WaveAttenuationConst;

		// geology
		float InitialGradient, CliffElevation, MaxElevation, MinElevation, Resistance, 
				WeatheringRate, SubtidalEfficacy, CliffFailureDepth;

		// tectonics
		float UpliftMagnitude, UpliftFrequency;

		// time control
		double StartTime, EndTime, TimeStep, PrintInterval;

		// MCMC Params
		// variable parameters and their distribution metrics for MCMC
		float Resistance_Mean, Resistance_Std, Resistance_Min, Resistance_Max;
		float WeatheringRate_Mean, WeatheringRate_Std, WeatheringRate_Min, WeatheringRate_Max;
		float WaveAttenuation_Mean, WaveAttenuation_Std, WaveAttenuation_Min, WaveAttenuation_Max;
		// weightings
		float TopoWeighting, CRNWeighting;
		
		// input files
		string TopoFilename, CRNFilename;
		string SeaLevelFilename;
		
		// output files
		string Folder, Filename, ProjectName;
		string ProfileOutFilename, ConcentrationsOutFilename;
		string ParameterOutFilename;

		Parameters()
		{
			Initialise();
		}

		Parameters(string Folder, string ParameterFilename)
		{
			Initialise(Folder, ParameterFilename);
		}
};

#endif