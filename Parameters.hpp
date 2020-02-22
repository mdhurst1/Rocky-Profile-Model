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
 		map<string,int> int_DefaultParams;
 		map<string,float> float_DefaultParams;
  		map<string,bool> bool_DefaultParams;
  		map<string,string> string_DefaultParams;

		void Parameters::Initialise();	
		void Parameters::Initialise(string ParameterFilename);
		void Parameters::SetDefaultValues();
	
	public:
		Parameters()
		{
			Initialise();
		}

		Parameters(string ParameterFilename)
		{
			Initialise(string ParameterFilename);
		}
}
