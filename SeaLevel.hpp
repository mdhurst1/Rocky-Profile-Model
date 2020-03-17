/*------------------------------------------------------------------------

	SeaLevel.hpp
	
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

#ifndef SeaLevel_HPP
#define SeaLevel_HPP

#include <vector>
#include <cstring>

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

///@brief Main sea level object
class SeaLevel
{
  	private:
  	
  		//double Time;
  		double MeanSeaLevel;
		double SeaLevelRise;
  		int NTimes;

  		vector<double> Times;
  		vector<double> MeanSeaLevels;
  		
  		void Initialise();
  		void Initialise(double SLR, double StartTime, double EndTime, double TimeStep);
  		void Initialise(string SeaLevelDataFile);
  	
  	protected:

	public:
	
		SeaLevel()
		{
			Initialise();
		}
		
		SeaLevel(double SLR, double StartTime, double EndTime, double TimeStep)
		{
			Initialise(SLR);
		}
		
		SeaLevel(string SeaLevelDataFileIn)
		{
			Initialise(SeaLevelDataFileIn);
		}
				
		double get_SeaLevel(double Time);
};
#endif
		

