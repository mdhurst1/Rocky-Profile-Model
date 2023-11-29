/*------------------------------------------------------------------------

	RPM_dakota_driver_2.cpp
	
	Driver file for running the shore platform model of Matsumoto et al. (2016)
	Updated following improvements by Matsumoto et al. (2018)

    This driver file takes model parameters as command line arguments to enable running
    ensembles on a HPC using TORQUE/PBS qsub command
	
	C++ implementation of Hiro Matsumoto's Shore Platform Model with coupling to Cosmogenic Isotope production by RockCoastCRN/RoBoCoP.

	Matsumoto, H., Dickson, M. E., & Kench, P. S. (2016a)
	An exploratory numerical model of rocky shore profile evolution. 
	Geomorphology http://doi.org/10.1016/j.geomorph.2016.05.017
	
	Matsumoto, H., Dickson, M.E., and Kench, P.S. (2016b)
	Modelling the Development of Varied Shore Profile Geometry on Rocky Coasts.
	Journal of Coastal Research http://dx.doi.org/10.2112/SI75-120.1

	Matsumoto, H., Dickson, M. E., Kench, P. S., (2018)
	Modelling the relative dominance of wave erosion and weathering processes in shore platform development in micro- to mega-tidal settings
	Earth Surface Processes and Landforms  http://dx.doi.org/10.1002/esp.4422
	
    Hurst, M.D., Rood, D.H., Ellis, M.A., Anderson, R.S., and Dornbusch, U. (2016) Recent acceleration in coastal cliff retreat rates on the south coast of Great Britain. Proceedings of the National Academy of Sciences, http://dx.doi.org/10.1073/PNAS.1613044113

	Hurst, M.D., Rood, D.H., and Ellis, M.A. (2017)
	Controls on the distribution of cosmogenic 10 Be across shore platforms
	Earth Surface Dynamics http://dx.doi.org/10.5194/esurf-5-67-2017

	Martin D. Hurst, University of Glasgow
	Hironori Matsumoto, University of Auckland
	Mark Dickson, University of Auckland
	
	July 2019
	
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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <stdio.h>
#include <time.h>
#include "../src/RPM.hpp"
#include "../src/MCMC_RPM.hpp"
#include "../src/RockyCoastCRN.hpp"
#include "../src/SeaLevel.hpp"
#include "../src/FastExp.hpp"

using namespace std;

int main(int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "--------------------------------------------------------------" << endl;
	cout << "|  MCMC RPM CRN                                                |" << endl;
	cout << "|  This program peformas a Metropolis-Hastings Markov Chain    |" << endl;
    cout << "|  to find the best parameter combinations that fit the        |" << endl;
    cout << "|  topography and CRN concentrations measured at a field site. |" << endl;
    cout << "|                                                              |" << endl;
    cout << "|  The model simulates the development of shore platforms      |" << endl;
	cout << "|  following model developed by Matsumoto et al. (2016), and   |" << endl;
	cout << "|  the accumulation of cosmogenic isotopes following Hurst     |" << endl;
	cout << "|  et al. (2016; 2017).                                        |" << endl;
	cout << "|                                                              |" << endl;
	cout << "|  Implemented in C++ by Martin Hurst, University of Glasgow   |" << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << endl;

	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "Error: This program requires two inputs: " << endl;
		cout << " * First a path to the folder where the model will be run" << endl;
		cout << " * The name of the input parameter file (must be in the above folder)" << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_CRN_MCMC.out /ProjectFolder/ MyInputParameterFile.in" << endl;
		cout << "------------------------------------------------------" << endl;
		exit(EXIT_SUCCESS);
	}
    
    // retrieve prameter file and workspace from arguments
    string Folder = argv[1];
	string TempParamFilename = argv[2];
	string InputParamFilename = Folder+TempParamFilename;
		
	// load parameter parser object
  	Parameters Params(Folder,InputParamFilename);
    
    // initiate MCMC object
    // pass it the parameters object
    MCMC_RPM My_MCMC_RPM = MCMC_RPM(Params);
      
    //Sample next parameters from Metropolis Chain
	My_MCMC_Coast.RunMetropolisChain();
}
 