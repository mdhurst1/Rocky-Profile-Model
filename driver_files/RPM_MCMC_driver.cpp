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
#include "../../src/RPM.hpp"
#include "../../src/RockyCoastCRN.hpp"
#include "../../src/SeaLevel.hpp"
#include "../../src/FastExp.hpp"

using namespace std;

template <typename T> string tostr(const T& t)
{ 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}


int main(int nNumberofArgs,char *argv[])
{

	//clock 
	time_t begin, end;
	time(&begin);

    // retrieve prameter file and workspace from arguments
    string Folder = argv[1];
	string TempParamFilename = argv[2];
	string InputParamFilename = Folder+TempParamFilename;
		
	// load parameter parser object
  	Parameters Params(Folder,InputParamFilename);
    
    // initiate MCMC object
    // pass it the parameters object
    MCMC_RPM My_MCMC_RPM = MCMC_RPM(Params);
      
    // declare iteration counter
    while (i=0; i<My_MCMC_Coast.NIterations; ++i)
    {
        //Sample next parameters from Metropolis Chain
	    My_MCMC_Coast.SampleMetropolisChain();

        // Initiate RPM object
        RPM PlatformModel = RPM(dZ, dX, Params.InitialGradient, Params.CliffElevation, Params.MaxElevation, Params.MinElevation);

        //initialise RockyCoastCRN friend object
        RockyCoastCRN PlatformCRN = RockyCoastCRN();

        if (CRNFlag)
        {
            //Which Nuclides to track 10Be, 14C, 26Al, 36Cl?
            vector<int> Nuclides;
            Nuclides.push_back(10);
            
            //initialise RockyCoastCRN friend object
            PlatformCRN = RockyCoastCRN(PlatformModel, Nuclides);
        }
        
        //Initialise Tides
        double TidalPeriod = 12.42;
        PlatformModel.InitialiseTides(TidalRange);
        if (CRNFlag) PlatformCRN.InitialiseTides(TidalRange/2.,TidalPeriod);
            
        //Initialise Waves
        //Single Wave for now but could use the waveclimate object from COVE!?
        double WaveHeight_Mean = 3.;
        double WaveHeight_StD = 0.;
        double WavePeriod_Mean = 6.;
        double WavePeriod_StD = 0;
        PlatformModel.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

        // Wave coefficient constant
        double StandingCoefficient = 0.1;
        double BreakingCoefficient = 10.;
        double BrokenCoefficient = 1.;
        PlatformModel.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient, WaveAttenuationConst);

        //reset the geology
        double CliffFailureDepth = 0.1;
        PlatformModel.InitialiseGeology(CliffElevation, CliffFailureDepth, Resistance, WeatheringRate, SubtidalEfficacy);
        
        // Pass Metropolis Chain Parameters

        // Run the model
        //Loop through time
        while (Time >= EndTime)
        {
            //Do an earthquake?
            if (Time < UpliftTime)
            {
                // perform uplift
                PlatformModel.TectonicUplift(UpliftMagnitude);
                UpliftTime -= UpliftFrequency;
                
                // Update the Morphology 
                PlatformModel.UpdateMorphology();
            }		

            //Update Sea Level
            InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);
            PlatformModel.UpdateSeaLevel(InstantSeaLevel);

            //Get the wave conditions
            PlatformModel.GetWave();

            //Calculate forces acting on the platform
            PlatformModel.CalculateBackwearing();
            PlatformModel.CalculateDownwearing();

            //Do erosion
            PlatformModel.ErodeBackwearing();
            PlatformModel.ErodeDownwearing();

            //Update the Morphology 
            PlatformModel.UpdateMorphology();	
            
            //Implement Weathering
            PlatformModel.IntertidalWeathering();
            PlatformModel.SubtidalWeathering();
            
            //Update the Morphology 
            PlatformModel.UpdateMorphology();

            //Check for Mass Failure
            PlatformModel.MassFailure();
            
            //Update the Morphology 
            PlatformModel.UpdateMorphology();

            //Update the morphology inside RockyCoastCRN
            if (CRNFlag) PlatformCRN.UpdateMorphology(PlatformModel);

            //Update the CRN concentrations
            if (CRNFlag) PlatformCRN.UpdateCRNs();
                
            //print?
            if (Time <= PrintTime)
            {
                cout.flush();
                cout << "RPM: Time " << setprecision(2) << fixed << Time << " years\r";
                //PlatformModel.WriteProfile(OutputFileName, Time);  //This is for testing - need to remove
                        //if (CRNFlag) PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
                PrintTime -= PrintInterval;
            }

            //update time
            Time -= TimeInterval;
        }
        // calculate resulting likelihoods
        // calculate CRN likelihood
        // calculate topo likelihood
        // accept or reject new parameters

        // save results

    }
 
}