/*

comments

*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "./RPM.hpp"
#include "./SeaLevel.hpp"

using namespace std;

#ifndef MCMC_RPM_HPP
#define MCMC_RPM_HPP

class MCMC_RPM
{
    private:
      //Initialise Function
      void Initialise();
      void Initialise(char* ProfileDatafile);
      //void Initialise(char* ProfileDatafile, char* PlatformXSectionDatafile);

      //Vectors to hold extracted profile data
      int NProfileData;
      vector<double> ProfileXData;
      vector<double> ProfileZData;
      double ZStd;

      //Vector to hold StartTime input from parameter file
      double StartTime;

      //Vectors to hold Topographic Model results
      int NTopoData;
      vector<double> TopoXData;
      vector<double> TopoZData;

      //Declare RPM object ??
      RPM MCMCPlatform;

      //Declare a SeaLevel Object
      SeaLevel MCMCSeaLevel;

      //calculates the likelihood using measured and modelled data
      long double CalculateLikelihood();

      //runs a single iteration of the RPM model, then reports the likelihood of the parameters - this read in from driver?
      long double RunCoastIteration();

    public:

      /* MCMC_RPM Initialisation functions */

      /// @brief Empty initialisation function, will throw an error
      ///

      MCMC_RPM()
      {
          MCMC_RPM::Initialise();
      }

      /// @brief Initialisation function for MCMC_RPM objects
      /// @param ProfileDatafile Data file containing extracted platform morphology
      /// @param CRNDatafile containing observed CRN concentrations (not currently used).

      MCMC_RPM(char* ProfileDatafile)
      {
            
             MCMC_RPM::Initialise(ProfileDatafile);
      }

      /// this runs the metropolis algorithm along a chain with NIterations it prints to the file 'OutFilename'
	  /// @brief Launch the main MCMC program loop to search for most likely parameters
	  /// @details This function runs a Markov Chain Monte Carlo analysis to find the most likely
	  ///   combination of parameters to fit observed concentrations of 10Be from a coastal platform.
	  ///   The MCMC algorithm runs the RockyCoastCRN model many times (NIterations c. 200k), checks the 
	  ///   likelihood between modelled and measured concentrations and either accepts or rejects the 
	  ///   new parameters. The model sometimes accepts less likely results to allow exploration of the
	  ///   parameter space. 
	  /// @param NIterations Number of times to run the RockyCoastCRN model in the Markov Chain
	  /// @param ParamFilename Filename of file containing parameters for the MCMC analysis (see example file)
	  /// @param OutFilename File to write the results of each iteration of the chain to.
	  
      void RunMetropolisChain(int NIterations, char* ParamFilename, char* OutFilename);

  
      
//void ResetModel() - where to place reset model & update geology?
		//{
			//Initialise(dZ,dX, InitialGradient, CliffHeight, MinimumElevation);
		//}

};

#endif