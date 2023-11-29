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
#include "./Parameters.hpp"
#include "./RockyCoastCRN.hpp"

using namespace std;

#ifndef MCMC_RPM_HPP
#define MCMC_RPM_HPP

class MCMC_RPM
{
    private:
      //Initialise Function
      void Initialise();
      void Initialise(Parameters InitialParams);
      
      //Vectors to hold extracted profile and CRN data
      int NProfileData, NCRNData;
      vector<double> ProfileXData, ProfileZData, ProfileZStdData;
      vector<double> CRNXData, CRNNData, CRNNErrorData;
      

      //Vector to hold StartTime input from parameter file
      double StartTime;

      //Vectors to hold Model results for comparison to data
      vector <vector<double>> CRNModelArray;
      vector<double> XModel, ZModel, CRNModel;
      vector<double> ZModelData;
      vector<double> CRNModelData;
      vector<double> BlankTopoDataVec, BlankCRNDataVec;

      //other useful stuff
      int NXModel;
      double CliffPositionX, XPos, InterpScale;
      
      //Declare parameters object
      Parameters Params;

      //Declare RPM object
      RPM MCMC_RPM;

      //declare RockyCoastCRN object
      RockyCoastCRN MCMC_RockyCoastCRN;

      //Declare a SeaLevel Object
      SeaLevel MCMC_Sealevel;

      //calculates the likelihood using measured and modelled data
      long double TopoLikelihood, CRNLikelihood;
      long double CalculateTopoLikelihood();
      long double CalculateCRNLikelihood();

      //runs a single iteration of the RPM model, then reports the likelihood of the parameters - this read in from driver?
      void RunCoastIteration();
      void ResetModel();

    public:

      /* MCMC_RPM Initialisation functions */

      /// @brief Empty initialisation function, will throw an error
      ///

      MCMC_RPM()
      {
          MCMC_RPM::Initialise();
      }

      /// @brief Initialisation function for MCMC_RPM objects
      /// @param Params A parameters object read from an parameter file
      
      MCMC_RPM(Parameters InitialParams)
      {
            
             MCMC_RPM::Initialise(InitialParams);
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
	  /// @param OutFilename File to write the results of each iteration of the chain to.
	  
      void RunMetropolisChain(int NIterations, char* ParamFilename, char* OutFilename);

  
      
      

};

#endif

