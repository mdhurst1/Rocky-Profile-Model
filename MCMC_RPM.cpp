/*

Comments

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
//#include <omp.h>
#include "RPM.hpp"
#include "MCMC_RPM.hpp"

using namespace std;

#ifndef MCMC_RPM_CPP
#define MCMC_RPM_CPP

void MCMC_RPM::Initialise()
{
    /* initialise an empty Markov Chain Monte Carlo object
        */
   cout << "unable to initialise, TestObject is an empty object" << endl;
   exit(EXIT_SUCCESS);
}

void MCMC_RPM::Initialise(char* ProfileDatafile)
{
    /*initialise a Markov Chain object with extracted platform profile
    */

   //Declare temp variables 
   char Dummy[32];
   float TempProfileXData, TempProfileZData;

   //Generate input filestream and read data into vectors
   ifstream READProfileDatafile(ProfileDatafile);
   if (!READProfileDatafile);
   { 
       printf("MCMC_Coast::%s line %d: Input Profile data file \"%s\" doesn't exist\n\n", __func__, __LINE__, ProfileDatafile);
       exit(EXIT_SUCCESS);
   }

    // ignore header lines by reading to Dummy
    // file format is...
    // X_header | Z_header
    //   X[0]   |   Z[0]
    //   X[1]   |   Z[1]
    //  X[...]  |  Z[...]
    //   X[n]   |   Z[n]
   ReadProfileDatafile >> Dummy >>  Dummy;
   while(ReadProfileDatafile >> TempProfileXData >> TempProfileZData)
   {
       ProfileXData.push_back(TempProfileXData);
       ProfileZData.push_back(TempProfileZData);
   }
   ProfileData = ProfileXData.size();

   RPM MCMCRPM = RPM();

   ////?????
}

long double MCMC_RPM::CalculateLikelihood()
{
    /* Function to calculate the likelihood by comparing measured and modelled data (dsm extracted and modelled)
    */

   //declarations
   double DiffX, Scale;
   long double Likelihood = 1.L;

     //Work out the modelled morphology
   XModel = MCMCPlatform.get_X();  
   ZModel = MCMCPlatform.get_Z();
   vector<double> TopoData(ProfileData);
   vector<double> Residuals(ProfileData);

   //Interpolate to extracted morphology X positions
   for (int i=0; i<ProfileData; ++i)
   {
       //Take X value of extracted morph position and interpolate to get model results at this point
       int j=0;
       while ((XModel[j]-ProfileXData[i]) <0) ++j;
       DiffX = XModel[j]-XData[i];
         Scale = DiffX/(XModel[j]-XModel[j-1]);
        
        //Get Interpolated Z value
        TopoData[i] = ZModel[j]-Scale*(ZModel[j]-Zmodel[j-1]);
   }
   
   //Calculate likelihood
   for (int i=0; i<ProfileData; ++i)
   {
       Residuals[i] = (ProfileZData[i]-TopoData[i])*(ProfileZData[i]-TopoData[i]);
       Likelihood *= exp(-(fabs(Residuals[i]))/(//find std of topodata to use for error)
   }
   return Likelihood;
}



long double MCMC_RPM::RunCoastIteration(double Resistance, double WeatheringRate)  ////?????
{
    /* runs a single instance of the RPM Model, then reportd the likelihood of the parameters
    */

   //Run a coastal iteration
   int WriteResultsFlag = 0;
    string OutfileName = "emptyfilename";
	MCMCPlatform.UpdateParameters(Resistance, WeatheringRate);
	bool Flag = MCMCPlatform.RunModel(OutfileName,WriteResultsFlag);

    //Calculate likelihood
    if (Flag == false) return CalculateLikelihood();
    else return -9999;
}

void MCMC_RPM::RunMetropolisChain(int NIterations, char* ParameterFilename, char* OutFilename)   ///??
{
    /* Run the metropolis algorithm along a chain with NIterations

      ParameterFilename is the name of the parameter file containing 
      Retreat Type, Minimum, Maximum, Standard Deviation (for Metropolis search) and 
      Initial Values for each parameter, the Platform Gradient, 
      Cliff Height and Tidal Amplitude (see example)
     
      Prints to the chain results to file 'OutFilename'
     
      Martin Hurst, January 2015 */

    //Declarations
	long double LastLikelihood = 0.L;			//Last accepted likelihood
	long double NewLikelihood = 0.L;			//New likelihood
	long double LikelihoodRatio = 0.L;			//Ratio between last and new likelihoods
	double AcceptanceProbability; //New iteration is accepted if likelihood ratio exceeds

    //int RetreatType;        //Style of cliff retreat 0 = single rate, 1 = step change in rates, 2 = linear change in rates
	//int BeachType = 0;      // Beach type is fixed width

	int NAccepted = 0;      //count accepted parameters
	int NRejected = 0;      //count rejected parameters

	double Rand1, Rand2;    //For generating random numbers

    //Holders to define parameter space	
	double  Resistance_New, Resistance_Old, Resistance_Min, Resistance_Max, Resistance_Std, Resistance_Init,
            WeatheringRate_New, WeatheringRate_Old, WeatheringRate_Min, WeatheringRate_Max, WeatheringRate_Std, WeatheringRate_Init,
           
	        //Parameters included in driver??? BermHeight, BeachSteepness, JunctionElevation, PlatformGradient, CliffHeight, CliffGradient, TidalAmplitude, SLR;

    double dFR, dK; //change in parameter values for Resistance (FR) and WeatheringRate (K)
    //double MeanChange = 0.; //Change in parameter values centred on zero allow changes in both directions(pos and neg) ???
  
    char Dummy[32];
    string RSLFilename, ScalingFilename;
        SLR = -9999;

    //Initialise seed for random number generation
    int RandomSeed = 1;
    srand(RandomSeed);

    //Create datafile out and write ParameterFilename
	ofstream ChainFileOut(OutFilename);
	ChainFileOut << "ParameterFile: " << ParameterFilename << endl;
	ChainFileOut  << "i Resistance_New WeatheingRate_New NewLikelihood LastLikelihood NAccepted NRejected" << endl;

    //Read in parameters for monte carlo run from parameter file
	//Min and max values for paramters from param file
	ifstream ParamFileIn(ParameterFilename);
	if (!ParameterFilename)
	{
	  printf("MCMC_Coast::%s: line %d Input parameter data file \"%s\" doesn't exist\n\n", __func__, __LINE__, ParameterFilename);
	  exit(EXIT_SUCCESS);
	}

    ParamFileIn //>> Dummy >> RetreatType
	            >> Dummy >> Resistance_Min >> Dummy >> Resistance_Max >> Dummy >> Resistance_Std >> Dummy >> Resistance_Init
	            >> Dummy >> WeatheringRate_Min >> Dummy >> WeatheringRate_Max >> Dummy >> WeatheringRate_Std >> Dummy >> WeatheringRate_Init
	            //>> Dummy >> BermHeight >> Dummy >> BeachSteepness >> Dummy >> JunctionElevation >> Dummy >> PlatformGradient 
                //>> Dummy >> CliffHeight >> Dummy >> CliffGradient >> Dummy >> TidalAmplitude
                //>> Dummy >> RSLFilename >> Dummy >> ScalingFilename >> Dummy >> WhichNuclideTemp;

    ParamFileIn.close();

    //Initialise RPM object
	MCMCPlatform = RPM(//RPM driver?? RetreatRate1_Init, RetreatRate2_Init, RetreatType, ChangeTime_Init, 
                                    //BeachWidth_Init, BeachType, BermHeight, BeachSteepness, PlatformGradient, 
                                    //CliffHeight, CliffGradient, JunctionElevation, TidalAmplitude,
                                    //SLR, SteppedPlatform, StepSize, WhichNuclides);

    //Initialise Sea Level history
    MCMCPlatform.InitialiseRSLData(RSLFilename);

    //Initialise Geomag Scaling
    MCMCPlatform.InitialiseScalingData(ScalingFilename);

    /*  start the chain with a guess this guess is a very coarse approximation of what the 'real' values 
	    might be. The Metropolis algorithm will sample around this */
	Resistance_New = Resistance_Init;
    WeatheringRate_New = WeatheringRate_Init;

    //Run a single coastal iteration to get the initial Likelihood for the initial parameters
	LastLikelihood = RunCoastIteration(Resistance_New, WeatheringRate_New);
	LastLikelihood = CalculateLikelihood();

    //set old parameters for comparison and updating
	Resistance_Old = Resistance_New;
	WeatheringRate_Old = WeatheringRate_New;

    int Accept = 0;

    //Do the metropolis algorithm
	for (int j=0; j<NIterations; ++j)
	{
		fflush(stdout);
 	  	printf("Iteration %d\r",j+1);  
		
		//Update the variables following a normal distribution

	  //First Update Resistance
		Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
			dFR = MeanChange + RetreatRate1_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
			RetreatRate1_New = RetreatRate1_Old + dFR;
			if ((Resistance_New < Resistance_Min) || (Resistance_New > Resistance_Max)) continue;
	  		else Accept = 1;
		}

        // Update WeatheringRate
	 	Accept = 0;
		while (Accept == 0)
        {	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
	  		dK = MeanChange + WeatheringRate_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
	  		WeatheringRate_New = WeatheringRate_Old + dK;
			if ((WeatheringRate_New < WeatheringRate_Min) || (WeatheringRate_New > WeatheringRate_Max)) continue; ////????
			else Accept = 1;
	  	}

      //Run a model iteration with new parameters
		NewLikelihood = RunCoastIteration(Resistance_New, WeatheringRate_New); 

      //Get the likelihood ratio
		LikelihoodRatio = NewLikelihood/LastLikelihood;

      //Get acceptance probability (from uniform distribution between 0 and 1)
		//This allows some false results to be accepted in order to explore the parameter space.
		AcceptanceProbability = (double)rand()/RAND_MAX;

      //Test for acceptance
		if (LikelihoodRatio > AcceptanceProbability)
		{
			LastLikelihood = NewLikelihood;
			++NAccepted;
			
			//Last variables become equal to new variables
			Resistance_Old = Resistance_New;
			WeatheringRate_Old = WeatheringRate_New;
		}
		else ++NRejected;

        //write the result to the output file
		ChainFileOut  << j << " " 
		              << Resistance_New << " " << WeatheringRate_New << " " 
		              << NewLikelihood << " " << LastLikelihood << " " 
		              << NAccepted << " " << NRejected << endl;
        }
    
    
	ChainFileOut.close();
    }
	
#endif



