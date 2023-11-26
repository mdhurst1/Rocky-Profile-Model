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
#include "RPM.hpp"
#include "MCMC_RPM.hpp"
#include "SeaLevel.hpp"

using namespace std;

#ifndef MCMC_RPM_CPP
#define MCMC_RPM_CPP

void MCMC_RPM::Initialise()
{
    /* initialise an empty Markov Chain Monte Carlo object
        */
   cout << "unable to initialise, MCMC_RPM is an empty object" << endl;
   cout << "you need to specify a Parameters object when initialising the MCMC object" << endl;
   exit(EXIT_SUCCESS);
}

void MCMC_RPM::Initialise(Parameters Params)
{
    /*initialise a Markov Chain object with extracted platform profile and concentrations
    */


    // file name params
    string ProfileDatafile(Params.TopoFilename);
    string ConcentrationDatafile(Params.CRNFilename);

    // some dummy variables for reading from file
    float TempXData, TempZData, TempNData;
    
    //Generate input filestream and read data into vectors
    ifstream READProfileDatafile(ProfileDatafile);
    if (!READProfileDatafile)
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
    READProfileDatafile >> Dummy >>  Dummy;
    while(READProfileDatafile >> TempProfileXData >> TempProfileZData)
    {
        ProfileXData.push_back(TempProfileXData);
        ProfileZData.push_back(TempProfileZData);
    }

    READProfileDatafile.close();

    // get size of the profile data vectors
    NProfileData = ProfileXData.size();

    //Generate input filestream and read data into vectors
    ifstream READCRNDatafile(CRNDatafile);
    if (!READCRNDatafile)
    { 
        printf("MCMC_Coast::%s line %d: Input CRN data file \"%s\" doesn't exist\n\n", __func__, __LINE__, ProfileDatafile);
        exit(EXIT_SUCCESS);
    }

    // ignore header lines by reading to Dummy
    // file format is...
    // X_header | N_header
    //   X[0]   |   N[0]
    //   X[1]   |   N[1]
    //  X[...]  |  N[...]
    //   X[n]   |   N[n]
    READCRNDatafile >> Dummy >>  Dummy;
    while(READCRNDatafile >> TempXData >> TempNData)
    {
        CRNXData.push_back(TempXData);
        CRNNData.push_back(TempNData);
    }

    READCRNDatafile.close();

    // get size of the profile data vectors
    NCRNData = CRNXData.size();

    // initialise RPM object with default morphology
    RPM MCMC_RPM = RPM(Params.dZ, Params.dX, Params.InitialGradient, Params.CliffElevation, Params.MaxElevation, Params.MinElevation);

    //initialise RockyCoastCRN friend object
	RockyCoastCRN MCMC_RockyCoastCRN;

	// THIS SHOULD BE IN PARAMETER FILE
	if (Params.CRN_Predictions)
	{
		//Which Nuclides to track 10Be, 14C, 26Al, 36Cl?
		vector<int> Nuclides;
        if (Params.Berylium) Nuclides.push_back(10);
        if (Params.Carbon) Nuclides.push_back(14);
        if (Params.Aluminium) Nuclides.push_back(26);
		
		//initialise RockyCoastCRN friend object
		MCMC_RockyCoastCRN = RockyCoastCRN(PlatformModel, Nuclides);
	}

    // initialise sea level object
    if (Params.ReadSeaLevelFromFile) MCMC_Sealevel = SeaLevel(Params.SeaLevelFilename);
	else MCMC_Sealevel = SeaLevel(Params.SeaLevelRise, Params.StartTime, Params.EndTime, Params.TimeStep);

    // Get initial sea level
	float InstantSeaLevel = MCMC_SeaLevel.get_SeaLevel(Params.StartTime);
	MCMC_RPM.UpdateSeaLevel(InstantSeaLevel);

	//Initialise Tides
	MCMC_RPM.InitialiseTides(Params.TidalRange);
    if (Params.CRN_Predictions) MCMC_RockyCoastCRN.InitialiseTides(Params.TidalRange/2.,Params.TidalPeriod);
	
	//Initialise Waves
	MCMC_RPM.InitialiseWaves(Params.WaveHeight_Mean, Params.WaveHeight_StD, Params.WavePeriod_Mean, Params.WavePeriod_StD);
}

void MCMC_RPM::RunMetropolisChain(int NIterations, char* ParameterFilename, char* OutFilename)
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

	int NAccepted = 0;      //count accepted parameters
	int NRejected = 0;      //count rejected parameters

	double Rand1, Rand2;    //For generating random numbers

    //Holders to define parameter space	
	double  Resistance_New, Resistance_Old, Resistance_Min, Resistance_Max, Resistance_Std, Resistance_Init,
            WeatheringRate_New, WeatheringRate_Old, WeatheringRate_Min, WeatheringRate_Max, WeatheringRate_Std, WeatheringRate_Init,
            SubmarineDecayConst;
           
	        //Parameters included in driver BermHeight, BeachSteepness, JunctionElevation, PlatformGradient, CliffHeight, CliffGradient, TidalAmplitude, SLR;
    double TidalRange;
    double dFR, dK; //change in parameter values for Resistance (FR) and WeatheringRate (K)
    double MeanChange = 0.; //Change in parameter values centred on zero allow changes in both directions(pos and neg)
  
    // morphology parameters
    double dZ, dX, Gradient, CliffHeight, CliffFailureDepth, MinElevation;
    
    char Dummy[32];
    string RSLFilename, ScalingFilename;
    
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

    ParamFileIn >> Dummy >> Resistance_Min >> Dummy >> Resistance_Max >> Dummy >> Resistance_Std >> Dummy >> Resistance_Init
	            >> Dummy >> WeatheringRate_Min >> Dummy >> WeatheringRate_Max >> Dummy >> WeatheringRate_Std >> Dummy >> WeatheringRate_Init
	            >> Dummy >> TidalRange >> Dummy >> dZ >> Dummy >> dX >> Dummy >> Gradient >> Dummy >> CliffHeight >> Dummy >> MinElevation
                >> Dummy >> SubmarineDecayConst >> Dummy >> RSLFilename >> Dummy >> StartTime >> Dummy >> ZStd;

    ParamFileIn.close();

    //cout << "Resistance_Min: " << Resistance_Min << endl;


    //Initialise RPM object
	MCMCPlatform = RPM(dZ, dX, Gradient, CliffHeight, MinElevation);

    //Initialise Sea Level history
    MCMCSeaLevel = SeaLevel(RSLFilename);
    
    //
    double InitialSeaLevel = MCMCSeaLevel.get_SeaLevel(StartTime);

    //
    MCMCPlatform.UpdateSeaLevel(InitialSeaLevel);

    //Initialise Tides
	MCMCPlatform.InitialiseTides(TidalRange);

    //Initialise Waves
	//Single Wave for now but could use the waveclimate object from COVE!?
	double WaveHeight_Mean = 3.;
	double WaveHeight_StD = 0.;
	double WavePeriod_Mean = 6.;
	double WavePeriod_StD = 0;
	MCMCPlatform.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

    // Wave coefficient constant
	double StandingCoefficient = 0.1;
	double BreakingCoefficient = 10.;
	double BrokenCoefficient = 1.;
	double WaveAttenuationConst = 0.01;
	MCMCPlatform.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient, WaveAttenuationConst);

    //initialise the geology
	MCMCPlatform.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistance_Init, WeatheringRate_Init, MinElevation);

    /*  start the chain with a guess this guess is a very coarse approximation of what the 'real' values 
	    might be. The Metropolis algorithm will sample around this */
	Resistance_New = Resistance_Init;
    WeatheringRate_New = WeatheringRate_Init;

    //Run a single coastal iteration to get the initial Likelihood for the initial parameters
	LastLikelihood = RunCoastIteration();
	LastLikelihood = CalculateLikelihood();

    //set old parameters for comparison and updating
	Resistance_Old = Resistance_New;
	WeatheringRate_Old = WeatheringRate_New;

    int Accept = 0;

    //Do the metropolis algorithm
	for (int j=0; j<NIterations; ++j)
	{
		//fflush(stdout);
 	  	printf("Iteration %d\n",j+1);  
		
		//Update the variables following a normal distribution

	  //First Update Resistance
		Accept = 0;
		while (Accept == 0)
		{	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
			dFR = MeanChange + Resistance_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
			Resistance_New = Resistance_Old + dFR;
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
			if ((WeatheringRate_New < WeatheringRate_Min) || (WeatheringRate_New > WeatheringRate_Max)) continue;
			else Accept = 1;
	  	}
        
        //Reset parameters to be read from input file
        //reset model
        MCMCPlatform.ResetModel();
        //reset morphology (input parameters)
        //MCMCPlatform.ResetMorphology(dZ, dX, Gradient, CliffHeight, CliffFailureDepth, MinElevation);
        //reset sea level
        MCMCPlatform.UpdateSeaLevel(InitialSeaLevel);
        //reset tides
        MCMCPlatform.InitialiseTides(TidalRange);
        // initialise waves
        double WaveHeight_Mean = 3.;
	    double WaveHeight_StD = 0.;
	    double WavePeriod_Mean = 6.;
	    double WavePeriod_StD = 0;
	    MCMCPlatform.InitialiseWaves(WaveHeight_Mean, WaveHeight_StD, WavePeriod_Mean, WavePeriod_StD);

        // Wave coefficient constant
	    double StandingCoefficient = 0.1;
	    double BreakingCoefficient = 10.;
	    double BrokenCoefficient = 1.;
	    double WaveAttenuationConst = 0.01;
	    MCMCPlatform.Set_WaveCoefficients(StandingCoefficient, BreakingCoefficient, BrokenCoefficient, WaveAttenuationConst);

        // reset model parameters with new values
        MCMCPlatform.InitialiseGeology(CliffHeight, CliffFailureDepth, 
                                        Resistance_New,  WeatheringRate_New, SubmarineDecayConst);

        //Run a model iteration with new parameters
		NewLikelihood = RunCoastIteration();

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


long double MCMC_RPM::RunCoastIteration()  
{
    /* runs a single instance of the RPM Model, then reported the likelihood of the parameters
    */

    //Time control parameters
	//Time runs in yrs bp
	double EndTime = 0;
	double Time = StartTime;
	double TimeInterval = 1;
    double InstantSeaLevel;
    
    double PrintTime = StartTime;
    double PrintInterval = 100;

    //reset the model domain
	//MCMCPlatform.ResetModel();

    //Loop through time
	while (Time >= EndTime)
	{
		// print time to screen
        if (Time < PrintTime)
        {
            printf("Time %4.f\n",Time);
            PrintTime -= PrintInterval;
        }

        //set up if statement to only print every 100/1000 years? 
        

        //Update Sea Level
		InstantSeaLevel = MCMCSeaLevel.get_SeaLevel(Time);
		MCMCPlatform.UpdateSeaLevel(InstantSeaLevel);

		//Get the wave conditions
		MCMCPlatform.GetWave();

        //Implement Weathering
		MCMCPlatform.IntertidalWeathering();

		//Calculate forces acting on the platform
		MCMCPlatform.CalculateBackwearing();
		MCMCPlatform.CalculateDownwearing();

		//Do erosion
		MCMCPlatform.ErodeBackwearing();
		MCMCPlatform.ErodeDownwearing();
	
		//Check for Mass Failure
		MCMCPlatform.MassFailure();
		
		//Update the Morphology 
		MCMCPlatform.UpdateMorphology();
		
		//update time
		Time -= TimeInterval;
	}
       
    //Calculate likelihood
    return CalculateLikelihood();    
}

long double MCMC_RPM::CalculateLikelihood()
{
    /* Function to calculate the likelihood by comparing measured and modelled data (dsm extracted and modelled)
    */

   //declarations
   double Scale;
   long double Likelihood = 1.L;
   vector<double> XModel, ZModel;

     //Work out the modelled morphology
   XModel = MCMCPlatform.get_X(); 
   ZModel = MCMCPlatform.get_Elevations();
   int XSize = XModel.size();
   double CliffPositionX = XModel[XSize-1];
   
   vector<double> XPos(NProfileData);
   vector<double> TopoData(NProfileData);
   vector<double> Residuals(NProfileData);
   vector<double> DiffX(NProfileData);

   //Interpolate to extracted morphology X positions
   for (int i=0; i<NProfileData; ++i)
   {
       
       //Normalising profile data to modelled cliff position - using Swath profile data where cliff position = 0
       XPos[i] = CliffPositionX - ProfileXData[i];


       //Take X value of extracted morph position and interpolate to get model results at this point
       int j=0;
       while ((XModel[j]- XPos[i]) <0) ++j;
       DiffX[i] = XModel[j] - XPos[i];
         Scale = DiffX[i]/(XModel[j]-XModel[j-1]);
        
        //Get Interpolated Z value
        TopoData[i] = ZModel[j]-Scale*(ZModel[j]-ZModel[j-1]);
   }
   
   //Calculate likelihood
   for (int i=0; i<NProfileData; ++i)
   {
       Residuals[i] = (ProfileZData[i]-TopoData[i])*(ProfileZData[i]-TopoData[i]);
       Likelihood *= exp(-(fabs(Residuals[i]))/(ZStd*ZStd));    //ZStd read in from parameter file?
   }
   return Likelihood;
}


#endif



