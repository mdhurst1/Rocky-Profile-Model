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

void MCMC_RPM::Initialise(Parameters InitialParams)
{
    /*initialise a Markov Chain object with extracted platform profile and concentrations
    */
    Params = InitialParams;

    // some dummy variables for reading from file
    float TempProfileXData, TempProfileZData, TempProfileZStdData;
    float TempXData, TempNData, TempNErrorData;
    char Dummy[32];

    //Generate input filestream and read data into vectors
    ifstream READProfileDatafile(Params.TopoFilename);
    if (!READProfileDatafile)
    { 
        printf("MCMC_Coast::%s line %d: Input Profile data file \"%s\" doesn't exist\n\n", __func__, __LINE__, Params.TopoFilename.c_str());
        exit(EXIT_SUCCESS);
    }

    // ignore header lines by reading to Dummy
    // file format is...
    // X_header | Z_header | Zstd_header
    //   X[0]   |   Z[0]   |  Zstd[0]
    //   X[1]   |   Z[1]   |  Zstd[1]
    //  X[...]  |  Z[...]  | Zstd[...]
    //   X[n]   |   Z[n]   |  Zstd[n]

    READProfileDatafile >> Dummy >>  Dummy >> Dummy;
    while(READProfileDatafile >> TempProfileXData >> TempProfileZData >> TempProfileZStdData)
    {
        ProfileXData.push_back(TempProfileXData);
        ProfileZData.push_back(TempProfileZData);
        ProfileZStdData.push_back(TempProfileZStdData);
    }

    READProfileDatafile.close();

    // get size of the profile data vectors
    NProfileData = ProfileXData.size();
    vector<double> BlankTopo(NProfileData);
    BlankTopoDataVec = BlankTopo;

    //Generate input filestream and read data into vectors
    ifstream READCRNDatafile(Params.CRNFilename);
    if (!READCRNDatafile)
    { 
        printf("MCMC_Coast::%s line %d: Input CRN data file \"%s\" doesn't exist\n\n", __func__, __LINE__, Params.CRNFilename.c_str());
        exit(EXIT_SUCCESS);
    }

    // ignore header lines by reading to Dummy
    // file format is...
    // X_header | N_header |   N_Error
    //   X[0]   |   N[0]   |  N_Err[0]
    //   X[1]   |   N[1]   |  N_Err[1]
    //  X[...]  |  N[...]  | N_Err[...]
    //   X[n]   |   N[n]   |  N_Err[n]
    READCRNDatafile >> Dummy >>  Dummy >> Dummy;
    while(READCRNDatafile >> TempXData >> TempNData >> TempNErrorData)
    {
        CRNXData.push_back(TempXData);
        CRNNData.push_back(TempNData);
        CRNNErrorData.push_back(TempNErrorData);
    }

    READCRNDatafile.close();

    // get size of the profile data vectors
    NCRNData = CRNXData.size();
    vector<double> BlankCRN(NCRNData);
    BlankCRNDataVec = BlankCRN;

    // initialise RPM object with default morphology
    RPM MCMC_ProfileModel = RPM(Params.dZ, Params.dX, Params.InitialGradient, Params.CliffElevation, Params.MaxElevation, Params.MinElevation);

    //initialise RockyCoastCRN friend object
	RockyCoastCRN MCMC_RockyCoastCRN;

	// THIS SHOULD BE IN PARAMETER FILE
	if (Params.CRN_Predictions)
	{
		//initialise RockyCoastCRN friend object
		MCMC_RockyCoastCRN = RockyCoastCRN(MCMC_ProfileModel, Params.Nuclides);
	}

    // initialise sea level object
    if (Params.ReadSeaLevelFromFile) MCMC_Sealevel = SeaLevel(Params.SeaLevelFilename);
	else MCMC_Sealevel = SeaLevel(Params.SeaLevelRise, Params.StartTime, Params.EndTime, Params.TimeStep);

    // Get initial sea level
	float InstantSeaLevel = MCMC_Sealevel.get_SeaLevel(Params.StartTime);
	MCMC_ProfileModel.UpdateSeaLevel(InstantSeaLevel);

	//Initialise Tides
	MCMC_ProfileModel.InitialiseTides(Params.TidalRange);
    if (Params.CRN_Predictions) MCMC_RockyCoastCRN.InitialiseTides(Params.TidalRange/2.,Params.TidalPeriod);
	
	//Initialise Waves
	MCMC_ProfileModel.InitialiseWaves(Params.WaveHeight_Mean, Params.WaveHeight_StD, Params.WavePeriod_Mean, Params.WavePeriod_StD);

    // Wave coefficient constant
	MCMC_ProfileModel.Set_WaveCoefficients(Params.StandingWaveCoef, Params.BreakingWaveCoef, 
										Params.BrokenWaveCoef, Params.WaveAttenuationConst);

	//reset the geology
	MCMC_ProfileModel.InitialiseGeology(Params.CliffElevation, Params.CliffFailureDepth, Params.Resistance, 
									Params.WeatheringRate, Params.SubtidalEfficacy);	

}

void MCMC_RPM::RunMetropolisChain(int NIterations, char* ParameterFilename, char* OutFilename)
{
    /* Run the metropolis algorithm along a chain with NIterations

      Prints to the chain results to file 'OutFilename'
     
      Martin Hurst, January 2015 */

    //Declarations
	long double LastLikelihood = 0.L;			//Last accepted likelihood
	long double NewLikelihood = 0.L;			//New likelihood
	long double LikelihoodRatio = 0.L;			//Ratio between last and new likelihoods
	long double CombinedLikelihood = 0.L;
    double AcceptanceProbability;               //New iteration is accepted if likelihood ratio exceeds

	int NAccepted = 0;      //count accepted parameters
	int NRejected = 0;      //count rejected parameters

	double Rand1, Rand2;    //For generating random numbers

    //Holders to define parameter space	
	double  Resistance_New, Resistance_Old, Resistance_Min, Resistance_Max, Resistance_Std,
            WeatheringRate_New, WeatheringRate_Old, WeatheringRate_Min, WeatheringRate_Max, WeatheringRate_Std,
            WaveAttenuation_New, WaveAttenuation_Old, WaveAttenuation_Min, WaveAttenuation_Max, WaveAttenuation_Std;
           
	//Initialise seed for random number generation
    int RandomSeed = 1;
    srand(RandomSeed);

    //Create datafile out and write ParameterFilename
	ofstream ChainFileOut(OutFilename);
	ChainFileOut << "ParameterFile: " << ParameterFilename << endl;
	ChainFileOut  << "i Resistance_New WeatheingRate_New WaveAttenuation_New TopoLikelihood CRNLikelihood NewLikelihood LastLikelihood NAccepted NRejected" << endl;

    /*  start the chain with a random guess this guess is a very coarse approximation of what the 'real' values 
	    might be. The Metropolis algorithm will sample around this */
	Resistance_New = Params.Resistance_Min + ((double)rand()/RAND_MAX)*(Params.Resistance_Max - Params.Resistance_Min);
    WeatheringRate_New = Params.WeatheringRate_Min + ((double)rand()/RAND_MAX)*(Params.WeatheringRate_Max - Params.WeatheringRate_Min);
    WaveAttenuation_New = Params.WaveAttenuation_Min + ((double)rand()/RAND_MAX)*(Params.WaveAttenuation_Max - Params.WaveAttenuation_Min);

    //Run a single coastal iteration to get the initial Likelihoods for the initial parameters
	RunCoastIteration();
	CalculateTopoLikelihood();
    CalculateCRNLikelihood();
    CombinedLikelihood = Params.TopoWeighting*TopoLikelihood + Params.CRNWeighting*CRNLikelihood;
    LastLikelihood = CombinedLikelihood;

    //set old parameters for comparison and updating
	Resistance_Old = Resistance_New;
	WeatheringRate_Old = WeatheringRate_New;
    WaveAttenuation_Old = WaveAttenuation_New;

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
			Resistance_New = Resistance_Old + Resistance_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));
			if ((Resistance_New < Resistance_Min) || (Resistance_New > Resistance_Max)) continue;
	  		else Accept = 1;
		}

        // Update WeatheringRate
	 	Accept = 0;
		while (Accept == 0)
        {	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
	  		WeatheringRate_New = WeatheringRate_Old + WeatheringRate_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));;
			if ((WeatheringRate_New < WeatheringRate_Min) || (WeatheringRate_New > WeatheringRate_Max)) continue;
			else Accept = 1;
	  	}
        
        // Update WaveAttenuation
	 	Accept = 0;
		while (Accept == 0)
        {	
			Rand1 = (double)rand()/RAND_MAX; Rand2 = (double)rand()/RAND_MAX;
	  		WaveAttenuation_New = WaveAttenuation_Old + WaveAttenuation_Std*sqrt(-2.*log(Rand1))*cos(2.*M_PI*(Rand2));;
			if ((WaveAttenuation_New < WaveAttenuation_Min) || (WaveAttenuation_New > WaveAttenuation_Max)) continue;
			else Accept = 1;
	  	}

        // add here tectonics if free parameters?!

        //Reset parameters to be read from input file
        //reset model
        MCMC_ProfileModel.ResetModel();

        //Run a single coastal iteration to get the initial Likelihoods for the initial parameters
        RunCoastIteration();
        CalculateTopoLikelihood();
        CalculateCRNLikelihood();
        CombinedLikelihood = Params.TopoWeighting*TopoLikelihood + Params.CRNWeighting*CRNLikelihood;
        
        //Get the likelihood ratio
		LikelihoodRatio = CombinedLikelihood/LastLikelihood;

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
            WaveAttenuation_Old = WaveAttenuation_New;
		}
		else ++NRejected;

        //write the result to the output file
		ChainFileOut  << j << " " 
		              << Resistance_New << " " << WeatheringRate_New << " " << WaveAttenuation_New << " " 
		              << TopoLikelihood << " " << CRNLikelihood << " " << CombinedLikelihood << " " << LastLikelihood << " " 
		              << NAccepted << " " << NRejected << endl;
        }
    
    
	ChainFileOut.close();
    }


void MCMC_RPM::RunCoastIteration()  
{
    /* runs a single instance of the RPM Model, then reported the likelihood of the parameters
    */

    //Time control parameters
	//Time runs in yrs bp
	double EndTime = Params.EndTime;
	double Time = Params.StartTime;
	double TimeInterval = Params.TimeStep;
    double InstantSeaLevel;
    double UpliftTime = Time - Params.UpliftFrequency;
    double PrintTime = Time;

    //reset the model domain
	//MCMC_ProfileModel.ResetModel();

    //Loop through time
	while (Time >= EndTime)
	{
		//set up if statement to only print every 100/1000 years? 
        //Do an earthquake?
		if (Params.Earthquakes && Time < UpliftTime)
		{
			MCMC_ProfileModel.TectonicUplift(Params.UpliftMagnitude);
			UpliftTime -= Params.UpliftFrequency;
			
			//Update the Morphology 
			MCMC_ProfileModel.UpdateMorphology();
		}		

        //Update Sea Level
		InstantSeaLevel = MCMC_Sealevel.get_SeaLevel(Time);
		MCMC_ProfileModel.UpdateSeaLevel(InstantSeaLevel);

		//Get the wave conditions
		MCMC_ProfileModel.GetWave();

        //Implement Weathering
		MCMC_ProfileModel.IntertidalWeathering();

		//Calculate forces acting on the platform
		MCMC_ProfileModel.CalculateBackwearing();
		MCMC_ProfileModel.CalculateDownwearing();

		//Do erosion
		MCMC_ProfileModel.ErodeBackwearing();
		MCMC_ProfileModel.ErodeDownwearing();
	
		//Check for Mass Failure
		MCMC_ProfileModel.MassFailure();
		
		//Update the Morphology 
		MCMC_ProfileModel.UpdateMorphology();

        //Update the morphology and CRNs inside RockyCoastCRN
		if (Params.CRN_Predictions) 
		{
			MCMC_RockyCoastCRN.UpdateMorphology(MCMC_ProfileModel);
			MCMC_RockyCoastCRN.UpdateCRNs();
		}
		
        //print?
		if (Time <= PrintTime)
		{
			MCMC_ProfileModel.WriteProfile(Params.ProfileOutFilename, Time);
			if (Params.CRN_Predictions) MCMC_RockyCoastCRN.WriteCRNProfile(Params.ConcentrationsOutFilename, Time);
			PrintTime -= Params.PrintInterval;
		}

		//update time
		Time -= TimeInterval;
	}
}

void MCMC_RPM::ResetModel()
{
    // will need to think about how to use new parameters during reset.
    
    // reinitialise RPM object with default morphology
    MCMC_ProfileModel = RPM(Params.dZ, Params.dX, Params.InitialGradient, Params.CliffElevation, Params.MaxElevation, Params.MinElevation);

    //reinitialise RockyCoastCRN friend object
	MCMC_RockyCoastCRN = RockyCoastCRN(MCMC_ProfileModel, Params.Nuclides);

    // Get initial sea level
	float InstantSeaLevel = MCMC_Sealevel.get_SeaLevel(Params.StartTime);
	MCMC_ProfileModel.UpdateSeaLevel(InstantSeaLevel);

	//Initialise Tides
	MCMC_ProfileModel.InitialiseTides(Params.TidalRange);
    if (Params.CRN_Predictions) MCMC_RockyCoastCRN.InitialiseTides(Params.TidalRange/2.,Params.TidalPeriod);
	
	//Initialise Waves
	MCMC_ProfileModel.InitialiseWaves(Params.WaveHeight_Mean, Params.WaveHeight_StD, Params.WavePeriod_Mean, Params.WavePeriod_StD);

    // Wave coefficient constant
	MCMC_ProfileModel.Set_WaveCoefficients(Params.StandingWaveCoef, Params.BreakingWaveCoef, 
										Params.BrokenWaveCoef, Params.WaveAttenuationConst);

	//reset the geology
	MCMC_ProfileModel.InitialiseGeology(Params.CliffElevation, Params.CliffFailureDepth, Params.Resistance, 
									Params.WeatheringRate, Params.SubtidalEfficacy);
}

long double MCMC_RPM::CalculateTopoLikelihood()
{
    /* Function to calculate the likelihood by comparing measured and modelled data (dsm extracted and modelled)
    */

    // reset likelihood
    TopoLikelihood = 1.L;
    
    //declarations
    XModel = MCMC_ProfileModel.get_X(); 
    ZModel = MCMC_ProfileModel.get_Elevations();
    ZModelData = BlankTopoDataVec;
    NXModel = XModel.size();
    CliffPositionX = XModel[NXModel-1];
           
    //Interpolate to extracted morphology X positions
    for (int i=0; i<NProfileData; ++i)
    {
        //Normalising profile data to modelled cliff position (using Swath profile data where cliff position of measured data = 0)
        XPos = CliffPositionX - ProfileXData[i];
        
        //Take X value of extracted morph position and interpolate to get model results at this point
        int j=0;
        while ((XModel[j]- XPos) <0) ++j; //XPos starts at point nearest cliff and works offshore - starts at 40m from cliff so will never =0
        InterpScale = (XModel[j] - XPos)/(XModel[j]-XModel[j-1]);

        //Get Interpolated Z value
        ZModelData[i] = ZModel[j]-InterpScale*(ZModel[j]-ZModel[j-1]);
        
    }

    //calculate likelihood from topo residuals
    for (int i=0; i<NProfileData; ++i)
    {
        //this was Jen's calcs for Dakota, which read RMSE
        //Residuals[i] = fabs(ProfileZData[i]-TopoData[i]);
        //TotalResiduals += pow(ProfileZData[i]-TopoData[i],2);
        TopoLikelihood *= exp(-(fabs((ProfileZData[i]-ZModelData[i])*(ProfileZData[i]-ZModelData[i])))/(ProfileZStdData[i]*ProfileZStdData[i]));    //ZStd read in from parameter file?
    }
    return TopoLikelihood;
}

long double MCMC_RPM::CalculateCRNLikelihood()
{
    ///////////////////////////////////////                                
    //                                   //
    //Calculations for CRN concentrations//
    //                                   //
    ///////////////////////////////////////

    // reset likelihood
    CRNLikelihood = 1.L;
    
    //declarations
    XModel = MCMC_RockyCoastCRN.get_X();  // does this exist!?
    NXModel = XModel.size();
    
    // get CRN model array and just sample first element (i.e. only 1 nuclide)
    // this will need to be revamped if we are ever to work with multiple nuclides
    CRNModelArray = MCMC_RockyCoastCRN.get_SurfaceN(); // does this exist?!
    vector<double> BlankCRNModel(NXModel);
    CRNModel = BlankCRNModel;
    for (int j=0; j<NXModel; ++j) CRNModel.push_back(CRNModelArray[j][0]);
    
    CRNModelData = BlankCRNDataVec;
    CliffPositionX = XModel[NXModel-1];

    //declarations CRN
    vector<double> XPosCRN(NCRNData);
    vector<double> NModel(NCRNData);
    vector<double> ResidualsCRN(NCRNData); 
    
    //Interpolate to sample locations
    for (int i=0; i<NCRNData; ++i)
    {
        //Normalising CRN data to modelled cliff position (CRN data file: cliff position = 0)
        XPos = CliffPositionX - CRNXData[i];   //testing with CliffX not CliffPositionCRNX
        
        //Take X value of sample and interpolate to get model results at this point
        int j=0;
        while ((XModel[j]-XPos) < 0) ++j;
       InterpScale = (XModel[j]-XPos)/(XModel[j]-XModel[j-1]);

        //Get Interpolated N value
        CRNModel[i] = CRNModelData[j]-ScaleCRN*(CRNModelData[j]-CRNModelData[j-1]);
    }


    //calculate CRN likelihood
    for (int i=0; i<NCRNData; ++i)
    {
        CRNLikelihood *= exp(-(fabs((CRNNData[i]-CRNModel[i])*(CRNNData[i]-CRNModel[i])))/(CRNNErrorData[i]*CRNNErrorData[i]));
    }

    //need a return statement.
    //separate functions for CRN And Topo likelihoods?
    return CRNLikelihood;
}

#endif
