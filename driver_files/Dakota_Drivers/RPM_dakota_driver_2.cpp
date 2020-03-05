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
#include "../../RPM.hpp"
#include "../../RoBoCoP_CRN/RockyCoastCRN.hpp"
#include "../../SeaLevel.hpp"
#include "../../FastExp.hpp"

using namespace std;

template <typename T> string tostr(const T& t)
{ 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}


int main(int nNumberofArgs,char *argv[])
{
	cout << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "|  Rocky Profile Model (RPM)                                                     |" << endl;
	cout << "|  This program models the development of shore platforms                        |" << endl;
	cout << "|  following model developed by Matsumoto et al. (2016)                          |" << endl;
	cout << "|                                                                                |" << endl;
	cout << "|  Implemented in C++ by Martin Hurst, University of Glasgow                     |" << endl;
	cout << "|  for coupling to RockyCoastCRN; model for predicting                           |" << endl;
	cout << "|  cosmogenic radionuclide concentrations in shore platforms                     |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << endl;

	//Test for correct input arguments
	if (nNumberofArgs!=12)
	{
		cout << "Error: This program requires 10 (YES TEN, one-zero) command line inputs: " << endl;
		cout << " * First a path to the folder where the model will be run" << endl;
		cout << " * The name of the project/model run" << endl;
        cout << " * The name of the topo profile data file" << endl;
        cout << " * The name of the CRN conc data file" << endl;
		cout << " * A Flag to run with CRNs (1 = True)" << endl;
        cout << " * The initial topographic gradient" << endl;
        cout << " * The tidal range (m)" << endl;
        cout << " * The subtidal weathering efficacy (multiplier)" << endl;
        cout << " * The wave attenuation constant" << endl;
        cout << " * The rock resistance (kg/m2)" << endl;
        cout << " * The Maximum weathering rate (kg/m2/yr)" << endl;
        cout << " * " << endl;
		cout << "-----------------------------------------------------------------------------" << endl;
		cout << "Then the command line argument will be: " << endl;
		cout << "In linux:" << endl;
		cout << "  ./RPM_dakota.out /ProjectFolder/ RPM_dakota_test CB_profile.txt CB_CRN.data 1 1. 4. 0.005 0.01 1000. 10." << endl;
		cout << "-----------------------------------------------------------------------------" << endl;
        cout << endl;
		exit(EXIT_SUCCESS);
	}

    // read parameters from command line arguments
	string Folder = argv[1];
	char* DakotaFilename = argv[2];
    char* ProfileDatafile = argv[3];
    char* CRNDatafile = argv[4];
	int CRNFlag = atoi(argv[5]);
	double Gradient = atof(argv[6]);
	double TidalRange = atof(argv[7]);
    double SubtidalEfficacy = atof(argv[8]);

	//Free parameters

    double WaveAttenuationConst = (atof(argv[9]));
	//double WaveAttenuationConst = pow(10,(atof(argv[9])));

    //double Resistance = pow(10,(atof(argv[10])));          //dakota varies FR on log scale
	double Resistance = atof(argv[10]);

    //double WeatheringRate = Resistance * pow(10,(atof(argv[11])));      //dakota varies K proportional to FR 0 - 0.5 range 
	//double WeatheringRate = pow(10,(atof(argv[11]))); 
	double WeatheringRate = atof(argv[11]);

	cout << "Resistance = " << Resistance << endl;
	cout << "WeatheringRate = " << WeatheringRate << endl;
	cout << "WaveAttenuationConstant = " << WaveAttenuationConst << endl;
	

	string Res = to_string(Resistance);
	string WRate = to_string(WeatheringRate);

    //initialisation parameters
	double dZ = 0.1;
	double dX = 0.1;
	double CliffHeight = 15.;
	double MinElevation = -15.;

	//Time control parameters
	//Time runs in yrs bp
	double EndTime = 0.;
	double Time = 8000.;
	double TimeInterval = 1.;

	//Print Control
	double PrintInterval = 100;
	double PrintTime = Time;

    //set up output file - used for visual when testing 
	//whwen using dakota use FR arguments etc. 
	string OutputFileName = Folder+DakotaFilename+"_ShoreProfile.xz";
	string OutputConcentrationFileName = Folder+DakotaFilename+"Concentrations.xn";
	

    // initialise sea level here and calculate MinElevation based on lowest sea level
	// Initialise Sea level from datafile
	string RelativeSeaLevelFile = "CB_RSL.data";
	SeaLevel RelativeSeaLevel = SeaLevel(RelativeSeaLevelFile);
	
	// Get initial sea level
	double InstantSeaLevel = RelativeSeaLevel.get_SeaLevel(Time);

	//MinElevation calculated from InitialRSL
	 if (MinElevation >= InstantSeaLevel)
	 { 
		MinElevation = round((InstantSeaLevel-10.)*10)/10;

		//dz needs to be rounded to 0.1 
	 }

    //initialise RPM Model
	RPM PlatformModel = RPM(dZ, dX, Gradient, CliffHeight, MinElevation);

    // Initialise sea level
	PlatformModel.UpdateSeaLevel(InstantSeaLevel);

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
	PlatformModel.InitialiseGeology(CliffHeight, CliffFailureDepth, Resistance, WeatheringRate, SubtidalEfficacy);

    // print initial condition to file - this is for testing - remove
	double TempTime = -9999;
    PlatformModel.WriteProfile(OutputFileName, TempTime);			
	if (CRNFlag) PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, TempTime);

    //Loop through time
	while (Time >= EndTime)
	{
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
			PlatformModel.WriteProfile(OutputFileName, Time);  //This is for testing - need to remove
            if (CRNFlag) PlatformCRN.WriteCRNProfile(OutputConcentrationFileName, Time);
			PrintTime -= PrintInterval;
		}

		//update time
		Time -= TimeInterval;
	}
    
    //declarations modelled results 
    vector<double> XModel, ZModel, XDataModel, CRNConcModel;

    //get morphology from model 
    XModel = PlatformModel.get_X();
    ZModel = PlatformModel.get_Elevations();
    int XSize = XModel.size();
    double CliffPositionX = XModel[XSize-1];

	cout << "CliffPositionX = " << CliffPositionX << endl;

    //get CRN concentrations from model 
    XDataModel = PlatformCRN.get_X();
	CRNConcModel = PlatformCRN.get_SurfaceN()[0];
    //int XCRNSize = XDataModel.size();
    //double CliffPositionCRNX = XDataModel[XCRNSize-1];

    //Vectors to hold extracted profile data
    int NProfileData;
    vector<double> ProfileXData;
    vector<double> ProfileZData;
    //double ZStd = 1.;   //where define ZStd?

    //Vectors to hold CRN concentration data
    int NData;
    vector<double> XData;
    vector<double> CRNConcData;
    vector<double> CRNConcErrorData;

   //////////////////////////////////// 
   //                                //
   //Read in topographic profile data// 
   //                                // 
   ////////////////////////////////////


   //Declare temp variables 
   char Dummy[32];
   float TempProfileXData, TempProfileZData;

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
   // get size of the profile data vectors
   NProfileData = ProfileXData.size();


   ////////////////////////////////// 
   //                              //
   //Read in CRN concentration data//
   //                              //
   //////////////////////////////////


   //Declare temp variables
   float TempXData, TempCRNConcData, TempCRNConcErrorData;
  
   //Generate input filestream and read data into vectors

   
	ifstream ReadCRNDataFile(CRNDatafile);
	if (!ReadCRNDataFile)
	{
	  printf("MCMC_Coast::%s line %d: Input CRN data file \"%s\" doesn't exist\n\n", __func__, __LINE__, CRNDatafile);
	  exit(EXIT_SUCCESS);
	}
	
    // ignore header lines by reading to Dummy
    // ignore first column - CRN sample name 
    // file format is...
    // Sample_header | X_header | CRN_header | Error_header |
    //   CB[0]       |   X[0]   |   CRN[0]   |     e[0]     |
    //   CB[1]       |   X[1]   |   CRN[1]   |     e[1]     |
    //  CB[...]      |  X[...]  |  CRN[...]  |    e[...]    |
    //   CB[n]       |   X[n]   |   CRN[n]   |     e[n]     |

	ReadCRNDataFile >> Dummy >> Dummy >> Dummy >> Dummy;
	while(ReadCRNDataFile >> Dummy >> TempXData >> TempCRNConcData >> TempCRNConcErrorData)
	{
    XData.push_back(TempXData);
    CRNConcData.push_back(TempCRNConcData);
    CRNConcErrorData.push_back(TempCRNConcErrorData);
	}

    NData = XData.size();

   ///////////////////////////////
   //                           //
   //Calculations for topography//
   //                           //
   ///////////////////////////////

   //declarations
   vector<double> XPos(NProfileData);
   vector<double> TopoData(NProfileData);
   vector<double> DiffX(NProfileData);
   double RMSE;
   double Scale;
   bool FailFlag = false;
   

   //Interpolate to extracted morphology X positions
   for (int i=0; i<NProfileData; ++i)
   {
       //Normalising profile data to modelled cliff position (using Swath profile data where cliff position of measured data = 0)
       XPos[i] = CliffPositionX - ProfileXData[i];

	   //if statement XPos [i] < 0 fail flag?
	   if (XPos[i]<0)
	   {
		   FailFlag = true;
		   break;
	   }

       //Take X value of extracted morph position and interpolate to get model results at this point
       int j=0;
       while ((XModel[j]- XPos[i]) <0) ++j; //XPos starts at point nearest cliff and works offshore - starts at 40m from cliff so will never =0
	   
       DiffX[i] = XModel[j] - XPos[i];
       Scale = DiffX[i]/(XModel[j]-XModel[j-1]);

       //Get Interpolated Z value
       TopoData[i] = ZModel[j]-Scale*(ZModel[j]-ZModel[j-1]);
	   
   }

   //Calculate Residuals and likelihood
   //Fail flag for inf numbers (topo profile)
   //Calculating max and min Residual for topo data 
   //Declarations for normalised residuals 

   
   double TotalResiduals = 0;
   double TotalNResiduals = 0;
   vector<double> Residuals(NProfileData); 
   vector<double> NResiduals(NProfileData);
   vector<double> LResiduals(NProfileData);
   long double Likelihood = 1.L;
   double ZStd = 0.1;
   //double MaxTopo = Residuals[0];
   //double MinTopo = Residuals[0];
   

   //standardise topo residuals
   for (int i=0; i<NProfileData; ++i)
   {
       Residuals[i] = fabs(ProfileZData[i]-TopoData[i]);
	   TotalResiduals += pow(ProfileZData[i]-TopoData[i],2);
	   
	   //Fail Flag

	   if (isinf(TotalResiduals))
	   {
		   FailFlag = true;
		   break;
	   }

	   //MinTopo = *min_element(begin(Residuals), end(Residuals));
	   //MaxTopo = *max_element(begin(Residuals), end(Residuals));  
    }
	
	for (int i=0; i<NProfileData; ++i)
	{
		//Residuals calc for Likelihood
	   LResiduals[i] = pow(ProfileZData[i]-TopoData[i],2);
	   Likelihood *= fastexp(-(fabs(LResiduals[i])/(ZStd*ZStd)));     //(ZStd*ZStd));
	}
	//return Likelihood;

	for (int i=0; i<NProfileData; ++i)
	{
	   //Feature scaling - min-max normalisation (distribution between 0 and 1)
	   //NResiduals[i] = (Residuals[i]-MinTopo)/(MaxTopo-MinTopo);

	   //normalise topo to tidalrange 
	   NResiduals[i] = (Residuals[i]/TidalRange);

	   //Likelihood *= exp(-(fabs(LResiduals[i]))/(ZStd*ZStd));
	   
	   TotalNResiduals += pow(NResiduals[i],2);
	   //Likelihood *= exp(-(fabs(Residuals[i]))/(ZStd*ZStd));    //ZStd read in from parameter file?
	}

   
   //cout << " MaxTopo = " << setprecision(10) << MaxTopo << endl;
   //cout << " MinTopo = " << setprecision(10) << MinTopo << endl;
   cout << " Total residuals Topo = " << TotalResiduals << endl;
   //cout << " TotalNResiduals = " << TotalNResiduals << endl;
   

   cout << " Likelihood = " << scientific << Likelihood << endl;
   

   ///////////////////////////////////////                                
   //                                   //
   //Calculations for CRN concentrations//
   //                                   //
   ///////////////////////////////////////

   //declarations CRN
   vector<double> DiffCRNX(NData);
   vector<double> XPosCRN(NData);
   vector<double> NModel(NData);
   double ScaleCRN;
   //double CRN_RMSE;
   //long double LikelihoodCRN = 1.L;  


   //Interpolate to sample locations
	for (int i=0; i<NData; ++i)
	{
	    //Normalising CRN data to modelled cliff position (CRN data file: cliff position = 0)
        XPosCRN[i] = CliffPositionX - XData[i];   //testing with CliffX not CliffPositionCRNX
        
        //Take X value of sample and interpolate to get model results at this point
	    int j=0;
	    while ((XDataModel[j]-XPosCRN[i]) < 0) ++j;

	    DiffCRNX[i] = XDataModel[j]-XPosCRN[i];
        ScaleCRN = DiffCRNX[i]/(XDataModel[j]-XDataModel[j-1]);
  
        //Get Interpolated N value
        NModel[i] = CRNConcModel[j]-ScaleCRN*(CRNConcModel[j]-CRNConcModel[j-1]);
	}
	
	//Calculate likelihood
    //double TotalResidualsCRN = 0;
	vector<double> ResidualsCRN(NData); 
	//double MaxCRN = ResidualsCRN[0];
	//double MinCRN = ResidualsCRN[0];

    //Declarations for normalised Residuals 
    vector<double> NResidualsCRN(NData);
    double TotalNResidualsCRN = 0;
	double TotalResidualsCRN = 0;
	double MaxCRNCB = 16162;

	//Standardise CRN residuals
	for (int i=0; i<NData; ++i)
   {
       ResidualsCRN[i] = fabs(CRNConcData[i]-NModel[i]);

	   TotalResidualsCRN += pow(CRNConcData[i]-NModel[i],2);

	   if (isinf(TotalResidualsCRN))
	   {
		   FailFlag = true;
		   break;
	   }

	   //MinCRN = *min_element(begin(ResidualsCRN), end(ResidualsCRN));
	   //MaxCRN = *max_element(begin(ResidualsCRN), end(ResidualsCRN)); 
   }

   for (int i=0; i<NData; ++i)
   { 
	   //Feature scaling - min-max normalisation (distribution between 0 and 1)
	   //NResidualsCRN[i] = (ResidualsCRN[i]-MinCRN)/(MaxCRN-MinCRN);

	   //Normalise CRN to max measured CRN conc
	   NResidualsCRN[i] = (ResidualsCRN[i]/MaxCRNCB);

	   TotalNResidualsCRN += pow(NResidualsCRN[i],2);    

	   //TotalResidualsCRN += pow(CRNConcData[i]-NModel[i],2);
	   //LikelihoodCRN *= exp(-(fabs(ResidualsCRN[i]))/(CRNConcErrorData[i]*CRNConcErrorData[i]));
   }
	
    
   //cout << " MaxCRN = " << setprecision(10) << MaxCRN << endl;
   //cout << " MinCRN = " << setprecision(10) << MinCRN << endl;
   cout << " Total Residuals CRN = " << TotalResidualsCRN << endl;
   cout << " TotalNResidualsCRN = " << TotalNResidualsCRN << endl;

    ////////////////////////////////////
    //                                //
    //Write results to file for dakota//
    //                                //
    ////////////////////////////////////
   
   //Output residuals/ likelihood to file 
   ofstream outfile;
   outfile.open(DakotaFilename);

   //Weightings - eqaul to 1
   //double TopoWeighting = 0.5;
   //double CRNWeighting = 0.5;
   //double WeightedRMSE;
   double RMSE_N;
   //double CRN_RMSE_N;
   long double Likelihood_N = 1.L;
   long double Neg_Log_Likelihood = 1.L;
   long double Neg_Log_Likelihood_N = 1.L;

   //RMSE calculations 

   RMSE = sqrt(TotalResiduals/NProfileData);
   //CRN_RMSE = sqrt(TotalResidualsCRN/NData);

   //Normalise RMSE
   RMSE_N = RMSE/TidalRange;  // min-max rather than tidal range?
   //CRN_RMSE_N = CRN_RMSE/MaxCRNCB;

   //WeightedRMSE = (RMSE_N*TopoWeighting)+(CRN_RMSE_N*CRNWeighting);

   //Normalise Likelihood
   Likelihood_N = Likelihood/TidalRange;

   //negative log likelihood
   //if normalised correct, take -ve log of normalised likelihood
   Neg_Log_Likelihood = -log(Likelihood);

   //normalised neegative log likelihood 
   Neg_Log_Likelihood_N = Neg_Log_Likelihood/TidalRange;


   cout << " RMSE = " << RMSE << endl;
   //cout << " CRN RMSE = " << CRN_RMSE << endl;
   cout << " Normalised RMSE = " << RMSE_N << endl;
   //cout << " Normalised RMSE CRN = " << CRN_RMSE_N << endl;
   //cout << " Weighted RMSE = " << WeightedRMSE << endl; 
   //cout << " Likelihood = " << setprecision(10) << Likelihood << endl;
   cout << " Normalised Likelihood = " << scientific << Likelihood_N << endl;
   cout << " -ve log likelihood = " << scientific << Neg_Log_Likelihood << endl;
   cout << " Normalised -log likelihood = " << scientific << Neg_Log_Likelihood_N << endl;

   if (isinf(Neg_Log_Likelihood))
	   {
		   FailFlag = true;
	   }

   //add another fail flag for nan?
   
   //Check outfile is open

	if (outfile)
	{
		cout << "Filestream open" << endl;
	}
	else
	{ 
		cout << "Filestream failed to open" << endl;   
	}
    
	//Write to outfile

	if (FailFlag)
	{
		outfile << "FAIL" << endl;
	}
	else if (!CRNFlag)
	{
		outfile << Neg_Log_Likelihood_N << endl;
	}
	else
	{
	    outfile << Neg_Log_Likelihood_N << endl;
	}


	outfile.close();

	cout << endl << "Done!" << endl << endl;

}
