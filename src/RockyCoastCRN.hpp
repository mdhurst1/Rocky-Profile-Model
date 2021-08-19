/*------------------------------------------------------------------------

	PlatformCRN.hpp
	
	Object to evolve the CRN concentration across a coastal profile as a 
	function of tides, topographic shielding, relative sea level rise, 
	cliff retreat and platform downwear. For comparison with measured
	10Be CRN measurements across a coastal platform.
	
	Martin D. Hurst, British Geological Survey, December 2014

	Copyright (C) 2015, Martin Hurst

  Developer can be contacted:
  mhurst@bgs.ac.uk
  
  Martin D. Hurst
  British Geological Survey,
  Environmental Science Centre,
  Nicker Hill,
  Keyworth,
  Nottingham,
  UK,
  NG12 5GG
  	
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

/** @file RockyCoastCRN.hpp
@author Martin D. Hurst, British Geological Survey

@version Version 0.0.1
@brief RockyCoastCRN object for predicting CRN concentrations on a coastal platform
@details This object contains a simple geomorphic model for the evolution of a rocky
  coastal profile and routines to predict the accumulation of cosmogenic radionuclides
  (10Be only at the moment) in the coastal platform. Model is intended to allow inversion
  to determine rates of cliff retreat from the distribution of CRNs by linkage with the
  MCMC_RockyCoast object which calls it.
*/

/**
@mainpage
This is the documentation for the RockyCoastCRN model
These pages describe the software.
@author Martin D. Hurst, British Geological Survey

------------------------------------------------------------------------*/

#ifndef RockyCoastCRN_HPP
#define RockyCoastCRN_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "RoBoCoP.hpp"
#include "../RPM.hpp"

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
//TEMPLATES
/////////////////////////////////////////////////////////////////////////////////////////*/

///@brief Main coastal platform object.
class RockyCoastCRN
{
	friend class RoBoCoP;
	friend class RPM;
	
	private:
		vector<double> X;	  //cross shore distance (m)
		vector<double> Z;		//elevation (m)
		vector<double> Zx;	//elevation (m) at each X
		vector<int> SurfaceInd; //Positions in X of surface Z
		
		int NXNodes;	//Number of nodes across the coastline
		int NZNodes;	//Number of nodes depth
		
		double dX;    //Nodes spacing in cross shore (m)
		double dZ;    //Node spacing in vertical (m)
		
		//Sea level high latitude production rates
		///@brief   10Be spallation reference production rate.
		///@details Spallation (a/g/yr) calibrated 10Be production rate (Lifton et al. 2014).
		float Po_10Be_Spal;

		///@brief   14C spallation reference production rate.
		///@details Spallation (a/g/yr) calibrated 14C production rate (add ref).
		float Po_14C_Spal;

		///@brief   26Al spallation reference production rate.
		///@details Spallation (a/g/yr) calibrated 26Al production rate (add ref).
		float Po_26Al_Spal;

		///@brief   36Cl spallation reference production rate.
		///@details Spallation (a/g/yr) calibrated 36Cl production rate (add_ref).
		float Po_36Cl_Ca_Spal;
		float Po_36Cl_K_Spal;

		///@brief   10Be muogneic reference produciton rate.
		///@details Total muogenic production rate (a/g/yr) following Braucher et al. (2013).
		float Po_10Be_Muon_Fast;
		float Po_10Be_Muon_Slow;

		///@brief   14C muogneic reference produciton rate.
		///@details Total muogenic production rate (a/g/yr)
		float Po_14C_Muon_Fast;
		float Po_14C_Muon_Slow;

		///@brief   26Al muogneic reference produciton rate.
		///@details Total muogenic production rate (a/g/yr)
		float Po_26Al_Muon_Fast;
		float Po_26Al_Muon_Slow;

		///@brief   36Cl muogneic reference produciton rate.
		///@details Total muogenic production rate (a/g/yr)
		float Po_36Cl_Muon_Fast;
		float Po_36Cl_Muon_Slow;

		// Attenuation Lengths
		///@brief Spallogenic attenuation length (kg/m^2).
		float Lamb_Spal;

		///@brief Muogenic attenuation length (kg/m^2) following Braucher et al. (2013).
		float Lamb_Muon;

		// density of rock and water respectively
		float rho_r;
		float rho_w;

		//Decay length scales
		float z_rs;
		float z_ws;
		float z_rm;
		float z_wm;

		///@brief Half life of 10Be (Korschineck et al. 2010).
		float HalfLife_10Be;
		float Lambda_10Be;

		///@brief Half life of 14C (ref).
		float HalfLife_14C;
		float Lambda_14C;

		///@brief Half life of 26Al (ref).
		float HalfLife_26Al;
		float Lambda_26Al;

		///@brief Half life of 36Cl (ref).
		float HalfLife_36Cl;
		float Lambda_36Cl;

		int NoNuclides;								//How many nuclides to track
		vector<int> Nuclides;						//Which nuclides to track, labelled by atomic number
		vector<double>	Po_Spal;						//Holder for all spalation surface production rates
		vector<double> Po_Muon_Fast;				//Holder for all fast muon surface production rates
		vector<double> Po_Muon_Slow;				//Holder for all slow muon surface production rates
		vector<double> Lambda;						//Holder for all decay constants
		vector< vector<double> > SurfaceN;		//CRN surface concentrations	(a/g)
		vector< vector< vector<double> > > N;	//concentration of nuclides (a/g) as a function of position and depth and nuclide
		
		vector<double> PlatformElevation;		  //Platform Surface Elevations (m)
		vector<double> PlatformElevationOld;	//Platform Surface Elevations Old (m)
		vector<double> SurfaceElevation;		  //Surface Elevations (including beach cover) (m)
				
		vector<double> GeoMagTime;
		vector<double> GeoMagScalingFactors;  //Holds data from Lifton et al. (2014) Geomag model
		vector<double> RSLTime;
		vector<double> RSLRate;               //Holds data from Bradley et al. (2011) GIA model
	    bool RetreatRateSLRFlag;              // Flag to track if retreat rate has had to be modified due to SLR being too fast
		double PlatformWidth;				  //Width of model domain
		double PlatformDepth;				  //Depth of model domain
		double NDV;							  //No Data Value placeholder
		
		//CRN
		vector <vector <double> >P_Spal, P_Muon_Fast, P_Muon_Slow; //local production rates for spallation and muons for any nuclide
		
		//TIDES
		double TidalPeriod;
		double TidalAmplitude;
		double TidalRange;
		vector<double> TideLevels;
		vector<double> WaterLevels;
		vector<double> WaterDepths;
		double NTidalValues;
		
		//RetreatRates
		double RetreatRate1,RetreatRate2;
		double RetreatRate;
		double RetreatType;
		double ChangeTime;
		double PlatformGradient;
		double CliffHeight;
		double CliffGradient;
		double JunctionElevation;
		double CliffPositionX;        //tracks the cliff position in X (m)
		int CliffPositionInd;         //tracks the index of the cliff position in X
		double XMin, XMax, OldXMax;            //Extent of the model domain in X (m)
		double ZMin, ZMax;
		int ZTrackInd;
		
		vector<int> SampleInd;				//Indices for resampling morphology from another model
		int SamplingInterval;				//Spacing of SampleInd

		//sea level parameters
		double SeaLevel, SLRRate, SeaLevelChange;
		
		//Scaling parameters
		double GeoMagScalingFactor, TopoShieldingFactor;
		
		//Run time control parameters
		double Time, MaxTime, dt;
		
		//RetreatSyle
		int SteppedPlatform;
		double StepSize;
	
	    //Beach profile stuff
	    vector<double> BeachThickness;
	    double InitialBeachWidth; //Initial value of beach width when using a thinning beach width
	    double BeachWidth;      //Width of beach at BermHeight
	    double MeanBeachWidth;  //Mean Width of Beach when using a variable beach width
		double BermHeight;      //Height of Berm
        double BruunA;          //Beach steepness factor of bruun profile
		int BeachType;          //Style of beach evolution, 0 = fixed beach, 1 = sinusoidal beach width, 2 = thinning beach width
		double A;               //Sediment Scale Parameter in Bruun Profile (m^1/3)

	    string OutFileName;
	  
		//Initialise Function
		void Initialise();
		void Initialise(vector<int> WhichNuclides);
		void Initialise(double retreatrate, double beachwidth, int beachtype, double bermheight, double beachsteepness, double platformgradient, double cliffheight, double cliffgradient, double junctionelevation, double tidalamplitude, double slr, int steppedplatform, double stepsize, vector<int> WhichNuclides);
		void Initialise(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double beachsteepness, double platformgradient, double cliffheight, double cliffgradient, double junctionelevation, double tidalamplitude, double SLR, int steppedplatform, double stepsize, vector<int> WhichNuclides);
		
		void Initialise(RoBoCoP RoBoCoPCoast, vector<int> WhichNuclides);
		void Initialise(RPM RPMCoast, vector<int> WhichNuclides);
		
		//function to initialise production schematics
		void InitialiseNuclides(vector<int> WhichNuclides);
		
		//functions to initialise platform morphology
		void InitialisePlanarPlatformMorphology();
				
		//function to retrieve topographic shielding factor
		double GetTopographicShieldingFactor(double Distance, double CliffHeight);
		
		//retrieve geomagnetic scaling factor
		double GetGeoMagScalingFactor(double Time);
		
		//factor to retrieve a rate of sea level rise
		double GetSeaLevelRise(double Time);
		
		//function to get sinusoid beachwidth through time
		void GetSinWaveBeachWidth(double Time);
		void GetThinningBeachWidth(double Time);
		
	protected:

	public:
	
		/* PlatformCRN Initialisation functions */

		/// @brief Initialisation function for two retreat rate scenario.
		/// @param retreatrate1 First rate of cliff retreat (m/yr)
		/// @param retreatrate2 Second rate of cliff retreat (m/yr)
		/// @param changetime Time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param junctionelevation Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
	  ///	@author Martin D. Hurst 
    /// @date 14/09/2015
		RockyCoastCRN(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double beachsteepness, double platformgradient, double cliffheight, double cliffgradient, double junctionelevation, double tidalamplitude, double slr, int steppedplatform, double stepsize, vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(retreatrate1, retreatrate2, retreattype, changetime, beachwidth, beachtype, bermheight, beachsteepness, platformgradient, cliffheight, cliffgradient, junctionelevation, tidalamplitude, slr, steppedplatform, stepsize, WhichNuclides);
		}
		
		/// @brief Update function for two retreat rate scenario.
		/// @param retreatrate1 First rate of cliff retreat (m/yr)
		/// @param retreatrate2 Second rate of cliff retreat (m/yr)
		/// @param changetime Time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param junctionelevation Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
	  ///	@author Martin D. Hurst 
    /// @date 14/09/2015
		void UpdateParameters(double retreatrate1, double retreatrate2, int retreattype, double changetime, double beachwidth, int beachtype, double bermheight, double beachsteepness, double platformgradient, double cliffheight, double cliffgradient, double junctionelevation, double tidalamplitude, double slr, int steppedplatform, double stepsize, vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(retreatrate1, retreatrate2, retreattype, changetime, beachwidth, beachtype, bermheight, beachsteepness, platformgradient, cliffheight, cliffgradient, junctionelevation,  tidalamplitude, slr, steppedplatform, stepsize, WhichNuclides);
		}
		
		
		/// @brief Initialisation function for a single retreat rate scenario.
		/// @param retreatrate Rate of cliff retreat (m/yr)
		/// @param beachwidth Width of the beach (constant) blocking CRN production
		/// @param platformgradient Gradient (dz/dx) of coastal platform
		/// @param cliffheight Height of retreating cliff (constant)
		/// @param junctionelevation Elevation of platform at the cliff
		/// @param tidalamplitude Amplitude of diurnal tides
		/// @author Martin D. Hurst 
		/// @date 14/09/2015
		RockyCoastCRN(double retreatrate, double beachwidth, int beachtype, double bermheight, double beachsteepness, double platformgradient, double cliffheight, double cliffgradient, double junctionelevation, double tidalamplitude, double slr, int steppedplatform, double stepsize, vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(retreatrate, beachwidth, beachtype, bermheight, beachsteepness, platformgradient, cliffheight, cliffgradient, junctionelevation, tidalamplitude, slr, steppedplatform, stepsize, WhichNuclides);
		}
		
		/// @brief Initialisation function with friend class RoBoCoP as the morphological model
		/// @param RoBoCoP RoBoCoPCoast a RoBoCoP coastal morphology object 
		/// @author Martin D. Hurst 
    	/// @date 8/3/2016
		RockyCoastCRN(RoBoCoP RoBoCoPCoast, vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(RoBoCoPCoast, WhichNuclides);
		}
		
		/// @brief Initialisation function with friend class RPM as the morphological model
		/// @param RPM RPMCoast a RPM coastal morphology object 
		/// @author Martin D. Hurst 
		/// @date 13/3/2017
		RockyCoastCRN(RPM RPMCoast, vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(RPMCoast, WhichNuclides);
		}
		
		/// @brief Empty initialisation function, will throw an error.
		/// @author Martin D. Hurst 
    	/// @date 14/09/2015
		RockyCoastCRN(vector<int> WhichNuclides)
		{
			RockyCoastCRN::Initialise(WhichNuclides);
		}
		
		/// @brief Empty initialisation function, will throw an error.
		///	@author Martin D. Hurst 
    	/// @date 14/09/2015
		RockyCoastCRN()
		{
			RockyCoastCRN::Initialise();
		}
		
		/// @brief function to initialise the tides
		/// @details This function initialises tides as a single cosine wave with fixed amplitude and
		///   period. 
		/// @param A Tidal Amplitude (metres)
		/// @param T Tidal Period (hours)
		///	@author Martin D. Hurst 
    /// @date 15/03/2016
		void InitialiseTides(double A, double T);
				
		/// @brief Update RockyCoastCRN object parameters for a new model run
		/// @param RetreatRate1_Test New first rate of cliff retreat (m/yr)
		/// @param RetreatRate2_Test second rate of cliff retreat (m/yr)
		/// @param ChangeTime_Test New time to switch from retreatrate1 to retreatrate2 (years BP)
		/// @param BeachWidth_Test New width of the beach (constant) blocking CRN production
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		void UpdateParameters(double RetreatRate1_Test, double RetreatRate2_Test, double ChangeTime_Test, double BeachWidth_Test);
		
		/// @brief Launch the main program loop to evolve coast and predict CRN concentrations
		/// @details This function evolves a rocky coastal platform through time assumming gradual
		///   cliff retreat and platform downwear proportional to cliff retreat (i.e. equillibrium
		///   cliff retreat; translating a constant crossshore profile landward through time). The
		///   model takes geomagnetic scaling as a text file (generated by Lifton et al. (2014) code)
		///   and a relative sea level history text file (generated from Bradley et al. (2011) model).
		/// @param RetreatType Single (=0), step change (=1), or gradual change (=2) retreat rate scenario
		/// @param WriteResultsFlag flag to write results to file (=1) or not (=0), default is on.
        /// @return boolean flag to signal RetreatRate has been overwritten
		///	@author Martin D. Hurst 
        /// @date 14/09/2015
		bool RunModel(string OutFileName, int WriteResultsFlag=1);
		
		/// @brief Determine current retreat rate
		/// @details This function determines the cliff retreat rate for the current timestep based on
		///   the selected scenario/mode of cliff retreat. Retreat type = 0 is constant retreat rate
		///   1 = step change in retreat rate at time ChangeTime and 2 = gradual change in reteat rate
		///   through time.
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void GetRetreatRate();
		
		/// @brief Updates the CRN concentrations at the platform surface and at depth
		/// @details This function calculates the accumulation of 10Be in the platform surface and
		///   at depth by both spallation and muogenic production. 
		/// @param TimeInterval the time step in years
		/// @author Martin D. Hurst 
		/// @date 09/02/2016
		void UpdateCRNs();
		void ParallelUpdateCRNs();

		/// @brief Updates the platform morphology
		/// @details This function calculates the amount of platform downwear and updates the elevations
		///   of the platform surface. Currently just does gradual uniform downwear or step-retreat.
		/// @author Martin D. Hurst 
		/// @date 09/02/2016
		void UpdateEquillibriumMorphology();
		
		/// @brief Updates the platform morphology
		/// @details This function calculates the amount of platform downwear and cliff retre updates the elevations
		///   of the platform surface based on an iteration of the RoBoCoP coast object.
		/// @author Martin D. Hurst
		/// @param RoBoCoPCoast A RoBoCoP Coastal model object
    	/// @date 14/03/2016
		void UpdateMorphology(RoBoCoP RoBoCoPCoast);
		
		/// @brief Updates the platform morphology
		/// @details This function calculates the amount of platform downwear and 
		/// 	cliff retre updates the elevations of the platform surface based on 
		/// 	an iteration of the RPM coast object.
		/// @author Martin D. Hurst
		/// @param RPMCoast A RPM Coastal model object
    	/// @date 13/03/2017
		void UpdateMorphology(RPM RPMCoast);
		
		/// @brief Updates the platform morphology based on recieving X and Z vectors
		/// @details This function updates the platform morphology based on a profile passed through
		/// vectors X and Z. These can be provided from model output from an alternative model such
		/// as that of Matsumoto et al. (2016)
		/// @author Martin D. Hurst
		/// @param X vector<double> of horizontal position
		/// @param Z vector<double> of platform elevation
    	/// @date 26/02/2017
		void UpdateMorphology(vector<double> X, vector<double> Z);
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the elevations of the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		///	@author Martin D. Hurst 
    /// @date 09/02/2016
		void WriteProfile(string OutputFileName, double Time);
		
		/// @brief Writes the platform morphology to file
		/// @details This function writes the CRN concnetrations on the platform surface at the current time to
		///   a file. If the file exists, this is appended.
		/// @author Martin D. Hurst 
		/// @date 15/03/2016
		void WriteCRNProfile(string OutputFileName, double Time);
		
		/// @brief Writes the array of concentrations N to file
		/// @details This function writes the CRN concnetrations at the current time to
		///   a file.
		/// @author Martin D. Hurst 
		/// @date 16/03/2017
		void WriteNuclideArray(string OutputFileName, double Time, int Nuclide);
		
		/// @brief Get X coordinates
		/// @return X coordinates
		///	@author Martin D. Hurst 
    /// @date 14/09/2015
		vector<double> get_X() const { return X; }
		
		/// @brief Get surface CRN concentration
		/// @return Surface CRN concentration
		///	@author Martin D. Hurst 
		/// @date 14/09/2015
		vector< vector<double> > get_SurfaceN() const { return SurfaceN; }
		
		/// @brief Initialise Relative Sea Level history fro ma file
		///	@author Martin D. Hurst 
		/// @date 14/06/2018
		void InitialiseRSLData(string RSLFilename);
        
        /// @brief Initialise Scaling timeseries from file
		///	@author Martin D. Hurst 
		/// @date 14/11/2018
		void InitialiseScalingData(string ScalingFilename);
};

#endif
