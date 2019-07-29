#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
using namespace std;

int main()
{
   //Declare temp variables 
   char Dummy[32];
   float TempProfileXData, TempProfileZData;

   string ProfileDatafile = "SY_profile.txt";

   //Generate input filestream and read data into vectors
   ifstream READProfileDatafile(ProfileDatafile.c_str());

   // check for problem with file
   if (!READProfileDatafile)
   { 
       printf("MCMC_RPM::%s line %d: Input Profile data file \"%s\" doesn't exist\n\n", __func__, __LINE__, ProfileDatafile.c_str());
       exit(EXIT_SUCCESS);
   }
   else
   {
       printf("MCMC_RPM::%s line %d: Input Profile data file \"%s\"\n\n", __func__, __LINE__, ProfileDatafile.c_str());
   }
   
    // ignore header lines by reading to Dummy
    // file format is...
    // X_header | Z_header
    //   X[0]   |   Z[0]
    //   X[1]   |   Z[1]
    //  X[...]  |  Z[...]
    //   X[n]   |   Z[n]
   READProfileDatafile >> Dummy;
   cout << Dummy << endl;
   READProfileDatafile >> Dummy;
   cout << Dummy << endl;
   READProfileDatafile >> TempProfileXData;
   cout << TempProfileXData << endl;
   READProfileDatafile >> TempProfileZData;
   cout << TempProfileZData << endl;
   

}