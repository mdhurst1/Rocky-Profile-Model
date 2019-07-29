// MCMC.hpp
#include "./SL.hpp"
#include "./RP.hpp"

#ifndef MCMC_HPP
#define MCMC_HPP

class MCMC
{
    private:
        //internal variables
        int NoIterations;
        SL MCMCSeaLevel;
        RP MCMCRockProfile;

        // private initialise function
        void Initialise();
        void Initialise(int InputNoIterations);
        void Initialise(int InputNoIterations, double InputSeaLevel, double InputGradient);

    public:
        // these are called constructors since it is a call to itself
        // we can have empty constructors or constructors with arguments
        MCMC(){Initialise();}
        MCMC(int InputNoIterations)
        {
            Initialise(InputNoIterations);
        }
        MCMC(int InputNoIterations, double InputSeaLevel, double InputGradient) 
        {
            Initialise(InputNoIterations, InputSeaLevel, InputGradient);
        }

        // function to get a value from inside the object
        int get_NoIterations() {return NoIterations;}
};

#endif