// MCMC.cpp
#include <cstdio>
#include "./SL.hpp"
#include "./RP.hpp"
#include "./MCMC.hpp"

#ifndef MCMC_CPP
#define MCMC_CPP

void MCMC::Initialise() {}

void MCMC::Initialise(int InputNoIterations)
{
    NoIterations = InputNoIterations;
}

void MCMC::Initialise(int InputNoIterations, double InputSeaLevel, double InputGradient)
{
    NoIterations = InputNoIterations;
    SL MySL = SL(InputSeaLevel);
    RP MyRP = RP(InputGradient);

    printf("Sea level is %1.1f m\n", MySL.get_SeaLevel());
    printf("Gradient is %1.1f m/m\n", MyRP.get_Gradient());
}

#endif