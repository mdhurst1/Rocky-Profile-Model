#include <cstdio>
#include "SL.hpp"
#include "RP.hpp"
#include "MCMC.hpp"

int main()
{
    double InitialSeaLevel = 0;
    SL MySL = SL(InitialSeaLevel);
    printf("Sea level is %1.1f m\n", MySL.get_SeaLevel());

    double InitialGradient = 0.5;
    RP MyRP = RP(InitialGradient);
    printf("Gradient is %1.1f m/m\n", MyRP.get_Gradient());

    int NoIterations = 200000;
    MCMC MyMCMC = MCMC(NoIterations);
    printf("No of iterations is %d\n", MyMCMC.get_NoIterations());
}