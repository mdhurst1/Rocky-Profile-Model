/* 
A fast implementation of the exp algorithm taken from
https://stackoverflow.com/questions/10552280/fast-exp-calculation-possible-to-improve-accuracy-without-losing-too-much-perfo
max. rel. error <= 1.73e-3 on [-87,88] */

#include <cmath> 


float fastexp (float x)
{
    volatile union
    {
        float f;
        unsigned int i;
    } cvt;

   /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
   float t = x * 1.442695041f;
   float fi = floorf(t);
   float f = t - fi;
   int i = (int)fi;
   cvt.f = (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
   cvt.i += (i << 23);                                          /* scale by 2^i */
   return cvt.f;
}


