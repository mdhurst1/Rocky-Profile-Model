/* 

Fast implementation of approximation to the exponential function taken from

N. N. Schraudolph. "A fast, compact approximation of the exponential function." 
Neural Computation, 11(4), May 1999, pp.853-862.

MDH, July 2018

*/

#include <cmath>
#include <iostream>


using namespace std;

static union 
{
  double d;
  struct {
#ifdef LITTLE_ENDIAN
    int j,i;
#else 
    int i,j;
#endif
  } n;
} _eco;

#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define fastexp(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)