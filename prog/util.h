#ifndef UTIL_H__
#define UTIL_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "vct.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif


/* update the MC move size according to the acceptance ratio */
__inline static void update_mcamp(double *amp, double acc,
    double *nacc, double *ntot)
{
  double x = sqrt( *nacc / *ntot / acc );
  if ( x > 2 ) x = 2;
  else if ( x < 0.5 ) x = 0.5;
  *amp *= x;
  fprintf(stderr, "acc %g%%, amp %g\n", 100*(*nacc)/(*ntot), *amp);
  *nacc = 0;
  *ntot = DBL_MIN;
}


#endif
