#ifndef UTIL_H__
#define UTIL_H__



/* utilities */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



#ifndef xrenew
#define xrenew(x, n) { \
  if ( (x = realloc(x, sizeof(*(x)) * (n))) == NULL ) { \
    fprintf(stderr, "no memory for " #x " x %d\n", (int) (n)); \
    exit(1); } }
#endif



/* return the larger of x and y */
__inline static int intmax(int x, int y)
{
  return x > y ? x : y;
}



/* return the smaller of x and y */
__inline static int intmin(int x, int y)
{
  return x < y ? x : y;
}



/* return the larger of a and b */
__inline double dblmax(double a, double b)
{
  return a > b ? a : b;
}



/* return the smaller of a and b */
__inline double dblmin(double a, double b)
{
  return a < b ? a : b;
}



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



/* String routines */



/* remove leading and trailing spaces */
__inline static char *strstrip(char *s)
{
  char *p, *q;

  /* remove trailing spaces */
  for ( p = s + strlen(s) - 1; isspace(*p); p-- )
    *p = '\0';

  /* remove leading spaces */
  for ( p = s; *p && isspace(*p); p++ )  ;
  if ( p != s )
    for ( q = s; (*q++ = *p++) != '\0'; ) ;
  return s;
}



#define strcmpfuzzy(s, t) strncmpfuzzy(s, t, INT_MAX)

/* comparison, ignoring cases, spaces and punctuations */
__inline static int strncmpfuzzy(const char *s, const char *t, int n)
{
  int is, it, i;
  const char cset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789()[]{}";

  for ( i = 0; i < n; s++, t++, i++ ) {
    while ( *s != '\0' && strchr(cset, *s) == NULL ) s++;
    while ( *t != '\0' && strchr(cset, *t) == NULL ) t++;
    is = tolower(*s);
    it = tolower(*t);
    if ( is != it ) return is - it;
    if ( *s == '\0' ) return 0;
  }
  return 0;
}



#endif
