#ifndef XYMAP_H__
#define XYMAP_H__



/* x-y mapping object */



#include "util.h"



/* object for computing average along a quantity x */
typedef struct {
  double xmin, xmax, dx;
  double dydxmin, dydxmax;
  int n;
  double (*arr)[2];
} xymap_t;



__inline static xymap_t *xymap_open(double xmin, double xmax, double dx,
    double dydxmin, double dydxmax)
{
  xymap_t *xy;

  xnew(xy, 1);
  xy->xmin = xmin;
  xy->dx = dx;
  xy->n = (int) ((xmax - xmin) / dx);
  xy->xmax = xmin + dx * xy->n;
  xy->dydxmin = dydxmin;
  xy->dydxmax = dydxmax;
  xnew(xy->arr, xy->n);
  return xy;
}



__inline static void xymap_close(xymap_t *xy)
{
  free(xy->arr);
  free(xy);
}



/* return the index of `x`, or -1 on error */
__inline static int xymap_getid(xymap_t *xy, double x)
{
  int i;

  if ( x < xy->xmin ) return -1;
  i = (int) ((x - xy->xmin) / xy->dx);
  return ( i < xy->n ) ? i : -1;
}


/* add a pair (x, y) to the database */
__inline static int xymap_add(xymap_t *xy, double x, double y)
{
  int i = xymap_getid(xy, x);

  if ( i < 0 ) return -1;
  xy->arr[i][0] += 1;
  xy->arr[i][1] += y;
  return 0;
}



/* return the average y value at x
 * or 0 on failure */
__inline static double xymap_gety(xymap_t *xy, double x)
{
  int i = xymap_getid(xy, x);

  if ( i < 0 || xy->arr[i][0] <= 0 ) return 0;
  return xy->arr[i][1] / xy->arr[i][0];
}



/* compute the derivative dy/dx at x
 * return 1.0 on error */
__inline static double xymap_getdydx(xymap_t *xy, double x)
{
  int i, i1, i2, n = xy->n;
  double side, dx = xy->dx, y1, y2, (*arr)[2] = xy->arr;
  double dydx;

  /* first compute i1 and i2 for two nearby bins */
  i = xymap_getid(xy, x);
  if ( i < 0 ) return 1.0;

  /* see if x lies in the lower or higher half of the bin
   * in the former case, we use (i-1, i)
   * in the latter case, we use (i, i+1) */
  side = (x - (xy->xmin + i * dx)) / dx;
  if ( side < 0.5 ) { /* lower half */
    if ( i > 0 ) {
      for ( i1 = i - 1; i1 >= 0 && arr[i1][0] <= 0; i1-- ) ;
      i2 = i;
    } else {
      i1 = 0;
      i2 = 1;
    }
  } else { /* higher half */
    if ( i < n - 1 ) {
      i1 = i;
      for ( i2 = i + 1; i2 < n && arr[i1][0] <= 0; i1++ ) ;
      i2 = i + 1;
    } else {
      i1 = n - 2;
      i2 = n - 1;
    }
  }

  if ( arr[i1][0] <= 0 || arr[i2][0] <= 0 )
    return 1.0;

  /* compute the derivative */
  y1 = arr[i1][1] / arr[i1][0];
  y2 = arr[i2][1] / arr[i2][0];

  dydx = ( y2 - y1 ) / ((i2 - i1) * dx);
  if ( dydx < xy->dydxmin ) {
    dydx = xy->dydxmin;
  } else if ( dydx > xy->dydxmax ) {
    dydx = xy->dydxmax;
  }

  return dydx;
}



/* save the xy profile to file */
__inline static int xymap_save(xymap_t *xy, const char *fn)
{
  FILE *fp;
  int i;
  double s, y;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  /* write the information line */
  fprintf(fp, "# %d %g %g\n",
      xy->n, xy->xmin, xy->dx);

  /* write data */
  for ( i = 0; i < xy->n; i++ ) {
    /* each line reads:  x   y   total */
    s = xy->arr[i][0];
    y = ( s > 0 ? xy->arr[i][1] / s : 0 );
    fprintf(fp, "%g %g %g\n",
        xy->xmin + (i + 0.5) * xy->dx,
        y, s);
  }

  fclose(fp);

  return 0;
}



#endif /* XYMAP_H__ */

