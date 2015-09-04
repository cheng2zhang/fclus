#ifndef HIST_H__
#define HIST_H__



/* simple histogram structure */



#include "util.h"



typedef struct {
  int n;
  double dx;
  double xmin;
  double xmax;
  double *arr;
} hist_t;



static hist_t *hist_open(double xmin, double xmax, double dx)
{
  hist_t *h;
  int i;

  xnew(h, 1);
  h->xmin = xmin;
  h->dx = dx;
  h->n = (int) ((xmax - xmin) / dx + 0.5);
  h->xmax = h->xmin + h->n * h->dx;
  xnew(h->arr, h->n);
  for ( i = 0; i < h->n; i++ )
    h->arr[i] = 0;
  return h;
}



static void hist_close(hist_t *h)
{
  free(h->arr);
  free(h);
}



static int hist_add(hist_t *h, double x)
{
  if ( x >= h->xmin ) {
    int i = ( x - h->xmin ) / h->dx;
    if ( i < h->n ) {
      h->arr[i] += 1;
      return 0;
    }
  }
  return -1;
}



static int hist_save(hist_t *h, const char *fn)
{
  int i, imin, imax;
  FILE *fp;

  /* find the minimal index */
  for ( i = 0; i < h->n; i++ ) {
    if ( h->arr[i] > 0 )
      break;
  }
  if ( i >= h->n ) return -1;
  imin = i;

  /* find the maximal index */
  for ( i = h->n - 1; i >= imin; i-- ) {
    if ( h->arr[i] > 0 )
      break;
  }
  if ( i < imin ) return -1;
  imax = i + 1;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %g %g %g %d\n", h->xmin, h->xmax, h->dx, h->n);
  for ( i = imin; i < imax; i++ ) {
    fprintf(fp, "%g %g\n", h->dx * (i + 0.5), h->arr[i]);
  }

  fprintf(stderr, "saving histogram file %s\n", fn);
  fclose(fp);
  return 0;
}



static int hist_load(hist_t *h, const char *fn)
{
  int i, n = 0;
  double xmin = 0, xmax = 0, dx = 0, x, y;
  char s[80];
  FILE *fp;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  fgets(s, sizeof s, fp);
  if ( s[0] != '#' ) {
    rewind(fp);
  } else {
    if ( sscanf(s + 1, "%lf%lf%lf%d", &xmin, &xmax, &dx, &n) != 4
      || fabs(h->xmin - xmin) > dx * 0.01
      || fabs(h->dx - dx) > dx * 0.01
      || h->n != n ) {
      fprintf(stderr, "incompatible setting, xmin %g vs %g, dx %g vs %g, n %d vs %d\n",
         h->xmin, xmin, h->dx, dx, h->n, n);
      fclose(fp);
      return -1;
    }
  }

  while ( fgets(s, sizeof s, fp) ) {
    sscanf(s, "%lf%lf", &x, &y);
    i = (int) ( (x - h->xmin) / h->dx );
    if ( i < 0 || i >= h->n ) {
      fprintf(stderr, "%s: bad x %g out of (%g, %g)\n",
          fn, x, h->xmin, h->xmax);
    }
    h->arr[i] = y;
  }

  fclose(fp);
  return 0;
}



#endif /* HIST_H__ */
