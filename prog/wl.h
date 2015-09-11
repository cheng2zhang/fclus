#ifndef WL_H__
#define WL_H__



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>





#define WL_FLATNESSABS 0x0010 /* use (hmax-hmin)/(hmax+hmin) for flatness */
#define WL_VERBOSE     0x0001

#define WL_VMAX        DBL_MAX



typedef struct {
  int n;
  int isfloat;
  int nmin;
  double xmin, xmax, dx;
  int isinvt; /* has entered the 1/t stage */
  unsigned flags;
  double *h; /* histogram */
  double *v; /* bias potential */
  double *cf; /* counts of force */
  double *sf; /* sum of force */
  double *vf; /* bias potential from integrating the mean force */
  double tot;
  double lnf0;
  double lnf; /* current updating factor lnf */
  double flatness; /* Wang-Landau threshold for the histogram flatness */
  double frac; /* Wang-Landau reduction factor for lnf */
  double c; /* c/t */
} wl_t;



__inline static wl_t *wl_open0(int n,
    double lnf0, double flatness, double frac,
    double c, unsigned flags)
{
  wl_t *wl;
  int i;

  if ( (wl = calloc(1, sizeof(*wl))) == NULL ) {
    fprintf(stderr, "no memory for WL\n");
    return NULL;
  }
  wl->n = n;
  if ( (wl->h = calloc(n, sizeof(double))) == NULL ) {
    fprintf(stderr, "no memory for the WL histogram\n");
    free(wl);
    return NULL;
  }
  if ( (wl->v = calloc(n, sizeof(double))) == NULL ) {
    fprintf(stderr, "no memory for the WL potential\n");
    free(wl);
    return NULL;
  }
  if ( (wl->cf = calloc(n, sizeof(double))) == NULL ) {
    fprintf(stderr, "no memory for the mean-force counts\n");
    free(wl);
    return NULL;
  }
  if ( (wl->sf = calloc(n + 1, sizeof(double))) == NULL ) {
    fprintf(stderr, "no memory for the mean force\n");
    free(wl);
    return NULL;
  }
  if ( (wl->vf = calloc(n + 1, sizeof(double))) == NULL ) {
    fprintf(stderr, "no memory for the potential of mean force\n");
    free(wl);
    return NULL;
  }
  for ( i = 0; i < n; i++ ) {
    wl->h[i] = 0.0;
    wl->v[i] = 0.0;
    wl->cf[i] = 0.0;
    wl->sf[i] = 0.0;
    wl->vf[i] = 0.0;
  }
  wl->lnf = wl->lnf0 = lnf0;
  wl->tot = 0;
  wl->isinvt = 0;
  wl->flatness = flatness;
  wl->frac = frac;
  wl->c = c;
  wl->flags = flags;
  return wl;
}



/* open a Wang-Landau object for an integer */
__inline static wl_t *wl_openi(int nmin, int nmax,
    double lnf0, double flatness, double frac,
    double c, unsigned flags)
{
  wl_t *wl;
  int n = nmax - nmin + 1;

  wl = wl_open0(n, lnf0, flatness, frac, c, flags);
  wl->nmin = nmin;
  wl->xmin = nmin;
  wl->dx = 1.0;
  wl->xmax = nmax;
  wl->isfloat = 0;
  return wl;
}



/* open a Wang-Landau object for a floating-point number */
__inline static wl_t *wl_openf(double xmin, double xmax, double dx,
    double lnf0, double flatness, double frac,
    double c, unsigned flags)
{
  wl_t *wl;
  int n;

  n = (int) ((xmax - xmin) / dx + 0.5);
  wl = wl_open0(n, lnf0, flatness, frac, c, flags);
  wl->xmin = xmin;
  wl->dx = dx;
  wl->xmax = xmin + dx * n;
  wl->nmin = 0;
  wl->isfloat = 1;
  return wl;
}



__inline static void wl_close(wl_t *wl)
{
  free(wl->h);
  free(wl->v);
  free(wl->cf);
  free(wl->sf);
  free(wl->vf);
  free(wl);
}



/* compute the total of the histogram */
__inline static double wl_gethtot(const double *h, int n)
{
  double htot = 0;
  int i;

  for ( i = 0; i < n; i++ ) {
    htot += h[i];
  }
  return htot;
}



/* clear the histogram */
__inline static void wl_clearh(double *h, int n)
{
  int i;

  for ( i = 0; i < n; i++ ) {
    h[i] = 0.0;
  }
}



/* trim the bottom of the potential */
__inline static void wl_trimv(double *v, int n)
{
  double vmin = v[0];
  int i;

  for ( i = 1; i < n; i++ ) {
    if ( v[i] < vmin ) {
      vmin = v[i];
    }
  }
  for ( i = 0; i < n; i++ ) {
    v[i] -= vmin;
  }
}



/* retrieve the bias potential */
__inline static double wl_getvi(wl_t *wl, int i)
{
  i -= wl->nmin;
  return ( i >= 0 && i < wl->n ) ? wl->v[i] : WL_VMAX;
}



/* retrieve the bias potential */
__inline static double wl_getvf(wl_t *wl, double x)
{
  if ( x < wl->xmin ) {
    return WL_VMAX;
  } else {
    int i = (int) ((x - wl->xmin) / wl->dx);
    return ( i < wl->n ) ? wl->v[i] : WL_VMAX;
  }
}



/* limit the mean force
 * if `x < wl->xmin`, limit the force in (flmin, flmax)
 * if `x > wl->xmax`, limit the force in (fhmin, fhmax) */
__inline static double wl_wrapmf(int flags, double f,
    double flmin, double flmax, double fhmin, double fhmax)
{
  if ( flags < 0 ) {
    if ( f < flmin ) {
      f = flmin;
    } else if ( f > flmax ) {
      f = flmax;
    }
  }
  if ( flags > 0 ) {
    if ( f < fhmin ) {
      f = fhmin;
    } else if ( f > fhmax ) {
      f = fhmax;
    }
  }

  return f;
}



/* retrieve the gradient dv/dx from the bias potential */
__inline static double wl_getdvdx_v(wl_t *wl, double x,
    double flmin, double flmax, double fhmin, double fhmax)
{
  int i, flags = 0;
  double f;

  if ( x < wl->xmin ) {
    i = 0;
    flags = -1;
  } else {
    i = (int) ((x - wl->xmin) / wl->dx + 0.5) - 1;
    if ( i < 0 ) {
      i = 0;
    } else if ( i >= wl->n - 2 ) {
      i = wl->n - 2;
    }
    if ( x >= wl->xmax ) {
      flags = 1;
    }
  }
  f = (wl->v[i + 1] - wl->v[i]) / wl->dx;

  return wl_wrapmf(flags, f, flmin, flmax, fhmin, fhmax);
}



__inline static double wl_getdelv_v(wl_t *wl, double x1, double x2)
{
  return wl_getvf(wl, x2) - wl_getvf(wl, x1);
}



/* compute the bin containing `x`
 * `flag` is +1 or -1 for overflow or underflow */
__inline static int wl_getid_mf(wl_t *wl, double x, int *flag)
{
  int i;

  if ( x < wl->xmin ) {
    i = 0;
    *flag = -1;
  } else {
    i = (int) ((x - wl->xmin) / wl->dx);
    if ( i >= wl->n ) {
      i = wl->n - 1;
      *flag = 1;
    } else {
      *flag = 0;
    }
  }
  return i;
}



/* retrieve the gradient dv/dx from the accumulated mean force */
__inline static double wl_getdvdx_mf(wl_t *wl, double x, double minh,
    double flmin, double flmax, double fhmin, double fhmax)
{
  int i, flag;
  double f;

  i = wl_getid_mf(wl, x, &flag);
  f = ( wl->cf[i] > minh ) ? wl->sf[i] / wl->cf[i] : 0;
  return wl_wrapmf(flag, f, flmin, flmax, fhmin, fhmax);
}



/* retrieve difference of the bias potential, v(x2) - v(x1),
 * from integrating the accumulated mean force
 * we use an approximate formula */
__inline static double wl_getdelv_mf(wl_t *wl,
    double x1, double x2, double minh,
    double flmin, double flmax, double fhmin, double fhmax)
{
  double f1 = wl_getdvdx_mf(wl, x1, minh, flmin, flmax, fhmin, fhmax);
  double f2 = wl_getdvdx_mf(wl, x2, minh, flmin, flmax, fhmin, fhmax);

  return (f1 + f2) / 2 * (x2 - x1);
}



/* integrate the mean force to get the bias potential */
__inline static void wl_intmf(wl_t *wl)
{
  int i;
  double mf;

  wl->vf[0] = 0;
  for ( i = 0; i < wl->n; i++ ) {
    mf = ( wl->cf[i] > 0 ) ? wl->sf[i] / wl->cf[i] : 0;
    wl->vf[i+1] = wl->vf[i] + mf * wl->dx;
  }
}



/* add an integral entry, update the histogram and potential */
__inline static int wl_addi(wl_t *wl, int i)
{
  i -= wl->nmin;
  if ( i < 0 || i >= wl->n ) {
    if ( wl->flags & WL_VERBOSE ) {
      fprintf(stderr, "wl: out of range i %d, n %d\n", i, wl->n);
    }
    return -1;
  }
  wl->h[i] += 1.0;
  wl->v[i] += wl->lnf;
  wl->tot += 1.0;
  return 0;
}




/* add a floating-point entry, update the histogram and potential */
__inline static int wl_addf(wl_t *wl, double x)
{
  int i;

  if ( x < wl->xmin ) {
    if ( wl->flags & WL_VERBOSE ) {
      fprintf(stderr, "wl: out of range x %g < %g\n", x, wl->xmin);
    }
    return -1;
  }
  i = (int) ((x - wl->xmin) / wl->dx);
  if ( i < 0 || i >= wl->n ) {
    if ( wl->flags & WL_VERBOSE ) {
      fprintf(stderr, "wl: out of range x %g >= %g\n", x, wl->xmin + wl->dx * wl->n);
    }
    return -1;
  }
  wl->h[i] += 1.0;
  wl->v[i] += wl->lnf;
  wl->tot += 1.0;
  return 0;
}



/* update mean force */
__inline static int wl_addforcef(wl_t *wl, double x, double f)
{
  int i;

  if ( x < wl->xmin ) {
    if ( wl->flags & WL_VERBOSE ) {
      fprintf(stderr, "wl: out of range x %g < %g\n", x, wl->xmin);
    }
    return -1;
  }
  i = (int) ((x - wl->xmin) / wl->dx);
  if ( i < 0 || i >= wl->n ) {
    if ( wl->flags & WL_VERBOSE ) {
      fprintf(stderr, "wl: out of range x %g >= %g\n", x, wl->xmin + wl->dx * wl->n);
    }
    return -1;
  }
  wl->sf[i] += f;
  wl->cf[i] += 1;
  return 0;
}



/* compute the histogram flatness from the absolute value */
__inline static double wl_getflatnessabs(const double *h, int n)
{
  double hmin, hmax;
  int i;

  hmin = hmax = h[0];
  for ( i = 1; i < n; i++ ) {
    if ( h[i] > hmax ) {
      hmax = h[i];
    } else if ( h[i] < hmin ) {
      hmin = h[i];
    }
  }
  return hmax > hmin ? (hmax - hmin) / (hmax + hmin) : 1.0;
}



/* compute the histogram flatness from the standard deviation */
__inline static double wl_getflatnessstd(const double *h, int n)
{
  double y, sh = 0, shh = 0;
  int i;

  for ( i = 0; i < n; i++ ) {
    y = h[i];
    sh += y;
    shh += y * y;
  }
  if ( sh <= 0 ) return 1.0;
  sh /= n;
  shh = shh / n - sh * sh;
  if ( shh < 0 ) shh = 0;
  return sqrt(shh) / sh;
}



__inline static double wl_getflatness(const wl_t *wl)
{
  if ( wl->flags & WL_FLATNESSABS ) {
    return wl_getflatnessabs(wl->h, wl->n);
  } else {
    return wl_getflatnessstd(wl->h, wl->n);
  }
}



/* lnf = 1/t */
__inline static double wl_lnfinvt(const wl_t *wl)
{
  return wl->c * wl->n / wl->tot;
}



/* update lnf, return 1 if the Wang-Landau stage is switched */
__inline static int wl_updatelnf(wl_t *wl)
{
  double flatness, nlnf, lnfinvt;

  if ( wl->lnf <= 0 ) {
    return 0;
  }
  if ( wl->isinvt ) {
    wl->lnf = wl_lnfinvt(wl);
    return 0;
  }

  flatness = wl_getflatness(wl);
  if ( flatness < wl->flatness ) {
    nlnf = wl->lnf * wl->frac;
    lnfinvt = wl_lnfinvt(wl);
    if ( nlnf < lnfinvt && lnfinvt < wl->lnf0 ) {
      fprintf(stderr, "changing lnf from %g to %g(1/t), flatness %g%%\n",
          wl->lnf, lnfinvt, flatness*100);
      wl->isinvt = 1;
      wl->lnf = lnfinvt;
    } else {
      fprintf(stderr, "changing lnf from %g to %g (1/t %g), flatness %g%%\n",
          wl->lnf, nlnf, lnfinvt, flatness*100);
      wl->lnf = nlnf;
      wl_clearh(wl->h, wl->n);
    }
    return 1;
  }
  return 0;
}



/* save data to file */
__inline static int wl_save(wl_t *wl, const char *fn)
{
  FILE *fp;
  int i;
  double htot, mf;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  htot = wl_gethtot(wl->h, wl->n);
  wl_trimv(wl->v, wl->n);
  /* integrate the mean force */
  wl_intmf(wl);

  fprintf(fp, "# %d %d %g %g %g %d %g\n",
      wl->isfloat, wl->n, wl->xmin, wl->dx, wl->tot,
      wl->isinvt, wl->lnf);
  for ( i = 0; i < wl->n; i++ ) {
    if ( wl->isfloat ) {
      /* floating-point version */
      fprintf(fp, "%g", wl->xmin + (i + 0.5) * wl->dx);
    } else {
      /* integer version */
      fprintf(fp, "%d", wl->nmin + i);
    }

    mf = (wl->cf[i] > 0) ? wl->sf[i] / wl->cf[i] : 0;
    fprintf(fp, " %g %g %g %g %g %g\n",
        wl->v[i], wl->h[i] / (htot * wl->dx), wl->h[i],
        (wl->vf[i] + wl->vf[i+1])/2, mf, wl->cf[i]);
  }
  fclose(fp);
  return 0;
}



/* load data from file */
__inline static int wl_load(wl_t *wl, const char *fn)
{
  FILE *fp;
  int i, n = wl->n, isfloat, next;
  int isinvt;
  double x1, x2, x, y, v, tot, xmin, dx;
  double vf, mf, cf;
  double lnf;
  char ln[64000] = "";

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  if ( fgets(ln, sizeof ln, fp) == NULL
    || ln[0] != '#'
    || sscanf(ln + 1, "%d%d%lf%lf%lf%n",
              &isfloat, &n, &xmin, &dx, &tot, &next) != 5 ) {
    fprintf(fp, "%s: bad information line!\n%s", fn, ln);
    fclose(fp);
    return -1;
  }
  if ( isfloat != wl->isfloat
    || n != wl->n
    || fabs(xmin - wl->xmin) > 1e-5
    || fabs(dx - wl->dx) > 1e-8 ) {
    fprintf(stderr, "%s: dimensions mismatch, ", fn);
    fprintf(stderr, "isfloat %d#%d, n %d#%d, xmin %g#%g, dx %g#%g\n",
        isfloat, wl->isfloat, n, wl->n, xmin, wl->xmin, dx, wl->dx);
    fclose(fp);
    return -1;
  }
  wl->tot = tot;

  if ( 2 == sscanf(ln + 1 + next, "%d%lf", &isinvt, &lnf) ) {
    wl->isinvt = isinvt;
    wl->lnf = lnf;
  }

  for ( i = 0; i < n; i++ ) {
    if ( fgets(ln, sizeof ln, fp) == NULL ) {
      fprintf(stderr, "%s: cannot read for n %d\n", fn, i);
      fclose(fp);
      return -1;
    }
    if ( wl->isfloat ) {
      x2 = wl->xmin + (i + 0.5) * wl->dx;
    } else {
      x2 = wl->nmin + i;
    }
    if ( 4 != sscanf(ln, "%lf%lf%lf%lf%n", &x1, &v, &x, &y, &next)
      || fabs(x1 - x2) > 1e-6 ) {
      fprintf(stderr, "%s: bad line %d, %g#%g\n%s", fn, i, x1, x2, ln);
      fclose(fp);
      return -1;
    }
    wl->v[i] = v;
    wl->h[i] = y;
    if ( 3 != sscanf(ln + next, "%lf%lf%lf", &vf, &mf, &cf) ) {
      wl->cf[i] = cf;
      wl->sf[i] = mf * cf;
    } else {
      wl->cf[i] = 0;
      wl->sf[i] = 0;
    }
  }
  wl_intmf(wl);
  fclose(fp);
  return 0;
}



#endif /* WL_H__ */

