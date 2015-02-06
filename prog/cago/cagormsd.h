#ifndef CAGORMSD_H__
#define CAGORMSD_H__



#include "cagocore.h"
#include "hmc.h"
#include "wl.h"



#define VMAX DBL_MAX



typedef struct {
  cago_t *go;
  double xmin;
  double dx;
  wl_t *wl;
} cagormsd_t;



__inline static cagormsd_t *cagormsd_open(cago_t *go, wl_t *wl,
    double xmin, double dx)
{
  cagormsd_t *r;

  xnew(r, 1);
  r->go = go;
  r->wl = wl;
  r->xmin = xmin;
  r->dx = dx;
  return r;
}



__inline static void cagormsd_close(cagormsd_t *r)
{
  free(r);
}



__inline static int cagormsd_id(cagormsd_t *r, double x)
{
  wl_t *wl = r->wl;
  int i;

  if ( x < r->xmin ) return -1;
  i = (int) ( (x - r->xmin) / r->dx );
  return ( i < wl->n ) ? i : -1;
}



__inline static double cagormsd_getv(cagormsd_t *r, double x)
{
  int i = cagormsd_id(r, x);
  return ( i < 0 ) ? VMAX : r->wl->v[i];
}



/* a step of HMC
 * `*rmsd` gives the current RMSD on return */
__inline static int cagormsd_hmc(cagormsd_t *r, hmc_t *hmc, double *rmsd1)
{
  cago_t *go = r->go;
  /* compute the current RMSD */
  double rmsd = cago_rmsd(go, go->x, NULL);
  int acc;
  double fdat[2];
  double dv = 0;

  dv  = cagormsd_getv(r, rmsd);
  dv -= cagormsd_getv(r, hmc->fdat[0]);

  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double rr = rand01();
    acc = ( rr < exp( -dv ) );
  }
  if ( acc ) {
    fdat[0] = rmsd;
    fdat[1] = go->epot;
    hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);
  } else {
    hmc_pop(hmc, go->x, go->v, go->f, NULL, fdat, 1);
    go->epot = fdat[1];
  }
  *rmsd1 = fdat[0];
  return acc;
}



/* randomly swap the velocities of k pairs of particles */
#define cago_vscramble(go, v, k) md_vscramble(v, NULL, go->n, k)



#endif /* CAGORMSD_H__ */

