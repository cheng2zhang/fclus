#ifndef CAGORMSD_H__
#define CAGORMSD_H__



#include "cagocore.h"
#include "hmc.h"
#include "wl.h"



#define VMAX DBL_MAX



typedef struct {
  cago_t *go;
  wl_t *wl;
} cagonc_t;



__inline static cagonc_t *cagonc_open(cago_t *go, wl_t *wl)
{
  cagonc_t *r;

  xnew(r, 1);
  r->go = go;
  r->wl = wl;
  return r;
}



__inline static void cagonc_close(cagonc_t *r)
{
  free(r);
}



__inline static double cagonc_getv(cagonc_t *r, int x)
{
  return wl_getvi(r->wl, x);
}



/* compute the number of contacts with a moved particle */
__inline static int cago_ncontacts2(cago_t *go,
  double (*x)[D], int i, double *xi, double gam,
  double *Q, int *mat)
{
  int j, n = go->n;

  for ( j = 0; j < n; j++ ) {
    vcopy(go->x1[j], x[j]);
  }
  vcopy(go->x1[i], xi);
  return cago_ncontacts(go, go->x1, gam, Q, mat);
}



/* Metropolis algorithm */
__inline static int cagonc_metro(cagonc_t *r,
    double amp, double bet, int *pnc)
{
  cago_t *go = r->go;
  int i, acc;
  double xi[D], du, dutot;
  int nc;
  double unc, dunc;

  i = (int) (go->n * rand01());
  xi[0] = amp * (rand01() * 2 - 1);
  xi[1] = amp * (rand01() * 2 - 1);
  xi[2] = amp * (rand01() * 2 - 1);
  vinc(xi, go->x[i]);
  du = cago_depot(go, go->x, i, xi);

  nc = cago_ncontacts2(go, go->x, i, xi, -1, NULL, NULL);
  unc = cagonc_getv(r, nc);
  dunc = unc - cagonc_getv(r, *pnc);

  dutot = bet * du + dunc;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    double rr = rand01();
    acc = ( rr < exp( -dutot ) );
  }
  if ( acc ) {
    vcopy(go->x[i], xi);
    go->epot += du;
    *pnc = nc;
    //printf("%d, %d\n", nc, cago_ncontacts(go, go->x, -1, NULL, NULL));
    //getchar();
    return 1;
  } else {
    return 0;
  }
}



/* a step of HMC
 * `*pnc` gives the current number of contact on return */
__inline static int cagonc_hmc(cagonc_t *r, hmc_t *hmc, int *pnc)
{
  cago_t *go = r->go;
  /* compute the current RMSD */
  int nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  int acc;
  double dv = 0;

  dv  = cagonc_getv(r, nc);
  dv -= cagonc_getv(r, hmc->idat[0]);

  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double rr = rand01();
    acc = ( rr < exp( -dv ) );
  }
  if ( acc ) {
    hmc_push(hmc, go->x, go->v, go->f, &nc, &go->epot);
  } else {
    hmc_pop(hmc, go->x, go->v, go->f, &nc, &go->epot, 1);
  }
  *pnc = nc;
  return acc;
}



/* randomly swap the velocities of k pairs of particles */
#define cago_vscramble(go, v, k) md_vscramble(v, NULL, go->n, k)



#endif /* CAGORMSD_H__ */

