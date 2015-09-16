#ifndef CAGONC_H__
#define CAGONC_H__



#include "cagocore.h"
#include "hmc.h"
#include "wl.h"



#define VMAX DBL_MAX



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
__inline static int cago_metro_nc(cago_t *go, wl_t *wl,
    double amp, double bet, int *pnc)
{
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
  unc = wl_getvi(wl, nc);
  dunc = unc - wl_getvi(wl, *pnc);

  dutot = bet * du + dunc;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dutot ) );
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
__inline static int cago_hmc_nc(cago_t *go, wl_t *wl,
    hmc_t *hmc, int *pnc)
{
  /* compute the current RMSD */
  int nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  int acc;
  double dv = 0;

  dv  = wl_getvi(wl, nc);
  dv -= wl_getvi(wl, hmc->idat[0]);

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





/* initialize an HMC object for implicit NC */
__inline hmc_t *cago_ihmc_nc_init(cago_t *go, int *idat, double **pfdat)
{
  hmc_t *hmc;
  double *fdat;
  int cnt = 1 + go->n * D;

  /* make a hybrid Monte-Carlo object
   * with 1 extra integer
   * and n*D + 1 extra floating-point numbers:
   * NC, potential energy, and the fit structure */
  hmc = hmc_open(go->n, 1, cnt);

  xnew(fdat, cnt);
  idat[0] = cago_ncontacts(go, go->x, -1, NULL, NULL);
  fdat[0] = cago_rmsd(go, go->x, go->x1);
  memcpy(fdat + 1, go->x1, go->n * D * sizeof(double));
  *pfdat = fdat;

  /* push the initial state */
  hmc_push(hmc, go->x, go->v, go->f, idat, fdat);

  return hmc;
}




#endif /* CAGONC_H__ */

