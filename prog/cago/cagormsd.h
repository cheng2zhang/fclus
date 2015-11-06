#ifndef CAGORMSD_H__
#define CAGORMSD_H__



#include "cagocore.h"
#include "hmc.h"
#include "wl.h"



#define VMAX DBL_MAX



/* compute the RMSD with a moved particle */
__inline static double cago_rmsd2(cago_t *go,
  double (*x)[D], int i, double *xi)
{
  int j, n = go->n;

  for ( j = 0; j < n; j++ ) {
    vcopy(go->x1[j], x[j]);
  }
  vcopy(go->x1[i], xi);
  return cago_rmsd(go, go->x1, NULL);
}



/* Metropolis algorithm (Monte Carlo) */
__inline static int cago_metro_rmsd(cago_t *go, wl_t *wl,
    double amp, double bet, double *prmsd)
{
  int i, acc;
  double xi[D], du, dutot;
  double rmsd, urmsd, durmsd;

  i = (int) (go->n * rand01());
  xi[0] = amp * (rand01() * 2 - 1);
  xi[1] = amp * (rand01() * 2 - 1);
  xi[2] = amp * (rand01() * 2 - 1);
  vinc(xi, go->x[i]);
  du = cago_depot(go, go->x, i, xi);

  rmsd = cago_rmsd2(go, go->x, i, xi);
  urmsd = wl_getvf(wl, rmsd);
  durmsd = urmsd - wl_getvf(wl, *prmsd);

  dutot = bet * du + durmsd;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    vcopy(go->x[i], xi);
    go->epot += du;
    *prmsd = rmsd;
    //printf("%g, %g\n", rmsd, cago_rmsd(go, go->x, NULL));
    //getchar();
    return 1;
  } else {
    return 0;
  }
}



/* initialize an HMC object for RMSD */
__inline hmc_t *cago_hmc_rmsd_init(cago_t *go)
{
  hmc_t *hmc;
  double fdat[2];

  /* make a hybrid Monte-Carlo object
   * and two extra floating-point numbers:
   * RMSD and potential energy */
  hmc = hmc_open(go->n, 0, 2);
  fdat[0] = cago_rmsd(go, go->x, NULL);
  fdat[1] = go->epot;
  /* push the initial state */
  hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);

  return hmc;
}



/* a step of HMC
 * `*rmsd` gives the current RMSD on return */
__inline static int cago_hmc_rmsd(cago_t *go, wl_t *wl,
    hmc_t *hmc, double *rmsd1)
{
  /* compute the current RMSD */
  double rmsd = cago_rmsd(go, go->x, NULL);
  int acc;
  double fdat[2];
  double dv = 0;

  /* get the potential of the new (current) rmsd */
  dv  = wl_getvf(wl, rmsd);
  /* get the potential of the old (stock) rmsd */
  dv -= wl_getvf(wl, hmc->fdat[0]);

  /* decide whether to accept the new rmsd */
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dv ) );
  }

  if ( acc ) {
    /* accept the new rmsd, deposit the new state */
    fdat[0] = rmsd;
    fdat[1] = go->epot;
    hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);
  } else {
    /* recover the old state */
    hmc_pop(hmc, go->x, go->v, go->f, NULL, fdat, 1);
    go->epot = fdat[1];
  }
  *rmsd1 = fdat[0];
  return acc;
}



/* randomly swap the velocities of k pairs of particles */
#define cago_vscramble(go, v, k) md_vscramble(v, NULL, go->n, k)





/* initialize an HMC object for implicit RMSD */
__inline static hmc_t *cago_ihmc_rmsd_init(cago_t *go, double **pfdat)
{
  hmc_t *hmc;
  double *fdat;
  int cnt = 2 + go->n * D;

  /* make a hybrid Monte-Carlo object
   * and n*D + 2 extra floating-point numbers:
   * RMSD, potential energy, and the fit structure */
  hmc = hmc_open(go->n, 0, cnt);

  xnew(fdat, cnt);
  fdat[0] = cago_rmsd(go, go->x, go->x1);
  fdat[1] = go->epot;
  memcpy(fdat + 2, go->x1, go->n * D * sizeof(double));
  *pfdat = fdat;

  /* push the initial state */
  hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);

  return hmc;
}



/* compute the RMSD of two structures */
__inline static double cago_rmsd_raw(cago_t *go,
    double (*x)[D], double (*xf)[D])
{
  int i, n = go->n;
  double wtot = 0, dev = 0, dx2;

  for ( i = 0; i < n; i++ ) {
    wtot += go->m[i];
    dx2 = vdist2(x[i], xf[i]);
    dev += go->m[i] * dx2;
  }

  return sqrt( dev / wtot );
}



/* velocity Verlet with RMSD bias (implicit hybrid MC) */
__inline static int cago_vv_rmsd(cago_t *go, double fs, double dt,
    double (*xf)[D], wl_t *wl,
    double mflmin, double mflmax, double mfhmin, double mfhmax,
    double kT, double step, int nsthmc, hmc_t *hmc, double *fdat,
    double *dvtot, double *hmctot)
{
  int i, n = go->n, acc;
  double dth = 0.5 * dt * fs;
  double rmsd, rmsd0, dv0, dv1, dv, dvdx;
  double dx[D], (*xf0)[D];

  for ( i = 0; i < n; i++ ) { /* VV part 1 */
    vsinc(go->v[i], go->f[i], dth / go->m[i] );
    vsinc(go->x[i], go->v[i], dt);
  }

  /* compute the normal force */
  go->epot = cago_force(go, go->x, go->f);

  /* compute the force from the bias potential */
  /* 1. rotate and translate `go->xref` to fit `go->x`
   *    save the transformed structure in `xf` */
  rmsd = cago_rmsdref(go, go->x, xf);

  /* apply the force from the bias */
  dvdx = kT * wl_getdvdx_v(wl, rmsd,
      mflmin, mflmax, mfhmin, mfhmax);
  dvdx /= rmsd * go->mtot;
  for ( i = 0; i < n; i++ ) {
    vdiff(dx, go->x[i], xf[i]);
    vsinc(go->f[i], dx, -dvdx * go->m[i]);
  }

  for ( i = 0; i < n; i++ ) { /* VV part 2 */
    vsinc(go->v[i], go->f[i], dth / go->m[i]);
  }

  /* 2. compute the RMSD from the old reference */
  xf0 = (vct *) (hmc->fdat + 2);
  rmsd0 = cago_rmsd_raw(go, go->x, xf0);

  /* 3. compute energy caused by changing
   *    the reference structure */
  /* get the potential of the new (current) rmsd */
  dv1 = wl_getvf(wl, rmsd);
  /* get the potential of the old (stock) rmsd */
  dv0 = wl_getvf(wl, rmsd0);
  dv = dv1 - dv0;

  /* accumulate this change */
  *dvtot += dv;

  if ( fmod(step, nsthmc) < 0.5 ) {
    /* 4. decide if the change of reference is acceptable */
    if ( *dvtot <= 0 ) {
      acc = 1;
    } else {
      double r = rand01();
      acc = ( r < exp( -*dvtot ) );
    }

    if ( acc ) {
      /* accept the state */
      fdat[0] = rmsd;
      fdat[1] = go->epot;
      memcpy(fdat + 2, xf, n * D * sizeof(double));
      /* push the current state */
      hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);
    } else {
      /* pop the old state, reverse the velocity */
      hmc_pop(hmc, go->x, go->v, go->f, NULL, fdat, 1);
      rmsd = fdat[0];
      go->epot = fdat[1];
      memcpy(xf, fdat + 2, n * D * sizeof(double));
    }

    ///* if out-of-range RMSD is not allowed,
    // * vscramble is needed, see ANCHOR 1 */
    //go->ekin = cago_vscramble(go, go->v, nvswaps);

    wl_addf(wl, rmsd);
    wl_updatelnf(wl);

    /* clean up the total work */
    *dvtot = 0;
    *hmctot += 1;
  } else {
    acc = 0;
  }

  return acc;
}



#endif /* CAGORMSD_H__ */

