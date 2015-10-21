#ifndef CAGONC_H__
#define CAGONC_H__



#include "cagocore.h"
#include "hmc.h"
#include "wl.h"
#include "xymap.h"



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



/* randomize the velocities
 * Here, we use a simple algorithm that only
 * randomly swap the velocities of k pairs of particles
 * In this way, the total momentum along each direction
 * is conserved. */
#define cago_vscramble(go, v, k) md_vscramble(v, NULL, go->n, k)



/* appproximate contact number */
__inline static double ncsmooth(vct a, vct b,
    double rc2, double eps, vct fa, vct fb)
{
  double dx[D], dr2, invr2, invr6, amp;

  vdiff(dx, a, b);
  dr2 = vsqr( dx );
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;

  if ( invr6 > 0.5 ) {
    /* r^6 < 2 */
    /* do nothing to the force */
    return eps;
  } else {
    /* invr6 < 0.5, r^6 > 2 */
    if ( fa ) {
      /* amp should be positive,
       * the force is the negative derivative
       * along a direction of decreasing the NC */
      amp = -eps * ( 48 * invr6 - 24 ) * invr6 / dr2;
      vsinc(fa, dx,  amp);
      vsinc(fb, dx, -amp);
    }
    return 4 * eps * (1 - invr6) * invr6;
  }
}



/* differentiable approximate number of contact
 * along with its force (negative gradient) */
__inline static double cago_forcesnc(cago_t *go, vct *x, vct *f)
{
  int i, j, id, n = go->n;
  double snc;

  for ( i = 0; i < n; i++ ) {
    vzero(f[i]);
  }

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 4; j < n; j++ ) {
      id = i * n + j;
      if ( !go->iscont[id] ) continue;
      snc += ncsmooth(x[i], x[j], go->r2ref[id],
          1.0, f[i], f[j]);
    }
  }
  return snc;
}



/* initialize an HMC object for implicit NC */
__inline static hmc_t *cago_ihmc_nc_init(cago_t *go, int *idat, double *fdat)
{
  hmc_t *hmc;

  /* make a hybrid Monte-Carlo object
   * with 1 extra integer: number of contacts
   * and 2 extra floating-point numbers:
   * bias potential energy and potential energy */
  hmc = hmc_open(go->n, 1, 2);

  idat[0] = cago_ncontacts(go, go->x, -1, NULL, NULL);
  fdat[0] = cago_forcesnc(go, go->x, go->x1);
  fdat[1] = go->epot;

  /* push the initial state */
  hmc_push(hmc, go->x, go->v, go->f, idat, fdat);

  return hmc;
}



/* velocity Verlet with NC bias (implicit hybrid MC) */
__inline static int cago_vv_nc(cago_t *go, double fs, double dt,
    wl_t *wl, xymap_t *xy,
    double mflmin, double mflmax, double mfhmin, double mfhmax,
    double kT, hmc_t *hmc, int *idat, double *fdat,
    double step, int nvswaps, int verbose)
{
  int i, n = go->n, acc;
  double dth = 0.5 * dt * fs;
  int nc1, nc0;
  double snc1, snc0;
  double nc1av, nc0av, dv0, dv1, dv;
  double dvdx1, dvdx2, dvdx;

  /* use implicit HMC, currently disabled
   * because it is not helpful */
  int usesnc = 0;

  for ( i = 0; i < n; i++ ) { /* VV part 1 */
    vsinc(go->v[i], go->f[i], dth / go->m[i] );
    vsinc(go->x[i], go->v[i], dt);
  }

  /* 1. compute the normal force */
  go->epot = cago_force(go, go->x, go->f);

  /* 2. compute the smooth number of contacts and the unscaled force,
   * the force is temporarily saved in `go->x1` */
  snc1 = cago_forcesnc(go, go->x, go->x1);

  /* retrieve the value of the stock state */
  snc0 = hmc->fdat[0];

  /* compute the current NC */
  nc1 = cago_ncontacts(go, go->x, -1, NULL, NULL);
  /* retrieve the NC of the stock state */
  nc0 = hmc->idat[0];

  /* add the NC into xymap
   * this should go before xymap_gety()
   * so we have get the current NC */
  xymap_add(xy, snc1, nc1);

  /* compute the average NC at the bias potential */
  nc1av = xymap_gety(xy, snc1);
  nc0av = xymap_gety(xy, snc0);

  /* compute the force scaling factor
   * The bias potential is V( <NC>( Q (x) ) )
   * where
   *    V(<NC>) is Wang-Landau bias potential
   *      that accepts the (average) NC as the input
   *    <NC>(Q) is the xymap object which accepts
   *      the smooth bias potential Q as input
   *      and returns the average number of contacts
   *    Q(x) is the differentiable approximate NC
   *      whose force, -dQ/dx, is currently
   *      saved in the array `go->x1`
   * So we have
   * d V( <NC>( Q(x) ) ) / dx
   * = (dV/d<NC>) (d<NC>/dQ) (dQ/dx)
   * */

  /* d V / d <NC> */
  dvdx1 = kT * wl_getdvdx_v(wl, nc1av,
      mflmin, mflmax, mfhmin, mfhmax);

  /* d <NC> / d Q */
  dvdx2 = xymap_getdydx(xy, snc1);

  dvdx = dvdx1 * dvdx2;

  /* apply the bias force */
  if ( usesnc ) {
    for ( i = 0; i < n; i++ ) {
      /* f += x1 * (-dvdx) */
      vsinc(go->f[i], go->x1[i], -dvdx);
    }
  }

  /* 3. compute the change of the WL potential
   * In implicit HMC, we need to compute the difference
   * of the difference of potential
   * The first difference is between the new and old states.
   * The second difference is between the potential at the
   * current NC and that at the average <NC> at the current
   * value of the smooth bias potential.
   * The second difference is due to that the average NC
   * is already taken care of by the smooth bias potential */

  /* get the potential of the new (current) state */
  dv1 = wl_getvi(wl, nc1);

  /* get the potential of the old (stock) state */
  dv0 = wl_getvi(wl, nc0);

  if ( usesnc ) {
    dv1 -= wl_getvi(wl, (int) (nc1av + 0.5) );
    dv0 -= wl_getvi(wl, (int) (nc0av + 0.5) );
  }

  dv = dv1 - dv0;

  /* 4. decide if the change is acceptable */
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dv ) );
  }

  if ( acc ) {
    /* accept the state */
    idat[0] = nc1;
    fdat[0] = snc1;
    fdat[1] = go->epot;

    /* push the current state */
    hmc_push(hmc, go->x, go->v, go->f, idat, fdat);

    if ( verbose >= 2 && fmod(step, 10000) < 0.5 ) {
      fprintf(stderr, "step %g, nc %d -> %d, snc %g -> %g, <nc> %g -> %g, "
          "dvdx %g(%g*%g), kT %g\n",
          step, nc0, nc1, snc0, snc1, nc0av, nc1av,
          dvdx, dvdx1, dvdx2, kT);
      if ( verbose >= 3 ) getchar();
    }
  } else {
    /* pop the old state, reverse the velocity */
    hmc_pop(hmc, go->x, go->v, go->f, idat, fdat, 1);
    nc1 = idat[0];
    snc1 = fdat[0];
    go->epot = fdat[1];
  }

  for ( i = 0; i < n; i++ ) { /* VV part 2 */
    vsinc(go->v[i], go->f[i], dth / go->m[i]);
  }

  // * vscramble is needed, see ANCHOR 1 */
  go->ekin = cago_vscramble(go, go->v, nvswaps);

  /* update the Wang-Landau object */
  wl_addi(wl, nc1);
  wl_updatelnf(wl);

  return acc;
}



#endif /* CAGONC_H__ */

