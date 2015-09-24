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



/* structure to map from the bias potential to
 * the number of contacts */
typedef struct {
  double xmin, xmax, dx;
  int n;
  double (*arr)[2];
} xymap_t;



__inline static xymap_t *xymap_open(double xmin, double xmax, double dx)
{
  xymap_t *xy;

  xnew(xy, 1);
  xy->xmin = xmin;
  xy->dx = dx;
  xy->n = (int) ((xmax - xmin) / dx);
  xy->xmax = xmin + dx * xy->n;
  xnew(xy->arr, xy->n);
  return xy;
}



__inline static void xymap_close(xymap_t *xy)
{
  free(xy->arr);
  free(xy);
}



__inline static int xymap_getid(xymap_t *xy, double x)
{
  int i;

  if ( x < xy->xmin ) return -1;
  i = (int) ((x - xy->xmin) / xy->dx);
  return ( i < xy->n ) ? i : -1;
}


__inline static int xymap_add(xymap_t *xy, double x, double y)
{
  int i = xymap_getid(xy, x);

  if ( i < 0 ) return -1;
  xy->arr[i][0] += 1;
  xy->arr[i][1] += y;
  return 0;
}



/* compute the y value */
__inline static double xymap_gety(xymap_t *xy, double x)
{
  int i = xymap_getid(xy, x);

  if ( i < 0 || xy->arr[i][0] <= 0) return 0;
  return xy->arr[i][1] / xy->arr[i][0];
}



/* compute the derivative */
__inline static double xymap_getdydx(xymap_t *xy, double x)
{
  int i, i1, i2;
  double dx, y1, y2;

  /* first compute i1 and i2 */
  i = xymap_getid(xy, x);
  if ( i < 0 ) return 0;
  dx = x - (xy->xmin + i * xy->dx);
  if ( dx < 0.5 * xy->dx ) { /* left half */
    if ( i > 0 ) {
      i1 = 0;
      i2 = 1;
    } else {
      i1 = i - 1;
      i2 = i;
    }
  } else { /* right half */
    if ( i < xy->n - 1 ) {
      i1 = i;
      i2 = i + 1;
    } else {
      i1 = xy->n - 2;
      i2 = xy->n - 1;
    }
  }

  /* compute the derivative */
  if ( xy->arr[i1][0] > 0 && xy->arr[i2][0] > 0 ) {
    y1 = xy->arr[i1][1] / xy->arr[i1][0];
    y2 = xy->arr[i2][1] / xy->arr[i2][0];
    return ( y2 - y1 ) / xy->dx;
  }

  return 0;
}



/* attractive part of the WCA de-composition */
__inline static double potattract(vct a, vct b,
    double rc2, double eps, vct fa, vct fb)
{
  double dx[D], dr2, invr2, invr6, amp;

  vdiff(dx, a, b);
  dr2 = vsqr( dx );
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;

  if ( invr6 > 0.5 ) {
    return -eps;
  } else {
    if ( fa ) {
      amp = eps * ( 48 * invr6 - 24 ) * invr6 / dr2;
      vsinc(fa, dx,  amp);
      vsinc(fb, dx, -amp);
    }
    return 4 * eps * (invr6 - 1) * invr6;
  }
}



/* bias potential */
__inline static double cago_vnc(cago_t *go, vct *x)
{
  int i, j, id, n = go->n;
  double ene;

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 4; j < n; j++ ) {
      id = i * n + j;
      if ( !go->iscont[id] ) continue;
      ene += potattract(x[i], x[j], go->r2ref[id],
          1.0, NULL, NULL);
    }
  }
  return ene;
}



/* bias potential with force */
__inline static double cago_forcenc(cago_t *go, vct *x, vct *f)
{
  int i, j, id, n = go->n;
  double ene;

  for ( i = 0; i < n; i++ ) {
    vzero(f[i]);
  }

  for ( i = 0; i < n; i++ ) {
    for ( j = i + 4; j < n; j++ ) {
      id = i * n + j;
      if ( !go->iscont[id] ) continue;
      ene += potattract(x[i], x[j], go->r2ref[id],
          1.0, f[i], f[j]);
    }
  }
  return ene;
}



/* initialize an HMC object for implicit NC */
__inline static hmc_t *cago_ihmc_nc_init(cago_t *go, int *idat, double **pfdat)
{
  hmc_t *hmc;
  double *fdat;
  int cnt = go->n * D + 2;

  /* make a hybrid Monte-Carlo object
   * with 1 extra integer
   * and n*D + 1 extra floating-point numbers:
   * NC, potential energy, and the fit structure */
  hmc = hmc_open(go->n, 1, cnt);

  xnew(fdat, cnt);
  idat[0] = cago_ncontacts(go, go->x, -1, NULL, NULL);
  fdat[0] = cago_vnc(go, go->x);
  fdat[1] = go->epot;
  memcpy(fdat + 2, go->x1, go->n * D * sizeof(double));
  *pfdat = fdat;

  /* push the initial state */
  hmc_push(hmc, go->x, go->v, go->f, idat, fdat);

  return hmc;
}



/* velocity Verlet with RMSD bias (implicit hybrid MC) */
__inline static int cago_vv_nc(cago_t *go, double fs, double dt,
    double (*xf)[D], wl_t *wl, xymap_t *xy,
    double mflmin, double mflmax, double mfhmin, double mfhmax,
    double kT, hmc_t *hmc, int *idat, double *fdat)
{
  int i, n = go->n, acc;
  double dth = 0.5 * dt * fs;
  int nc1, nc0;
  double vnc1, vnc0;
  double nc1av, nc0av, dv0, dv1, dv, dvdx;

  for ( i = 0; i < n; i++ ) { /* VV part 1 */
    vsinc(go->v[i], go->f[i], dth / go->m[i] );
    vsinc(go->x[i], go->v[i], dt);
  }

  /* compute the normal force */
  go->epot = cago_force(go, go->x, go->f);

  /* compute the bias energy and force,
   * the force is saved in `go->x1` */
  vnc1 = cago_forcenc(go, go->x, go->x1);
  nc1av = xymap_gety(xy, vnc1);

  vnc0 = hmc->fdat[0];
  nc0av = xymap_gety(xy, vnc0);

  nc1 = cago_ncontacts(go, go->x, -1, NULL, NULL);
  nc0 = hmc->idat[0];

  /* 3. compute energy caused by changing
   *    the reference structure */
  /* get the potential of the new (current) state */
  dv1 = wl_getvi(wl, nc1);
  dv1 -= wl_getvi(wl, (int) (nc1av + 0.5) );
  /* get the potential of the old (stock) state */
  dv0 = wl_getvi(wl, nc0);
  dv0 -= wl_getvi(wl, (int) (nc0av + 0.5) );
  dv = dv1 - dv0;

  /* 4. decide if the change of reference is acceptable */
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dv ) );
  }

  if ( acc ) {
    /* accept the state */
    idat[0] = nc1;
    fdat[0] = vnc1;
    fdat[1] = go->epot;
    memcpy(fdat + 2, xf, n * D * sizeof(double));

    /* compute the force scaling factor */
    dvdx = kT * wl_getdvdx_v(wl, nc1av,
        mflmin, mflmax, mfhmin, mfhmax);
    dvdx *= xymap_getdydx(xy, vnc1);
    /* TODO: apply the bias force */
    for ( i = 0; i < n; i++ ) {
      vsinc(go->f[i], go->x1[i], -dvdx);
    }

    /* push the current state */
    hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);
  } else {
    /* pop the old state, reverse the velocity */
    hmc_pop(hmc, go->x, go->v, go->f, NULL, fdat, 1);
    nc1 = idat[0];
    vnc1 = fdat[0];
    go->epot = fdat[1];
    memcpy(xf, fdat + 2, n * D * sizeof(double));
  }

  for ( i = 0; i < n; i++ ) { /* VV part 2 */
    vsinc(go->v[i], go->f[i], dth / go->m[i]);
  }

  ///* if out-of-range RMSD is not allowed,
  // * vscramble is needed, see ANCHOR 1 */
  //go->ekin = cago_vscramble(go, go->v, nvswaps);

  wl_addi(wl, nc1);
  wl_updatelnf(wl);

  return acc;
}



#endif /* CAGONC_H__ */

