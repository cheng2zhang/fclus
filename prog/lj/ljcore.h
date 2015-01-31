#ifndef LJCORE_H__
#define LJCORE_H__



/* this file collect routines common to all dimensions
 * define the dimension D before including this file
 * Note: coordinates are not reduced */



#include "mtrand.h"
#include "util.h"
#include "graph.h"
#include "hmc.h"
#include "wl.h"



typedef struct {
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  double rho;
  double l, vol;
  double rc2, rc;
  double rcdef; /* preferred cutoff */

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double (*x2)[D];
  double *r2ij; /* pair distances */
  double *r2i; /* pair distances from i for MC */
  double epot, ep0, eps;
  double vir;
  double ekin;
  double epot_shift;
  double epot_tail;
  double p_tail;

  graph_t *g, *g2;
  double rcls; /* cluster cutoff */
  double lamcls; /* coupling factor for the clustering energy */
  const double *vcls; /* bias potential */
  int cseed; /* seed of cluster */

  double ediv; /* division energy */
  double lamdiv; /* coupling factor the division energy */
  int *idiv; /* division id, either 0 or 1 */
  int *idiv2; /* trial division */

  /* spherical division */
  int divseed; /* center of the division sphere */
  double divr; /* radius of the division sphere */
  double divrmin;
  double divrmax;
} lj_t;



/* the following functions are dimension D dependent */
static void lj_initfcc(lj_t *lj);
static double lj_gettail(lj_t *lj, double rho, int n, double *ptail);
static void lj_shiftang(double (*x)[D], double (*v)[D], int n);



/* set density and compute tail corrections */
static void lj_setrho(lj_t *lj, double rho)
{
  double irc;

  lj->rho = rho;
  lj->vol = lj->n/rho;
  lj->l = pow(lj->vol, 1./D);
  if ((lj->rc = lj->rcdef) > lj->l/2) lj->rc = lj->l/2;
  lj->rc2 = lj->rc * lj->rc;
  irc = 1/lj->rc;
  irc *= irc * irc;
  irc *= irc;
  lj->epot_shift = 4*irc*(irc - 1);
  lj->epot_tail = lj_gettail(lj, rho, lj->n, &lj->p_tail);
}



/* remove the center of mass motion */
static void lj_rmcom(double (*x)[D], int n)
{
  int i;
  double rc[D] = {0};

  for ( i = 0; i < n; i++ )
    vinc(rc, x[i]);
  vsmul(rc, 1./n);
  for ( i = 0; i < n; i++ )
    vdec(x[i], rc);
}



/* open an LJ system */
static lj_t *lj_open(int n, double rho, double rcdef,
    double rcls, const double *vcls)
{
  lj_t *lj;
  int i, d;

  xnew(lj, 1);
  lj->n = n;
  lj->dof = n * D - D * (D+1)/2;
  lj->rcdef = rcdef;

  xnew(lj->x, n);
  xnew(lj->v, n);
  xnew(lj->f, n);
  xnew(lj->x2, n);
  xnew(lj->r2ij, n * n);
  xnew(lj->r2i, n);

  lj_setrho(lj, rho);

  lj_initfcc(lj);

  /* initialize random velocities */
  for (i = 0; i < n; i++)
    for ( d = 0; d < D; d++ )
      lj->v[i][d] = randgaus();

  lj_rmcom(lj->v, lj->n);
  lj_shiftang(lj->x, lj->v, lj->n);

  lj->g = graph_open(n);
  lj->g2 = graph_open(n);
  lj->rcls = rcls;
  lj->vcls = vcls;

  lj->cseed = 0; /* seed particle of the cluster */

  lj->lamcls = 1;
  lj->lamdiv = 0;
  xnew(lj->idiv, n);
  xnew(lj->idiv2, n);
  /* initially put everything within a single division */
  for ( i = 0; i < n; i++ ) {
    lj->idiv[i] = 0;
    lj->idiv2[i] = 0;
  }
  lj->ediv = 0;

  /* division sphere */
  lj->divseed = 0;
  lj->divrmin = 0.9;
  lj->divrmax = lj->l * sqrt(0.75) + 0.00001;
  lj->divr = 2; //lj->divrmax;
  return lj;
}



/* close the lj object */
static void lj_close(lj_t *lj)
{
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj->x2);
  free(lj->r2ij);
  free(lj->r2i);
  graph_close(lj->g);
  graph_close(lj->g2);
  free(lj->idiv);
  free(lj->idiv2);
  free(lj);
}



#define LJ_PBC(x, l, invl) { (x) -= ((int)((x)*invl + 1000.5) - 1000.)*l; }



static double *lj_vpbc(double *v, double l, double invl)
{
  int d;
  for ( d = 0; d < D; d++ )
    LJ_PBC(v[d], l, invl);
  return v;
}



static double lj_pbcdist2(double *dx, const double *a, const double *b,
    double l, double invl)
{
  lj_vpbc(vdiff(dx, a, b), l, invl);
  return vsqr( dx );
}



#define lj_mkgraph(lj, g) lj_mkgraph_low(lj, g, lj->rcls)

/* build a graph and do clustering */
static void lj_mkgraph_low(lj_t *lj, graph_t *g, double rm)
{
  int i, j, n = lj->n;
  double rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ )
      if ( lj->r2ij[i*n + j] < rm2 )
        graph_link(g, i, j);
  graph_clus(g);
}



#define lj_mkgraph2(lj, g, k) lj_mkgraph2_low(lj, g, k, lj->rcls)

/* build a graph with the distances from k computed r2i */
static void lj_mkgraph2_low(lj_t *lj, graph_t *g, int k, double rm)
{
  int i, j, n = lj->n;
  double r2, rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      if ( i == k ) {
        r2 = lj->r2i[j];
      } else if ( j == k ) {
        r2 = lj->r2i[i];
      } else {
        r2 = lj->r2ij[i*n + j];
      }
      if ( r2 < rm2 )
        graph_link(g, i, j);
    }
  }
  graph_clus(g);
}



#define lj_energy(lj) \
  lj->epot = lj_energy_low(lj, lj->x, lj->r2ij, \
      &lj->vir, &lj->ep0, &lj->eps)

/* compute force and virial, return energy */
__inline static double lj_energy_low(lj_t *lj, double (*x)[D],
    double *r2ij, double *virial, double *ep0, double *eps)
{
  double dx[D], dr2, ir6, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) continue;
      dr2 = 1 / dr2;
      ir6 = dr2 * dr2 * dr2;
      vir += ir6 * (48 * ir6 - 24); /* f.r */
      ep += 4 * ir6 * (ir6 - 1);
      npr++;
    }
  }
  if (virial) *virial = vir;
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  return ep + lj->epot_tail; /* unshifted energy */
}



#define lj_force(lj) \
  lj->epot = lj_force_low(lj, lj->x, lj->f, lj->r2ij, &lj->vir, &lj->ep0, &lj->eps)

/* compute force and virial, return energy
 * the pair distances are recomputed */
__inline static double lj_force_low(lj_t *lj, double (*x)[D], double (*f)[D],
    double *r2ij, double *virial, double *ep0, double *eps)
{
  double dx[D], fi[D], dr2, ir6, fs, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) continue;
      dr2 = 1 / dr2;
      ir6 = dr2 * dr2 * dr2;
      fs = ir6 * (48 * ir6 - 24); /* f.r */
      vir += fs; /* f.r */
      fs *= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += 4 * ir6 * (ir6 - 1);
      npr++;
    }
    vinc(f[i], fi);
  }
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - npr * lj->epot_shift; /* shifted energy */
  if (virial) *virial = vir;
  return ep + lj->epot_tail; /* unshifted energy */
}



/* compute pressure */
static double lj_calcp(lj_t *lj, double tp)
{
  return (lj->dof * tp + lj->vir) / (D * lj->vol) + lj->p_tail;
}



/* velocity-verlet */
__inline static void lj_vv(lj_t *lj, double dt)
{
  int i, n = lj->n;
  double dth = dt * .5, l = lj->l;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(lj->v[i], lj->f[i], dth);
    vwrap( vsinc(lj->x[i], lj->v[i], dt), l );
  }
  lj_force(lj);
  for (i = 0; i < n; i++) /* VV part 2 */
    vsinc(lj->v[i], lj->f[i], dth);
}



/* compute the kinetic energy */
static double lj_ekin(double (*v)[D], int n)
{
  int i;
  double ek = 0;
  for ( i = 0; i < n; i++ ) ek += vsqr( v[i] );
  return ek/2;
}



/* randomly swap the velocities of m pairs of particles */
__inline static double lj_vscramble(double (*v)[D], int n, int m)
{
  int im, i, j;
  double vt[D];

  for ( im = 0; im < m; im++ ) {
    i = (int) (rand01() * n);
    j = (i + 1 + (int) (rand01() * (n - 1))) % n;
    vcopy(vt, v[i]);
    vcopy(v[i], v[j]);
    vcopy(v[j], vt);
  }
  return lj_ekin(v, n);
}



#define lj_vrescale(lj, tp, dt) \
  lj_vrescale_low(lj->v, lj->n, lj->dof, tp, dt)

/* exact velocity rescaling thermostat */
__inline static double lj_vrescale_low(double (*v)[D], int n,
    int dof, double tp, double dt)
{
  int i;
  double ek1, ek2, s, c, r, r2;

  c = (dt < 700) ? exp(-dt) : 0;
  ek1 = lj_ekin(v, n);
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp / 2);
  if (ek2 < 0) ek2 = 0;
  s = sqrt(ek2/ek1);
  for (i = 0; i < n; i++) vsmul(v[i], s);
  return ek2;
}



/* position Langevin barostat, with coordinates only
 * NOTE: the first parameter is the degree of freedom
 * the scaling is r = r*s
 * set cutoff to half of the box */
__inline static void lj_langp0(lj_t *lj, double dt,
    double tp, double pext, int ensx)
{
  double pint, amp, s, dlnv;
  int i;

  pint = lj_calcp(lj, tp);
  amp = sqrt(2 * dt);
  dlnv = ((pint - pext) * lj->vol / tp + 1 - ensx) * dt + amp * randgaus();
  s = exp( dlnv / D );
  lj->vol *= exp( dlnv );
  lj_setrho(lj, lj->n / lj->vol);
  for ( i = 0; i < lj->n; i++ )
    vsmul(lj->x[i], s);
  lj_force(lj);
}



/* displace a random particle i, return i */
static int lj_randmv(lj_t *lj, double *xi, double amp)
{
  int i, d;

  i = (int) (rand01() * lj->n);
  for ( d = 0; d < D; d++ )
    xi[d] = lj->x[i][d] + (rand01() * 2 - 1) * amp;
  return i;
}



/* compute pair energy */
__inline static int lj_pair(double dr2, double rc2, double *u, double *vir)
{
  if ( dr2 < rc2 ) {
    double invr2 = 1 / dr2;
    double invr6 = invr2 * invr2 * invr2;
    *vir = invr6 * (48 * invr6 - 24); /* f.r */
    *u  = 4.f * invr6 * (invr6 - 1);
    return 1;
  } else {
    *vir = 0;
    *u = 0;
    return 0;
  }
}



/* return the energy change from displacing x[i] to xi */
__inline static double lj_depot(lj_t *lj, int i, double *xi, double *vir)
{
  int j, n = lj->n;
  double l = lj->l, invl = 1/l, rc2 = lj->rc2, u, du, dvir;
  double dx[D], r2;

  u = 0.0;
  *vir = 0.0;
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) continue;
    r2 = lj->r2ij[i*n + j];
    if ( lj_pair(r2, rc2, &du, &dvir) ) {
      u -= du;
      *vir -= dvir;
    }
    r2 = lj_pbcdist2(dx, xi, lj->x[j], l, invl);
    if ( lj_pair(r2, rc2, &du, &dvir) ) {
      u += du;
      *vir += dvir;
    }
    lj->r2i[j] = r2;
  }
  lj->r2i[i] = 0;
  return u;
}



/* commit a particle displacement */
__inline static void lj_commit(lj_t *lj, int i, const double *xi,
    double du, double dvir, double dudiv)
{
  int j, n = lj->n;

  vwrap( vcopy(lj->x[i], xi), lj->l );
  lj->ep0 += du;
  lj->epot += du;
  lj->vir += dvir;
  lj->ediv += dudiv;
  for ( j = 0; j < n; j++ ) {
    lj->r2ij[i*n + j] = lj->r2i[j];
    lj->r2ij[j*n + i] = lj->r2i[j];
  }
  for ( j = 0; j < n; j++ ) {
    lj->idiv2[j] = lj->idiv[j];
  }
  graph_copy(lj->g, lj->g2);

#if 0 /* debug code */
  printf("<<< ep %g, ", lj->epot);
  graph_clus_print(lj->g);
  lj_mkgraph(lj, lj->g);
  printf(">>> ep %g, ", lj->epot);
  graph_clus_print(lj->g);
  getchar();
#endif
}



/* compute the cluster energy
 * call lj_mkgraph() first */
__inline static double lj_eclus(const lj_t *lj, const graph_t *g)
{
  if ( lj->vcls ) return 0.0;
  else return lj->vcls[ graph_getcsize(g, lj->cseed) ];
}



__inline static int lj_getdivsize(lj_t *lj, const int *idiv)
{
  int i, c = 0, n = lj->n;

  for ( i = 0; i < n; i++ ) {
    c += (idiv[i] == 0);
  }
  return c;
}



/* compute the division */
__inline static int lj_getdiv(lj_t *lj, int *idiv,
    int seed, double divr)
{
  int i, n = lj->n, cnt = 0;
  double divr2 = divr * divr;

  for ( i = 0; i < n; i++ ) {
    idiv[i] = ( lj->r2ij[seed*n + i] > divr2 );
    cnt += (idiv[i] == 0);
  }
  return cnt;
}



/* compute the division */
__inline static void lj_getdiv2(lj_t *lj, int *idiv,
    int seed, double divr, int k)
{
  int i, n = lj->n;
  double r2, divr2 = divr * divr;

  for ( i = 0; i < n; i++ ) {
    if ( i == seed ) {
      idiv[i] = 0;
      continue;
    }
    if ( i == k ) {
      r2 = lj->r2i[seed];
    } else if ( seed == k ) {
      r2 = lj->r2i[i];
    } else {
      r2 = lj->r2ij[seed*n + i];
    }
    idiv[i] = ( r2 > divr2 );
  }
}



#define LJ_R2WCA 1.2599210498948732

/* compute the division energy
 * this functions uses lj->r2ij, so call lj_energy()/lj_force() first
 * idiv[i] gives the division id of particle i */
__inline static double lj_ediv(lj_t *lj, int *idiv,
    int divseed, double divr)
{
  int i, j, n = lj->n;
  double dr2, ir6, rc2 = lj->rc2, ep = 0;

  /* compute the division */
  lj_getdiv(lj, idiv, divseed, divr);

  for ( i = 0; i < n; i++ ) {
    /* i must be in division 0 */
    if ( idiv[i] != 0 ) continue;
    for ( j = 0; j < n; j++ ) {
      /* j must not be in division 1 */
      if ( idiv[j] == 0 ) continue;

      dr2 = lj->r2ij[i*n + j];
      if ( dr2 < LJ_R2WCA ) {
        ep += 1;
      } else if ( dr2 < rc2 ) {
        ir6 = 1/(dr2 * dr2 * dr2);
        ep -= 4 * ir6 * (ir6 - 1);
      }
    }
  }
  return ep;
}



/* compute the trial division energy
 * this functions uses lj->r2ij, lj->r2i */
__inline static double lj_ediv2(lj_t *lj, int *idiv,
    int divseed, double divr, int k)
{
  int i, j, n = lj->n;
  double dr2, ir6, rc2 = lj->rc2, ep = 0;

  /* compute the division */
  lj_getdiv2(lj, idiv, divseed, divr, k);
  for ( i = 0; i < n; i++ ) {
    /* i must be in division 0 */
    if ( idiv[i] != 0 ) continue;
    for ( j = 0; j < n; j++ ) {
      /* j must not be in division 0 */
      if ( idiv[j] == 0 ) continue;

      if ( i == k ) {
        dr2 = lj->r2i[j];
      } else if ( j == k ) {
        dr2 = lj->r2i[i];
      } else {
        dr2 = lj->r2ij[i*n + j];
      }

      if ( dr2 < LJ_R2WCA ) {
        ep += 1;
      } else if ( dr2 < rc2 ) {
        ir6 = 1/(dr2 * dr2 * dr2);
        ep -= 4 * ir6 * (ir6 - 1);
      }
    }
  }
  return ep;
}



/* Metropolis algorithm */
__inline static int lj_metro(lj_t *lj, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0, dutot, dvir = 0;
  double ucls, ducls = 0;
  double udiv, dudiv = 0;

  i = lj_randmv(lj, xi, amp);
  du = lj_depot(lj, i, xi, &dvir);

  lj_mkgraph2(lj, lj->g2, i);

  /* compute the clustering energy */
  if ( fabs( lj->lamcls ) > 0 ) {
    ucls = lj_eclus(lj, lj->g2);
    /* since the clustering potential might have been changed
     * we have to recompute the old clustering energy */
    ducls = ucls - lj_eclus(lj, lj->g);
  }

  /* compute the division energy */
  if ( fabs( lj->lamdiv ) > 0 ) {
    udiv = lj_ediv2(lj, lj->idiv2, lj->divseed, lj->divr, i);
    //lj->ediv = lj_ediv(lj, lj->idiv, lj->divseed, lj->divr);
    //printf("ediv %g, %g; seed %d, r %g\n", lj->ediv, lj_ediv(lj, lj->idiv, lj->divseed, lj->divr), lj->divseed, lj->divr); getchar();
    dudiv = udiv - lj->ediv;
  }

  dutot = bet * (du + lj->lamdiv * dudiv) + lj->lamcls * ducls;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    lj_commit(lj, i, xi, du, dvir, dudiv);
    return 1;
  }
  return 0;
}



/* attempt to change the seed for clustering */
__inline static int lj_changeseed(lj_t *lj, const graph_t *g)
{
  int i, sz0, sz1, n = lj->n, acc = 0;

  i = (lj->cseed + 1 + (int) (rand01() * (n - 1))) % n;
  sz0 = graph_getcsize(g, lj->cseed);
  sz1 = graph_getcsize(g, i);
  if ( sz0 == sz1 || lj->vcls == NULL ) {
    acc = 1;
  } else {
    double dv = lj->vcls[ sz1 ] - lj->vcls[ sz0 ];
    if ( dv < 0 ) {
      acc = 1;
    } else {
      double r = rand01();
      acc = (r < exp( -lj->lamcls * dv ));
    }
  }
  if ( acc ) {
    lj->cseed = i;
  }
  return acc;
}



/* attempt to change the division radius */
__inline static int lj_changediv(lj_t *lj, double amp, double bet)
{
  int i, acc = 0, n = lj->n;
  double rnew, ediv, dv;

  /* try to move i to the opposite division */
  rnew = lj->divr + amp * (2 * rand01() - 1);
  if ( rnew > lj->divrmax || rnew < lj->divrmin ) {
    return 0;
  }

  /* compute the energy change caused by the change */
  ediv = lj_ediv(lj, lj->idiv2, lj->divseed, rnew);
  //printf("ediv %g, %g; seed %d, r %g\n", lj->ediv, lj_ediv(lj, lj->idiv, lj->divseed, lj->divr), lj->divseed, lj->divr); getchar();
  dv = ediv - lj->ediv;
  if ( dv < 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = (r < exp( - bet * lj->lamdiv * dv ));
  }

  if ( acc ) {
    lj->divr = rnew;
    for ( i = 0; i < n; i++ ) {
      lj->idiv[i] = lj->idiv2[i];
    }
    lj->ediv = ediv;
  }
  return acc;
}



__inline static int lj_hmc(lj_t *lj, hmc_t *hmc, int *csize1)
{
  /* compute the current cluster size */
  int csize = graph_getcsize(lj->g, lj->cseed);
  int acc, idat[2] = { csize, lj->cseed };
  /* hmc->idat[0] is the previous size */
  double dv = 0;

  if ( lj->vcls ) {
    dv = lj->vcls[ csize ] - lj->vcls[ hmc->idat[0] ];
  }
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = ( r < exp( -dv ) );
  }
  if ( acc ) {
    hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);
  } else {
    hmc_pop(hmc, lj->x, lj->v, lj->f, idat, &lj->epot, 1);
    lj->cseed = idat[1];
    /* lj->r2ij should be refreshed in the next step */
    //lj_force(lj);
    //lj_mkgraph(lj, lj->g);
  }
  *csize1 = idat[0];
  return acc;
}



/* wrap coordinates such that particles stay in the box */
__inline static int lj_wrapbox(lj_t *lj,
    double (*xin)[D], double (*xout)[D])
{
  int i, n = lj->n;
  double l = lj->l;

  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }
  return 0;
}



/* wrap coordinates such that particles
 * in the same cluster stay close */
__inline static int lj_wrapclus(lj_t *lj,
    double (*xin)[D], double (*xout)[D], graph_t *g)
{
  int ic, i, j, n = lj->n, head, end;
  double l = lj->l, invl = 1 / l, dx[D];

  lj_mkgraph(lj, g);
  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }

  for ( ic = 0; ic < g->nc; ic++ ) {
    /* find the seed of this cluster */
    for ( i = 0; i < n; i++ )
      if ( g->cid[i] == ic )
        break;
    if ( i >= n ) {
      fprintf(stderr, "no particle belongs to cluster %d\n", ic);
      return -1;
    }

    g->queue[ head = 0 ] = i;
    g->cid[ i ] = -1;
    end = 1;
    for ( ; head < end; head++ ) {
      i = g->queue[head];
      /* add neighbors of i into the queue */
      for ( j = 0; j < n; j++ ) {
        if ( graph_linked(g, i, j) && g->cid[j] >= 0 ) {
          g->cid[j] = -1;
          g->queue[ end++ ] = j;
          /* align j with i */
          lj_vpbc(vdiff(dx, xout[j], xout[i]), l, invl);
          vinc( vcopy(xout[j], xout[i]), dx );
        }
      }
    }
    if ( end != g->csize[ic] ) {
      fprintf(stderr, "cluster %d: size %d vs %d\n", ic, end, g->csize[ic]);
    }
  }
  return 0;
}



/* write positions (and possibly velocities) */
__inline static int lj_writepos(lj_t *lj,
    double (*x)[D], double (*v)[D],
    const char *fn, int wrap)
{
  FILE *fp;
  int i, d, n = lj->n;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  /* wrap coordinates according to clusters */
  if ( wrap ) {
    lj_wrapclus(lj, x, lj->x2, lj->g2);
  } else {
    lj_wrapbox(lj, x, lj->x2);
  }

  fprintf(fp, "# %d %d %d %.14e\n", D, n, (v != NULL), lj->l);
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ )
      fprintf(fp, "%.14e ", lj->x2[i][d]);
    if ( v != NULL )
      for ( d = 0; d < D; d++ )
        fprintf(fp, "%.14e ", v[i][d]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



__inline static double lj_sphrvol(double r)
{
  if ( D == 3 ) return 4 * M_PI * r * r * r / 3;
  else return 4 * M_PI * r * r;
}



#endif /* LJCORE_H__ */
