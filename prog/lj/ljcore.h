#ifndef LJCORE_H__
#define LJCORE_H__



/* define the dimension D before including this file
 * Note: coordinates are not reduced */



#include "mtrand.h"
#include "util.h" /* vct.h and mat.h are included within */
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
} lj_t;



/* functions that are dimension D dependent */
#if D == 2
#include "lj2d.h"
#else
#include "lj3d.h"
#endif



/* set density and compute tail corrections */
static void lj_setrho(lj_t *lj, double rho)
{
  double irc;

  lj->rho = rho;
  lj->vol = lj->n/rho;
  lj->l = pow(lj->vol, 1.0 / D);
  lj->rc = dblmin( lj->rcdef, lj->l * 0.5 );
  lj->rc2 = lj->rc * lj->rc;
  irc = 1 / lj->rc;
  irc *= irc * irc;
  irc *= irc;
  lj->epot_shift = 4 * irc * (irc - 1);
  lj->epot_tail = lj_gettail(lj, rho, lj->n, &lj->p_tail);
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
  for ( i = 0; i < n; i++ )
    for ( d = 0; d < D; d++ )
      lj->v[i][d] = randgaus();

  rmcom(lj->v, NULL, n);
  shiftang(lj->x, lj->v, NULL, n);

  lj->g = graph_open(n);
  lj->g2 = graph_open(n);
  lj->rcls = rcls;
  lj->vcls = vcls;

  lj->cseed = 0; /* seed particle of the cluster */
  lj->lamcls = 1;

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
    double du, double dvir)
{
  int j, n = lj->n;

  vwrap( vcopy(lj->x[i], xi), lj->l );
  lj->ep0 += du;
  lj->epot += du;
  lj->vir += dvir;
  for ( j = 0; j < n; j++ ) {
    lj->r2ij[i*n + j] = lj->r2i[j];
    lj->r2ij[j*n + i] = lj->r2i[j];
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
  return lj->vcls ? lj->vcls[ graph_getcsize(g, lj->cseed) ] : 0.0;
}



/* Metropolis algorithm */
__inline static int lj_metro(lj_t *lj, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0, dutot, dvir = 0;
  double ucls, ducls = 0;

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

  dutot = bet * du + lj->lamcls * ducls;
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    lj_commit(lj, i, xi, du, dvir);
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
