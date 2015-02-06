#ifndef LJDIV_H__
#define LJDIV_H__



/* division based energy */



#include "ljcore.h"



typedef struct {
  lj_t *lj; /* reference lj struct */

  double ediv; /* division energy */
  double lamdiv; /* coupling factor the division energy */
  int *idiv; /* division id, either 0 or 1 */
  int *idiv2; /* trial division */

  /* spherical division */
  int divseed; /* center of the division sphere */
  double divr; /* radius of the division sphere */
  double divrmin;
  double divrmax;
} ljdiv_t;



static ljdiv_t *ljdiv_open(lj_t *lj)
{
  int i, n = lj->n;
  ljdiv_t *d;

  xnew(d, 1);
  d->lj = lj;
  d->lamdiv = 1;
  xnew(d->idiv, n);
  xnew(d->idiv2, n);
  /* initially put everything within a single division */
  for ( i = 0; i < n; i++ ) {
    d->idiv[i] = 0;
    d->idiv2[i] = 0;
  }
  d->ediv = 0;

  /* division sphere */
  d->divseed = 0;
  d->divrmin = 0.9;
  d->divrmax = lj->l * sqrt(0.75) + 0.00001;
  d->divr = 2; //lj->divrmax;

  return d;
}



static void ljdiv_close(ljdiv_t *d)
{
  free(d->idiv);
  free(d->idiv2);
  free(d);
}



static void ljdiv_commit(ljdiv_t *d, double dudiv)
{
  int i, n = d->lj->n;

  d->ediv += dudiv;
  for ( i = 0; i < n; i++ ) {
    d->idiv[i] = d->idiv2[i];
  }
}



__inline static int ljdiv_getdivsize(ljdiv_t *d, const int *idiv)
{
  int i, c = 0, n = d->lj->n;

  for ( i = 0; i < n; i++ ) {
    c += (idiv[i] == 0);
  }
  return c;
}



/* compute the division */
__inline static int ljdiv_getdiv(ljdiv_t *d, int *idiv,
    int seed, double divr)
{
  int i, n = d->lj->n, cnt = 0;
  double divr2 = divr * divr;

  for ( i = 0; i < n; i++ ) {
    idiv[i] = ( d->lj->r2ij[seed*n + i] > divr2 );
    cnt += (idiv[i] == 0);
  }
  return cnt;
}



/* compute the division */
__inline static void ljdiv_getdiv2(ljdiv_t *d, int *idiv,
    int seed, double divr, int k)
{
  int i, n = d->lj->n;
  double r2, divr2 = divr * divr;

  for ( i = 0; i < n; i++ ) {
    if ( i == seed ) {
      idiv[i] = 0;
      continue;
    }
    if ( i == k ) {
      r2 = d->lj->r2i[seed];
    } else if ( seed == k ) {
      r2 = d->lj->r2i[i];
    } else {
      r2 = d->lj->r2ij[seed*n + i];
    }
    idiv[i] = ( r2 > divr2 );
  }
}



#define LJ_R2WCA 1.2599210498948732

/* compute the division energy
 * this functions uses lj->r2ij, so call lj_energy()/lj_force() first
 * idiv[i] gives the division id of particle i */
__inline static double ljdiv_ediv(ljdiv_t *d, int *idiv,
    int divseed, double divr)
{
  int i, j, n = d->lj->n;
  double dr2, ir6, rc2 = d->lj->rc2, ep = 0;

  /* compute the division */
  ljdiv_getdiv(d, idiv, divseed, divr);

  for ( i = 0; i < n; i++ ) {
    /* i must be in division 0 */
    if ( idiv[i] != 0 ) continue;
    for ( j = 0; j < n; j++ ) {
      /* j must not be in division 1 */
      if ( idiv[j] == 0 ) continue;

      dr2 = d->lj->r2ij[i*n + j];
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
__inline static double ljdiv_ediv2(ljdiv_t *d, int *idiv,
    int divseed, double divr, int k)
{
  int i, j, n = d->lj->n;
  double dr2, ir6, rc2 = d->lj->rc2, ep = 0;

  /* compute the division */
  ljdiv_getdiv2(d, idiv, divseed, divr, k);
  for ( i = 0; i < n; i++ ) {
    /* i must be in division 0 */
    if ( idiv[i] != 0 ) continue;
    for ( j = 0; j < n; j++ ) {
      /* j must not be in division 0 */
      if ( idiv[j] == 0 ) continue;

      if ( i == k ) {
        dr2 = d->lj->r2i[j];
      } else if ( j == k ) {
        dr2 = d->lj->r2i[i];
      } else {
        dr2 = d->lj->r2ij[i*n + j];
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
__inline static int ljdiv_metro(ljdiv_t *d, double amp, double bet)
{
  lj_t *lj = d->lj;
  int i, acc = 0;
  double xi[D], r, du = 0, dutot, dvir = 0;
  double udiv, dudiv = 0;

  i = lj_randmv(lj, xi, amp);
  du = lj_depot(lj, i, xi, &dvir);

  /* compute the division energy */
  udiv = ljdiv_ediv2(d, d->idiv2, d->divseed, d->divr, i);
  //lj->ediv = ljdiv_ediv(d, d->idiv, d->divseed, d->divr);
  //printf("ediv %g, %g; seed %d, r %g\n", d->ediv, ljdiv_ediv(d, d->idiv, d->divseed, d->divr), d->divseed, d->divr); getchar();
  dudiv = udiv - d->ediv;

  dutot = bet * (du + d->lamdiv * dudiv);
  if ( dutot < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -dutot ) );
  }
  if ( acc ) {
    lj_commit(lj, i, xi, du, dvir);
    ljdiv_commit(d, dudiv);
    return 1;
  }
  return 0;
}



/* attempt to change the division radius */
__inline static int ljdiv_changediv(ljdiv_t *d, double amp, double bet)
{
  lj_t *lj = d->lj;
  int i, acc = 0, n = lj->n;
  double rnew, ediv, dv;

  /* try to move i to the opposite division */
  rnew = d->divr + amp * (2 * rand01() - 1);
  if ( rnew > d->divrmax || rnew < d->divrmin ) {
    return 0;
  }

  /* compute the energy change caused by the change */
  ediv = ljdiv_ediv(d, d->idiv2, d->divseed, rnew);
  //printf("ediv %g, %g; seed %d, r %g\n", d->ediv, ljdiv_ediv(d, d->idiv, d->divseed, d->divr), d->divseed, d->divr); getchar();
  dv = ediv - d->ediv;
  if ( dv < 0 ) {
    acc = 1;
  } else {
    double r = rand01();
    acc = (r < exp( - bet * d->lamdiv * dv ));
  }

  if ( acc ) {
    d->divr = rnew;
    for ( i = 0; i < n; i++ ) {
      d->idiv[i] = d->idiv2[i];
    }
    d->ediv = ediv;
  }
  return acc;
}



__inline static double lj_sphrvol(double r)
{
  if ( D == 3 ) return 4 * M_PI * r * r * r / 3;
  else return 4 * M_PI * r * r;
}



#endif /* LJDIV_H__ */
