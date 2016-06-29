#ifndef LJMIXCORE_H__
#define LJMIXCORE_H__



/* define the dimension D before including this file
 * Note: coordinates are not reduced */



#include "mtrand.h"
#include "util.h"
#include "mdutil.h"



//#define NSMAX 2 /* maximal number of species */



typedef struct {
  int n; /* total number of particles */
  int dof; /* degrees of freedom */
  int ns; /* number of species */
  int np[NSMAX]; /* number of particles */
  double sig[NSMAX]; /* diameter */
  double eps[NSMAX]; /* energy unit */
  double sigij[NSMAX][NSMAX];
  double epsij[NSMAX][NSMAX];
  double rho[NSMAX]; /* density, to be computed */
  int *type;
  double l, vol;
  double rc2, rc;
  double rcdef; /* preferred cutoff */

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double *r2ij; /* pair distances */
  double *r2i; /* pair distances from i */
  double epot, ep0;
  double vir;
  double ekin;
  double epot_tail;
  double p_tail;
} ljmix_t;



/* functions that are dimension D dependent */
#if D == 2



/* initialize a fcc lattice */
static void ljmix_initfcc(ljmix_t *lj)
{
  int i, j, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      if ( (i + j) % 2 != 0 ) continue;
      /* add some noise to prevent two atoms happened to
       * be separated by precisely some special cutoff distance,
       * which might be half of the box */
      lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
      lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
      id++;
    }
  }
}



/* get the tail correction */
static double ljmix_gettail(ljmix_t *lj, double *pptail)
{
  int is, js, ns = lj->ns;
  double rd, irc, irc3, irc6, sig, utail, ptail;

  rd = lj->rc * lj->rc;
  utail = 0;
  ptail = 0;
  for ( is = 0; is < ns; is++ ) {
    for ( js = 0; js < ns; js++ ) {
      sig = (lj->sig[is] + lj->sig[js]) / 2;
      irc = sig / lj->rc;
      irc3 = irc * irc * irc;
      irc6 = irc3 * irc3;
      utail += M_PI*lj->rho[is]*lj->rho[js]*(0.4*irc6 - 1)*irc6*rd;
      ptail += M_PI*lj->rho[is]*lj->rho[js]*(2.4*irc6 - 3)*irc6*rd;
    }
  }
  if ( pptail != NULL ) {
    *pptail = ptail;
  }
  utail *= lj->vol;
  return utail;
}



#else /* D == 3 */



/* initialize a fcc lattice */
static void ljmix_initfcc(ljmix_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for ( id = 0, i = 0; i < n1 && id < n; i++ ) {
    for ( j = 0; j < n1 && id < n; j++ ) {
      for ( k = 0; k < n1 && id < n; k++ ) {
        if ( (i + j + k) % 2 != 0 ) continue;
        /* add some noise to prevent two atoms happened to
         * be separated by precisely some special cutoff distance,
         * which might be half of the box */
        lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][2] = (k + .5) * a + noise * (2*rand01() - 1);
        id++;
      }
    }
  }
}



/* get the tail correction */
static double ljmix_gettail(ljmix_t *lj, double *pptail)
{
  int is, js, ns = lj->ns;
  double rd, irc, irc3, irc6, sig, utail, ptail;

  rd = lj->rc * lj->rc * lj->rc;
  utail = 0;
  ptail = 0;
  for ( is = 0; is < ns; is++ ) {
    for ( js = 0; js < ns; js++ ) {
      sig = (lj->sig[is] + lj->sig[js]) / 2;
      irc = sig / lj->rc;
      irc3 = irc * irc * irc;
      irc6 = irc3 * irc3;
      utail +=  8*M_PI*lj->rho[is]*lj->rho[js]/9*(irc6 - 3.0)*irc6*rd;
      ptail += 32*M_PI*lj->rho[is]*lj->rho[js]/9*(irc6 - 1.5)*irc6*rd;
    }
  }
  if ( pptail != NULL ) {
    *pptail = ptail;
  }
  utail *= lj->vol;
  return utail;
}



#endif /* D == 3 */



/* set volume and compute tail corrections */
static void ljmix_setvol(ljmix_t *lj, double vol)
{
  int is;

  lj->vol = vol;
  for ( is = 0; is < lj->ns; is++ ) {
    lj->rho[is] = lj->np[is] / lj->vol;
  }
  lj->l = pow(lj->vol, 1.0 / D);
  lj->rc = dblmin( lj->rcdef, lj->l * 0.5 );
  lj->rc2 = lj->rc * lj->rc;
  lj->epot_tail = ljmix_gettail(lj, &lj->p_tail);
}



/* open an LJ system */
static ljmix_t *ljmix_open(int ns,
    int *np, double *sig, double *eps,
    double rho, double rcdef)
{
  ljmix_t *lj;
  int i, j, t, n, is, js, d;
  double vol;

  xnew(lj, 1);
  lj->ns = ns;

  /* compute the total number of particles */
  for ( n = 0, is = 0; is < ns; is++ ) {
    lj->np[is] = np[is];
    lj->sig[is] = sig[is];
    lj->eps[is] = eps[is];
    n += np[is];
  }
  lj->n = n;

  /* compute the densities */
  vol = n / rho;
  for ( is = 0; is < ns; is++ ) {
    lj->rho[is] = lj->np[is] / vol;
  }

  /* compute the pairwise LJ parameters */
  for ( is = 0; is < ns; is++ ) {
    for ( js = 0; js < ns; js++ ) {
      lj->sigij[is][js] = (lj->sig[is] + lj->sig[js]) / 2;
      lj->epsij[is][js] = sqrt( lj->eps[is] * lj->eps[js] );
    }
  }

  /* assign the types */
  xnew(lj->type, n);
  i = 0;
  for ( is = 0; is < ns; is++ ) {
    for ( j = 0; j < lj->np[is]; j++ ) {
      lj->type[ i++ ] = is;
    }
  }

  lj->dof = n * D - D; // * (D+1)/2;
  lj->rcdef = rcdef;

  xnew(lj->x, n);
  xnew(lj->v, n);
  xnew(lj->f, n);
  xnew(lj->r2ij, n * n);
  xnew(lj->r2i, n);

  ljmix_setvol(lj, vol);

  ljmix_initfcc(lj);

  /* randomly swap coordinates to mix things up */
  for ( t = 0; t < n * n; t++ ) {
    i = (int) (rand01() * n);
    j = (i + 1 + (int) (rand01() * (n - 1))) % n;
    vswap(lj->x[i], lj->x[j]);
  }

  /* initialize random velocities */
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ ) {
      lj->v[i][d] = randgaus();
    }
  }

  rmcom(lj->v, NULL, n);
  shiftang(lj->x, lj->v, NULL, n);

  return lj;
}



/* close the lj object */
static void ljmix_close(ljmix_t *lj)
{
  free(lj->type);
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj->r2ij);
  free(lj->r2i);
  free(lj);
}



#define LJMIX_PBC(x, l, invl) { (x) -= ((int)((x)*invl + 1000.5) - 1000.)*l; }



static double *ljmix_vpbc(double *v, double l, double invl)
{
  int d;
  for ( d = 0; d < D; d++ ) {
    LJMIX_PBC(v[d], l, invl);
  }
  return v;
}



static double ljmix_pbcdist2(double *dx, const double *a, const double *b,
    double l, double invl)
{
  ljmix_vpbc(vdiff(dx, a, b), l, invl);
  return vsqr( dx );
}



#define ljmix_energy(lj) \
  lj->epot = ljmix_energy_low(lj, lj->x, lj->r2ij, \
      &lj->vir, &lj->ep0)

/* compute force and virial, return energy */
__inline static double ljmix_energy_low(ljmix_t *lj,
    double (*x)[D],
    double *r2ij, double *virial, double *ep0)
{
  double dx[D], dr2, ir2, ir6, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l, sig, eps;
  int i, j, itp, jtp, npr = 0, n = lj->n;

  ep = vir = 0;
  for ( i = 0; i < n - 1; i++ ) {
    itp = lj->type[i];
    for ( j = i + 1; j < n; j++ ) {
      jtp = lj->type[j];
      dr2 = ljmix_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) {
        continue;
      }
      sig = lj->sigij[itp][jtp];
      eps = lj->epsij[itp][jtp];
      ir2 = (sig * sig) / dr2;
      ir6 = ir2 * ir2 * ir2;
      vir += eps * ir6 * (48 * ir6 - 24); /* f.r */
      ep += eps * 4 * ir6 * (ir6 - 1);
      npr++;
    }
  }
  if ( virial != NULL ) {
    *virial = vir;
  }
  if ( ep0 != NULL ) {
    *ep0 = ep;
  }
  return ep + lj->epot_tail; /* unshifted energy */
}



#define ljmix_force(lj) \
  lj->epot = ljmix_force_low(lj, lj->x, lj->f, \
      lj->r2ij, &lj->vir, &lj->ep0)

/* compute force and virial, return energy
 * the pair distances are recomputed */
__inline static double ljmix_force_low(ljmix_t *lj,
    double (*x)[D], double (*f)[D],
    double *r2ij, double *virial, double *ep0)
{
  double dx[D], fi[D], dr2, ir2, ir6, fs, ep, vir, rc2 = lj->rc2;
  double l = lj->l, invl = 1/l, sig, eps;
  int i, j, itp, jtp, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) {
    vzero(f[i]);
  }

  ep = vir = 0;
  for ( i = 0; i < n - 1; i++ ) {
    itp = lj->type[i];
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      jtp = lj->type[j];
      dr2 = ljmix_pbcdist2(dx, x[i], x[j], l, invl);
      if ( r2ij != NULL ) {
        r2ij[i*n + j] = dr2;
        r2ij[j*n + i] = dr2;
      }
      if ( dr2 >= rc2 ) continue;
      sig = lj->sigij[itp][jtp];
      eps = lj->epsij[itp][jtp];
      ir2 = (sig * sig) / dr2;
      ir6 = ir2 * ir2 * ir2;
      fs = eps * ir6 * (48 * ir6 - 24); /* f.r */
      vir += fs; /* f.r */
      fs /= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += eps * 4 * ir6 * (ir6 - 1);
      npr++;
    }
    vinc(f[i], fi);
  }
  if ( virial != NULL ) {
    *virial = vir;
  }
  if ( ep0 != NULL ) {
    *ep0 = ep;
  }
  return ep + lj->epot_tail; /* unshifted energy */
}



/* compute pressure */
static double ljmix_calcp(ljmix_t *lj, double tp)
{
  return (lj->dof * tp + lj->vir) / (D * lj->vol) + lj->p_tail;
}



/* velocity-verlet */
__inline static void ljmix_vv(ljmix_t *lj, double dt)
{
  int i, n = lj->n;
  double dth = dt * .5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(lj->v[i], lj->f[i], dth);
    vsinc(lj->x[i], lj->v[i], dt);
  }
  ljmix_force(lj);
  for (i = 0; i < n; i++) { /* VV part 2 */
    vsinc(lj->v[i], lj->f[i], dth);
  }
}



/* compute the kinetic energy */
#define ljmix_ekin(v, n) md_ekin(v, NULL, n)



/* exact velocity-rescaling thermostat */
#define ljmix_vrescale(lj, tp, dt) \
  md_vrescale(lj->v, NULL, lj->n, lj->dof, tp, dt)



/* position Langevin barostat, with coordinates only
 * NOTE: the first parameter is the degree of freedom
 * the scaling is r = r*s
 * set cutoff to half of the box */
__inline static void ljmix_langp0(ljmix_t *lj, double dt,
    double tp, double pext, int ensx)
{
  double pint, amp, s, dlnv;
  int i;

  pint = ljmix_calcp(lj, tp);
  amp = sqrt(2 * dt);
  dlnv = ((pint - pext) * lj->vol / tp + 1 - ensx) * dt + amp * randgaus();
  s = exp( dlnv / D );
  lj->vol *= exp( dlnv );
  ljmix_setvol(lj, lj->vol);
  for ( i = 0; i < lj->n; i++ ) {
    vsmul(lj->x[i], s);
  }
  ljmix_force(lj);
}



/* displace a random particle i, return i */
static int ljmix_randmv(ljmix_t *lj, double *xi, double amp)
{
  int i, d;

  i = (int) (rand01() * lj->n);
  for ( d = 0; d < D; d++ ) {
    xi[d] = lj->x[i][d] + (rand01() * 2 - 1) * amp;
  }
  return i;
}



/* compute pair energy */
__inline static int ljmix_pair(double dr2,
    double sig, double eps,
    double rc2, double *u, double *vir)
{
  if ( dr2 < rc2 ) {
    double invr2 = (sig * sig) / dr2;
    double invr6 = invr2 * invr2 * invr2;
    *vir = eps * invr6 * (48 * invr6 - 24); /* f.r */
    *u  = eps * 4 * invr6 * (invr6 - 1);
    return 1;
  } else {
    *vir = 0;
    *u = 0;
    return 0;
  }
}



/* return the energy change from displacing x[i] to xi */
__inline static double ljmix_depot(ljmix_t *lj, int i, double *xi, double *vir)
{
  int j, itp, jtp, n = lj->n;
  double l = lj->l, invl = 1/l, rc2 = lj->rc2, u, du, dvir;
  double dx[D], r2, sig, eps;

  u = 0.0;
  *vir = 0.0;
  itp = lj->type[i];
  for ( j = 0; j < n; j++ ) { /* pair */
    if ( j == i ) continue;
    jtp = lj->type[j];
    r2 = lj->r2ij[i*n + j];
    sig = lj->sigij[itp][jtp];
    eps = lj->epsij[itp][jtp];
    if ( ljmix_pair(r2, sig, eps, rc2, &du, &dvir) ) {
      u -= du;
      *vir -= dvir;
    }
    r2 = ljmix_pbcdist2(dx, xi, lj->x[j], l, invl);
    if ( ljmix_pair(r2, sig, eps, rc2, &du, &dvir) ) {
      u += du;
      *vir += dvir;
    }
    lj->r2i[j] = r2;
  }
  lj->r2i[i] = 0;
  return u;
}



/* commit a particle displacement */
__inline static void ljmix_commit(ljmix_t *lj, int i, const double *xi,
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
}



/* Metropolis algorithm */
__inline static int ljmix_metro(ljmix_t *lj, double amp, double bet)
{
  int i, acc = 0;
  double xi[D], r, du = 0, dvir = 0;

  i = ljmix_randmv(lj, xi, amp);
  du = ljmix_depot(lj, i, xi, &dvir);
  if ( du < 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp( -bet * du ) );
  }
  if ( acc ) {
    ljmix_commit(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}



/* wrap coordinates such that particles stay in the box */
__inline static int ljmix_wrapbox(ljmix_t *lj,
    double (*xin)[D], double (*xout)[D])
{
  int i, n = lj->n;
  double l = lj->l;

  for ( i = 0; i < n; i++ ) {
    vwrap( vcopy(xout[i], xin[i]), l );
  }
  return 0;
}



/* write positions (and possibly velocities) */
__inline static int ljmix_writepos(ljmix_t *lj,
    double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d, n = lj->n;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %.14e\n", D, n, (v != NULL), lj->l);
  for ( i = 0; i < n; i++ ) {
    for ( d = 0; d < D; d++ ) {
      fprintf(fp, "%.14e ", x[i][d]);
    }
    if ( v != NULL ) {
      for ( d = 0; d < D; d++ ) {
        fprintf(fp, "%.14e ", v[i][d]);
      }
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



#endif /* LJMIXCORE_H__ */
