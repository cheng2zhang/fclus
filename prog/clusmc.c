/* Monte Carlo simulation that samples a flat histogram along the cluster size
 * The cluster under concern is the one that contain particle 0
 * Wang-Landau algorithm is used */
#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif



int n = 108;
double rho = 0.3;
double tp = 1.32;
double rcdef = 1e9;
double amp = 0.2;
double rcls = 1.6;

double nsteps = 1e10;
double nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";

double wl_flatness = 0.3;
double wl_frac = 0.5;

double tot = 0, acc = 0;



/* update the MC move size according to the acceptance ratio */
static void update_mcamp(double *nacc, double *ntot)
{
  double x = sqrt( *nacc / *ntot / 0.5 );
  if ( x > 2 ) x = 2;
  else if ( x < 0.5 ) x = 0.5;
  amp *= x;
  fprintf(stderr, "acc %g%%, amp %g\n", 100*(*nacc)/(*ntot), amp);
  *nacc = 0;
  *ntot = DBL_MIN;
}



int main(void)
{
  double t;
  lj_t *lj;
  double epsm = 0;
  double lnf = 0.001;
  int csize;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef, rcls);

  /* compute the potential energy */
  lj_energy(lj);
  /* build a graph */
  lj_mkgraph(lj, lj->g);

  for ( t = 1; t <= nsteps; t++ ) {
    /* a step of Metropolis algorithm
     * graph is implicitly computed */
    tot += 1;
    acc += lj_metro(lj, amp, 1/tp);
    /* try to change the seed particle */
    lj_changeseed(lj, lj->g);

    csize = graph_getcsize(lj->g, lj->cseed);
    /* add the cluster size to the histogram */
    lj_chist_add(lj, csize);
    /* update the adaptive potential */
    lj_update_vcls(lj, csize, lnf);
    /* change the updating magnitude */
    if ( lj_update_lnf(lj, &lnf, wl_flatness, wl_frac) != 0 ) {
      /* update the MC amplitude */
      update_mcamp(&acc, &tot);
    }

    if ( fmod(t, nstrep) < 0.1 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos);
      lj_chist_save(lj, fnchist);
      fprintf(stderr, "t %g, acc %.2f%% ep %g, seed %d, csize %d, lnf %g, flat %g%%, ",
          t, 100*acc/tot, lj->epot,
          lj->cseed, graph_getcsize(lj->g, 0),
          lnf, lj->hflatness*100);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

