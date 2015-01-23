/* implicit Gibbs-ensemble simulation */
#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif



int n = 108;
double rho = 0.12;
double tp = 1.32;
double rcdef = 1e9;
double amp = 0.2;
double rcls = 1.6;

double nsteps = 1e10;
double nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";

double mctot = 0, mcacc = 0;



int main(void)
{
  double t;
  lj_t *lj;
  int csize;
  double epsm = 0;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef, rcls);
  /* turn off the regular clustering energy */
  lj->lamcls = 0;
  /* turn on the division energy */
  lj->lamdiv = -1.5;

  /* compute the potential energy */
  lj_energy(lj);
  /* build a graph */
  lj_mkgraph(lj, lj->g);
  lj->ediv = lj_ediv(lj, lj->g, lj->cseed);

  for ( t = 1; t <= nsteps; t++ ) {
    /* a step of Metropolis algorithm
     * graph is implicitly computed */
    mctot += 1;
    mcacc += lj_metro(lj, amp, 1/tp);
    /* try to change the seed particle */
    lj_chseeddiv(lj, lj->g);
    //printf("cseed %d, %g, %g\n", lj->cseed, lj->ediv, lj_ediv(lj, lj->g, lj->cseed)); getchar();

    /* add the cluster size to the histogram */
    csize = graph_getcsize(lj->g, lj->cseed);
    lj_chist_add(lj, csize);

    if ( fmod(t, nstrep) < 0.1 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos, 0);
      lj_chist_save(lj, fnchist);
      fprintf(stderr, "t %g, mcacc %.2f%% ep %g, ediv %g, seed %d, csize %d, ",
          t, 100*mcacc/mctot, lj->epot, lj->ediv,
          lj->cseed, csize);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

