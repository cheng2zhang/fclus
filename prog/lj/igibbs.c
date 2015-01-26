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
double rho = 0.3;
double tp = 1.2;
double rcdef = 1e9;
double amp = 0.2;
double rcls = 1.6;

double nsteps = 1e10;
double nstrep = 10000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";

double mctot = 0, mcacc = 0;
double divtot = 0, divacc = 0;



int main(void)
{
  double t;
  lj_t *lj;
  int divsize;
  double epsm = 0;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef, rcls);
  /* turn off the regular clustering energy */
  lj->lamcls = 0;
  /* turn on the division energy */
  lj->lamdiv = 0.00001;

  /* compute the potential energy */
  lj_energy(lj);
  /* build a graph */
  lj_mkgraph(lj, lj->g);
  lj->ediv = lj_ediv(lj, lj->idiv, lj->divseed, lj->divr);

  for ( t = 1; t <= nsteps; t++ ) {
    /* a step of Metropolis algorithm
     * graph is implicitly computed */
    mctot += 1;
    mcacc += lj_metro(lj, amp, 1/tp);
    /* try to change the division */
    divtot += 1;
    divacc += lj_changediv(lj, lj->l * 0.1, 1.0/tp);
    //printf("%d, %g, %g\n", lj_getdivsize(lj, lj->idiv), lj->ediv, lj_ediv(lj, lj->idiv)); getchar();

    /* add the division size to histogram */
    divsize = lj_getdivsize(lj, lj->idiv);
    lj_chist_add(lj, divsize);

    if ( fmod(t, nstrep) < 0.1 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos, 0);
      lj_chist_save(lj, fnchist);
      fprintf(stderr, "t %g, mcacc %.2f%% ep %g, divacc %.2f%%, ediv %g, "
          "divsize %d, divr %g, divrho %g\n",
          t, 100*mcacc/mctot, lj->epot,
          100*divacc/divtot, lj->ediv,
          divsize, lj->divr, divsize/lj_sphrvol(lj->divr));
    }
    epsm += lj->epot;
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}
