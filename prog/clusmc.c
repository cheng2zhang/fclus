/* basic Monte Carlo simulation */
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

int nequil = 10000;
int nsteps = 2000000;
int nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";


int main(void)
{
  int t;
  lj_t *lj;
  double epsm = 0;

  lj = lj_open(n, rho, rcdef, rcls);
  /* potential energy */
  lj_energy(lj);
  /* cluster energy */
  lj_mkgraph(lj, lj->g, rcls);
  lj->ecls = lj_eclus(lj, lj->g);

  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_metro(lj, amp, 1/tp);
    if ( t <= nequil ) continue;
    lj_chist_add(lj, lj->g);
    if ( t % nstrep == 0 ) {
      lj_chist_save(lj, fnchist);
      printf("%d, ep %g, ecls %g, ", t, lj->epot, lj->ecls);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

