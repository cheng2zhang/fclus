/* basic Monte Carlo simulation */
#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif



int n = 256;
double rho = 0.3;
double tp = 1.32;
double rcdef = 1e9;
double amp = 0.2;
double rcls = 1.6;

int nequil = 100000;
int nsteps = 1000000;
int nstcls = 100;
int nstrep = 100000;

const char *fnpos = "lj.pos";



int main(void)
{
  int t;
  lj_t *lj;
  double epsm = 0;

  lj = lj_open(n, rho, rcdef);
  lj_energy(lj);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_metro(lj, amp, 1/tp);
    if ( t <= nequil ) continue;
    if ( t % nstcls == 0 )
      lj_clus(lj, lj->g, rcls);
    if ( t % nstrep == 0 ) {
      graph_chist_print(lj->g);
      printf("%d, ep %g\n", t, lj->epot);
    }
    epsm += lj->epot;
  }
  printf("%g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

