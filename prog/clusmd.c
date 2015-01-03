/* basic molecular dynamics simulation in the NVT or NVE ensemble */
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
double dt = 0.002;
double thdt = 0.02;
double rcls = 1.6;

int nequil = 1000;
int nsteps = 10000;
int nstrep = 1000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";



int main(void)
{
  int t;
  lj_t *lj;
  double epsm = 0;

  lj = lj_open(n, rho, rcdef, rcls);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    lj_vv(lj, dt);
    lj->ekin = lj_vrescale(lj, tp, thdt);
    if ( t <= nequil ) continue;
    lj_clus(lj, lj->g, rcls);
    if ( t % nstrep == 0 ) {
      lj_chist_print(lj);
      printf("%d, ep %g, ek %g, ", t, lj->epot, lj->ekin);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_chist_save(lj, fnchist);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

