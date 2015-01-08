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

double flatness = 0.3;
double frac = 0.5;



int main(void)
{
  double t;
  lj_t *lj;
  double epsm = 0;
  double lnf = 0.001 * tp;

  lj = lj_open(n, rho, rcdef, rcls);
  /* potential energy */
  lj_energy(lj);
  /* cluster energy */
  lj_mkgraph(lj, lj->g, rcls);
  lj->ecls = lj_eclus(lj, lj->g);
  for ( t = 1; t <= nsteps; t++ ) {
    lj_metro(lj, amp, 1/tp);

    lj_chist_add(lj, lj->g);
    lj_update_vcls(lj, lj->g, lnf);
    //lj_chist_save(lj, fnchist); getchar();
    lnf = lj_update_lnf(lj, lnf, flatness, frac);

    if ( fmod(t, nstrep) < 0.1 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos);
      lj_chist_save(lj, fnchist);
      printf("%g, ep %g, ecls %g, lnf %g, flat %g%%, ",
          t, lj->epot, lj->ecls, lnf, lj->hflatness*100);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  printf("rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

