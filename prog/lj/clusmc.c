/* Monte Carlo simulation that samples a flat histogram along the cluster size
 * Wang-Landau algorithm is used */
#include "ljcls.h"



int n = 108;
double rho = 0.25;
double tp = 1.32;
double rcdef = 1e9;
double amp = 0.2;
double rcls = 1.6;

int nstblk = 10;
double nsteps = 1e10;
double nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnvcls = "vcls.dat";

double wl_lnf0 = 0.001;
double wl_flatness = 0.3;
double wl_frac = 0.5;
double invt_c = 1.0;

double tot = 0, acc = 0;



int main(void)
{
  lj_t *lj;
  ljcls_t *c;
  wl_t *wl;
  int csize, it;
  double t;

  /* open a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef);

  /* open a Wang-Landau object */
  wl = wl_openi(1, n + 1, wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  c = ljcls_open(lj, rcls, wl->v - 1);

  /* compute the potential energy */
  lj_energy(lj);
  /* build a graph */
  ljcls_mkgraph(c, c->g);

  for ( t = nstblk; t <= nsteps; t += nstblk ) {
    for ( it = 0; it < nstblk; it++ ) {
      /* a step of Metropolis algorithm
       * graph is implicitly computed */
      tot += 1;
      acc += ljcls_metro(c, amp, 1/tp);

      /* try to change the seed particle */
      ljcls_changeseed(c, c->g);

      csize = graph_getcsize(c->g, c->cseed);
      /* add the cluster size to the histogram
       * and update the adaptive potential */
      wl_addi(wl, csize);
    }

    if ( wl_updatelnf(wl) ) {
      /* update the MC amplitude */
      update_mcamp(&amp, 0.5, &acc, &tot);
    }

    if ( fmod(t, nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      ljcls_writepos(c, lj->x, lj->v, fnpos, 1);
      wl_save(wl, fnvcls);
      fprintf(stderr, "t %g, acc %.2f%% ep %g, seed %d, csize %d, flatness %g%%, lnf %g ",
          t, 100 * acc / tot, lj->epot,
          c->cseed, csize,
          100.0 * flatness, wl->lnf);
      graph_clus_print(c->g);
    }
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  ljcls_close(c);
  wl_close(wl);
  fprintf(stderr, "rho %g, tp %g\n", rho, tp);
  return 0;
}

