/* Monte Carlo simulation that samples a flat histogram along the cluster size
 * Wang-Landau algorithm is used */
#include "ljmixmodel.h"
#include "ljmixcls.h"



int main(int argc, char **argv)
{
  ljmixmodel_t m[1];
  ljmix_t *lj;
  ljmixcls_t *c;
  wl_t *wl;
  int csize, it;
  double t;
  double tot = 0, acc = 0;

  ljmixmodel_default(m);
  ljmixmodel_doargs(m, argc, argv);

  /* open a Lennard-Jones object */
  lj = ljmix_open(m->ns, m->np, m->sig, m->rho, m->rcdef);

  /* open a Wang-Landau object */
  wl = wl_openi(1, m->np[0],
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);
  c = ljmixcls_open(lj, m->rcls, wl->v - 1);

  /* compute the potential energy */
  ljmix_energy(lj);
  /* build a graph */
  ljmixcls_mkgraph(c, c->g);

  for ( t = m->nstblk; t <= m->nsteps; t += m->nstblk ) {
    for ( it = 0; it < m->nstblk; it++ ) {
      /* a step of Metropolis algorithm
       * graph is implicitly computed */
      tot += 1;
      acc += ljmixcls_metro(c, m->mcamp, m->beta);

      /* try to change the seed particle */
      ljmixcls_changeseed(c, c->g);

      csize = graph_getcsize(c->g, c->cseed);
      /* add the cluster size to the histogram
       * and update the adaptive potential */
      wl_addi(wl, csize);
    }

    if ( wl_updatelnf(wl) ) {
      /* update the MC amplitude */
      update_mcamp(&m->mcamp, 0.5, &acc, &tot);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      ljmixcls_writepos(c, lj->x, lj->v, m->fnpos, 1);
      wl_save(wl, m->fnvcls);
      fprintf(stderr, "t %g, acc %.2f%% ep %g, seed %d, csize %d, flatness %g%%, lnf %g ",
          t, 100 * acc / tot, lj->epot,
          c->cseed, csize,
          100.0 * flatness, wl->lnf);
      graph_clus_print(c->g);
    }
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  ljmix_close(lj);
  ljmixcls_close(c);
  wl_close(wl);
  return 0;
}

