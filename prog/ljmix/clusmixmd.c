/* molecular dynamics simulation with hybrid MC
 * to sample a flat histogram along the cluster size */
#include "ljmixmodel.h"
#include "ljmixcls.h"



int main(int argc, char **argv)
{
  ljmixmodel_t m[1];
  ljmix_t *lj;
  ljmixcls_t *c;
  wl_t *wl;
  hmc_t *hmc;
  int idat[2], csize;
  double t;
  double hmcacc = 0, hmctot = DBL_MIN;

  ljmixmodel_default(m);
  ljmixmodel_doargs(m, argc, argv);

  /* make a Lennard-Jones object */
  lj = ljmix_open(m->ns, m->np, m->sig, m->eps, m->rho, m->rcdef);

  /* make a Wang-Landau object */
  wl = wl_openi(1, m->np[0],
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, NULL, 0);
  c = ljmixcls_open(lj, m->rcls, wl->v - 1);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( m->nvswaps > 0 ) {
    lj->dof = (lj->n - 1) * D;
  }

  /* compute the force and potential energy */
  ljmix_force(lj);
  /* compute the cluster energy */
  ljmixcls_mkgraph(c, c->g);

  /* make a hybrid Monte-Carlo object
   * with two extra integers: cluster size and seed
   * and one extra floating-point: potential energy */
  hmc = hmc_open(lj->n, 2, 1);
  idat[0] = graph_getcsize(c->g, c->cseed);
  idat[1] = c->cseed;
  hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);

  /* main molecular dynamics loop */
  for ( t = 1; t <= m->nsteps; t++ ) {
    /* velocity verlet */
    ljmix_vv(lj, m->mddt);
    /* velocity-rescaling thermostat */
    lj->ekin = ljmix_vrescale(lj, m->temp, m->thdt);
    /* build a graph */
    ljmixcls_mkgraph(c, c->g);
    if ( m->changeseed ) {
      ljmixcls_changeseed(c, c->g);
    }

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( fmod(t, m->nsthmc) < 0.1 ) {
      hmctot += 1;
      hmcacc += ljmixcls_hmc(c, hmc, &csize);
      lj->ekin = ljmix_vscramble(lj, lj->v, m->nvswaps);
    } else {
      csize = graph_getcsize(c->g, c->cseed);
    }

    /* add the cluster size to the histogram
     * and update the adaptive potential */
    wl_addi(wl, csize);

    if ( fmod(t, m->nstblk) < 0.1 ) {
      wl_updatelnf(wl);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      ljmixcls_writepos(c, lj->x, lj->v, m->fnpos, 1);
      wl_save(wl, m->fnvcls);
      fprintf(stderr, "t %g, ep %g, ek %g, csize %d, hmcacc %.2f%%, flatness %.2f%%, lnf %g, ",
          t, lj->epot, lj->ekin, graph_getcsize(c->g, c->cseed),
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
      graph_clus_print(c->g);
    }
  }
  ljmix_close(lj);
  ljmixcls_close(c);
  wl_close(wl);
  hmc_close(hmc);
  return 0;
}

