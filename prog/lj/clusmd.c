/* molecular dynamics simulation with hybrid MC
 * to sample a flat histogram along the cluster size */
#include "ljmodel.h"
#include "ljcls.h"



static void ljmodel_default_md(ljmodel_t *m)
{
  ljmodel_default(m);
  m->n = 55;
  m->wl_lnf0 = 4e-4;
}



int main(int argc, char **argv)
{
  ljmodel_t m[1];
  lj_t *lj;
  ljcls_t *c;
  wl_t *wl;
  hmc_t *hmc;
  int idat[2], csize;
  double t;
  double hmcacc = 0, hmctot = DBL_MIN;

  ljmodel_default_md(m);
  ljmodel_doargs(m, argc, argv);

  /* make a Lennard-Jones object */
  lj = lj_open(m->n, m->rho, m->rcdef);

  /* make a Wang-Landau object */
  wl = wl_openi(1, m->n + 1,
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);
  c = ljcls_open(lj, m->rcls, wl->v - 1);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( m->nvswaps > 0 ) {
    lj->dof = (m->n - 1) * D;
  }

  /* compute the force and potential energy */
  lj_force(lj);
  /* compute the cluster energy */
  ljcls_mkgraph(c, c->g);

  /* make a hybrid Monte-Carlo object
   * with two extra integers: cluster size and seed
   * and one extra floating-point: potential energy */
  hmc = hmc_open(m->n, 2, 1);
  idat[0] = graph_getcsize(c->g, c->cseed);
  idat[1] = c->cseed;
  hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);

  /* main molecular dynamics loop */
  for ( t = 1; t <= m->nsteps; t++ ) {
    /* velocity verlet */
    lj_vv(lj, m->mddt);
    /* velocity-rescaling thermostat */
    lj->ekin = lj_vrescale(lj, m->temp, m->thdt);
    /* build a graph */
    ljcls_mkgraph(c, c->g);
    if ( m->changeseed ) {
      ljcls_changeseed(c, c->g);
    }

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( fmod(t, m->nsthmc) < 0.1 ) {
      hmctot += 1;
      hmcacc += ljcls_hmc(c, hmc, &csize);
      lj->ekin = lj_vscramble(lj, lj->v, m->nvswaps);
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
      ljcls_writepos(c, lj->x, lj->v, m->fnpos, 1);
      wl_save(wl, m->fnvcls);
      fprintf(stderr, "t %g, ep %g, ek %g, csize %d, hmcacc %.2f%%, flatness %.2f%%, lnf %g, ",
          t, lj->epot, lj->ekin, graph_getcsize(c->g, c->cseed),
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
      graph_clus_print(c->g);
    }
  }
  lj_close(lj);
  ljcls_close(c);
  wl_close(wl);
  hmc_close(hmc);
  return 0;
}

