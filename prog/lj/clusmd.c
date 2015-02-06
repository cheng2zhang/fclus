/* molecular dynamics simulation with hybrid MC
 * to sample a flat histogram along the cluster size */
#include "ljcls.h"



int n = 55;
double rho = 0.25;
double tp = 2.0;
double rcdef = 1e9;
double dt = 0.002;
double thdt = 0.1;
double rcls = 1.6;

int nstblk = 10;
int nsthmc = 1; /* number of steps to do hybrid MC */
int nsteps = 1000000000;
int nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnvcls = "vcls.dat";

double wl_lnf0 = 4e-4;
double wl_flatness = 0.3;
double wl_frac = 0.5;
double invt_c = 1.0;

int changeseed = 0;
int nvswaps = 1;



int main(void)
{
  int t;
  lj_t *lj;
  ljcls_t *c;
  wl_t *wl;
  hmc_t *hmc;
  int idat[2], csize;
  double hmcacc = 0, hmctot = DBL_MIN;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef);

  /* make a Wang-Landau object */
  wl = wl_open(1, n + 1, wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  c = ljcls_open(lj, rcls, wl->v - 1);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( nvswaps > 0 ) {
    lj->dof = (n - 1) * D;
  }

  /* compute the force and potential energy */
  lj_force(lj);
  /* compute the cluster energy */
  ljcls_mkgraph(c, c->g);

  /* make a hybrid Monte-Carlo object
   * with two extra integers: cluster size and seed
   * and one extra floating-point: potential energy */
  hmc = hmc_open(n, 2, 1);
  idat[0] = graph_getcsize(c->g, c->cseed);
  idat[1] = c->cseed;
  hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);

  /* main molecular dynamics loop */
  for ( t = 1; t <= nsteps; t++ ) {
    /* velocity verlet */
    lj_vv(lj, dt);
    /* velocity-rescaling thermostat */
    lj->ekin = lj_vrescale(lj, tp, thdt);
    /* build a graph */
    ljcls_mkgraph(c, c->g);
    if ( changeseed ) {
      ljcls_changeseed(c, c->g);
    }

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( t % nsthmc == 0 ) {
      hmctot += 1;
      hmcacc += ljcls_hmc(c, hmc, &csize);
      lj->ekin = lj_vscramble(lj, lj->v, nvswaps);
    } else {
      csize = graph_getcsize(c->g, c->cseed);
    }

    /* add the cluster size to the histogram
     * and update the adaptive potential */
    wl_add(wl, csize);

    if ( t % nstblk == 0 ) {
      wl_updatelnf(wl);
    }

    if ( t % nstrep == 0 ) {
      double flatness = wl_getflatness(wl);
      ljcls_writepos(c, lj->x, lj->v, fnpos, 1);
      wl_save(wl, fnvcls);
      fprintf(stderr, "t %d, ep %g, ek %g, csize %d, hmcacc %.2f%%, flatness %.2f%%, lnf %g, ",
          t, lj->epot, lj->ekin, graph_getcsize(c->g, c->cseed),
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
      graph_clus_print(c->g);
    }
  }
  lj_close(lj);
  wl_close(wl);
  ljcls_close(c);
  hmc_close(hmc);
  fprintf(stderr, "rho %g, tp %g\n", rho, tp);
  return 0;
}

