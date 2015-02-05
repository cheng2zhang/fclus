/* molecular dynamics simulation with hybrid MC
 * to sample a flat histogram along the cluster size */
#define D 3
#include "ljcore.h"



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
  wl_t *wl;
  double epsm = 0;
  hmc_t *hmc;
  int idat[2], csize;
  double hmcacc = 0, hmctot = DBL_MIN;

  /* make a Wang-Landau object */
  wl = wl_open(1, n + 1, wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef, rcls, wl->v - 1);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( nvswaps > 0 ) {
    lj->dof = (n - 1) * D;
  }

  /* compute the force and potential energy */
  lj_force(lj);
  /* compute the cluster energy */
  lj_mkgraph(lj, lj->g);

  /* make a hybrid Monte-Carlo object
   * with two extra integers: cluster size and seed
   * and one extra floating-point: potential energy */
  hmc = hmc_open(n, 2, 1);
  idat[0] = graph_getcsize(lj->g, lj->cseed);
  idat[1] = lj->cseed;
  hmc_push(hmc, lj->x, lj->v, lj->f, idat, &lj->epot);

  /* main molecular dynamics loop */
  for ( t = 1; t <= nsteps; t++ ) {
    /* velocity verlet */
    lj_vv(lj, dt);
    /* velocity-rescaling thermostat */
    lj->ekin = lj_vrescale(lj, tp, thdt);
    /* build a graph */
    lj_mkgraph(lj, lj->g);
    if ( changeseed ) {
      lj_changeseed(lj, lj->g);
    }

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( t % nsthmc == 0 ) {
      hmctot += 1;
      hmcacc += lj_hmc(lj, hmc, &csize);
      lj->ekin = lj_vscramble(lj->v, lj->n, nvswaps);
    }

    csize = graph_getcsize(lj->g, lj->cseed);
    /* add the cluster size to the histogram
     * and update the adaptive potential */
    wl_add(wl, csize);

    if ( t % nstblk == 0 ) {
      wl_updatelnf(wl);
    }

    if ( t % nstrep == 0 ) {
      double flatness = wl_getflatness(wl);
      lj_writepos(lj, lj->x, lj->v, fnpos, 1);
      wl_save(wl, fnvcls);
      fprintf(stderr, "t %d, ep %g, ek %g, csize %d, hmcacc %.2f%%, flatness %.2f%%, lnf %g, ",
          t, lj->epot, lj->ekin, graph_getcsize(lj->g, lj->cseed),
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  lj_close(lj);
  hmc_close(hmc);
  wl_close(wl);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

