/* basic molecular dynamics simulation in the NVT or NVE ensemble */
#ifndef D
#define D 3
#endif

#if D == 3
#include "lj3d.h"
#else
#include "lj2d.h"
#endif

#include "hmc.h"



int n = 55;
double rho = 0.3;
double tp = 2.0;
double rcdef = 1e9;
double dt = 0.002;
double thdt = 0.1;
double rcls = 1.6;

int nsthmc = 1; /* number of steps to do hybrid MC */
int nsteps = 1000000000;
int nstrep = 1000000;

const char *fnpos = "lj.pos";
const char *fnchist = "chist.dat";

double wl_flatness = 0.3;
double wl_frac = 0.5;

int nvswaps = 1;



int main(void)
{
  int t;
  lj_t *lj;
  double epsm = 0;
  hmc_t *hmc;
  int hmc_csize = 0;
  double hmcacc = 0, hmctot = DBL_MIN;
  double lnf = 1e-4;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef, rcls);
  if ( nvswaps > 0 ) {
    lj->dof = (n - 1) * D;
  }

  /* compute the force and potential energy */
  lj_force(lj);
  /* compute the cluster energy */
  lj_mkgraph(lj, lj->g);

  /* make a hybrid Monte-Carlo object */
  hmc = hmc_open(n);
  hmc_push(hmc, lj->x, lj->v, lj->f);
  hmc_csize = graph_getcsize(lj->g, lj->cseed);

  for ( t = 1; t <= nsteps; t++ ) {
    /* velocity verlet */
    lj_vv(lj, dt);
    /* velocity-rescaling thermostat */
    lj->ekin = lj_vrescale(lj, tp, thdt);
    /* build a graph */
    lj_mkgraph(lj, lj->g);
    lj_changeseed(lj, lj->g);

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( t % nsthmc == 0 ) {
      hmctot += 1;
      hmcacc += lj_hmc(lj, hmc, &hmc_csize, nvswaps);
    }

    /* add the cluster size to the histogram */
    lj_chist_add(lj, lj->g);
    /* update the adaptive potential */
    lj_update_vcls(lj, lj->g, lnf);
    /* change the updating magnitude */
    lj_update_lnf(lj, &lnf, wl_flatness, wl_frac);

    if ( t % nstrep == 0 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos);
      //lj_chist_print(lj);
      lj_chist_save(lj, fnchist);
      fprintf(stderr, "t %d, ep %g, ek %g, csize %d, hmcacc %.2f%%, flatness %.2f%%, ",
          t, lj->epot, lj->ekin, graph_getcsize(lj->g, 0),
          100.0 * hmcacc / hmctot, 100.0 * lj->hflatness );
      graph_clus_print(lj->g);
    }
    epsm += lj->epot;
  }
  lj_close(lj);
  hmc_close(hmc);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

