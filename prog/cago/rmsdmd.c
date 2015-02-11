/* sample a flat histogram along the RMSD from the native structure */
#include "cagormsd.h"



const char *fnpdb = "pdb/1VII.pdb";
//const char *fnpdb = "pdb/1L2Y.pdb";
double kb = 200.0;
double ka = 40.0;
double kd1 = 1.0;
double kd3 = 0.5;
double nbe = 1.0;
double nbc = 4.0;
double rc = 6.0;

double tp = 1.1;
double mddt = 0.002;
double thdt = 0.1;
long nsteps = 10000000000L;
long nstrep = 100000;

double rmsdmin = 1.0;
double rmsdmax = 11.0;
//double rmsdmax = 6.0;
double rmsddel = 0.1; /* should be small enough */

const char *fnpos = "go.pos";
const char *fnvrmsd = "vrmsd.dat";

double wl_lnf0 = 1e-4;
double wl_flatness = 0.3;
double wl_frac = 0.5;
double invt_c = 5.0;

int nstblk = 10;
int nsthmc = 1;
int nvswaps = 1;



int main(void)
{
  cago_t *go;
  cagormsd_t *r;
  wl_t *wl;
  hmc_t *hmc;
  long t;
  double fdat[2], rmsd = 0;
  double hmcacc = 0, hmctot = DBL_MIN;

  /* open a Go model */
  if ( (go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, tp);

  /* open a Wang-Landau object */
  wl = wl_openf(rmsdmin, rmsdmax, rmsddel,
      wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  r = cagormsd_open(go, wl);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( nvswaps > 0 ) {
    go->dof = (go->n - 1) * D;
  }

  /* warm up the system such that RMSD > rmsdmin */
  for ( t = 1; rmsd < rmsdmin; t++ ) {
    cago_vv(go, 1.0, mddt);
    go->ekin = cago_vrescale(go, go->v, tp, thdt);
    rmsd = cago_rmsd(go, go->x, NULL);
  }
  fprintf(stderr, "warm-up t %ld, RMSD %g, ncont %d/%d\n",
      t, rmsd, cago_ncontacts(go, go->x, 1.2, NULL, NULL), go->ncont);

  /* make a hybrid Monte-Carlo object
   * and two extra floating-point numbers:
   * RMSD and potential energy */
  hmc = hmc_open(go->n, 0, 2);
  fdat[0] = rmsd;
  fdat[1] = go->epot;
  hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);

  for ( t = 1; t <= nsteps; t++ ) {
    cago_vv(go, 1.0, mddt);
    go->ekin = cago_vrescale(go, go->v, tp, thdt);

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( t % nsthmc == 0 ) {
      hmctot += 1;
      hmcacc += cagormsd_hmc(r, hmc, &rmsd);
      go->ekin = cago_vscramble(go, go->v, nvswaps);
    } else {
      rmsd = cago_rmsd(go, go->x, NULL);
    }

    /* Wang-Landau updating */
    wl_addf(wl, rmsd);

    if ( t % nstblk == 0 ) {
      wl_updatelnf(wl);
    }

    if ( t % nstrep == 0 ) {
      double flatness = wl_getflatness(wl);
      cago_writepos(go, go->x, go->v, fnpos);
      wl_save(wl, fnvrmsd);
      printf("%ld: ep %g, rmsd %g, hmcacc %.2f%%, flatness %.2f%%, lnf %g\n",
          t, go->epot, rmsd,
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
    }
  }

  cago_close(go);
  cagormsd_close(r);
  wl_close(wl);
  hmc_close(hmc);
  return 0;
}
