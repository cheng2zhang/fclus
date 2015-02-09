#include "cagormsd.h"



//const char *fnpdb = "pdb/1VII.pdb";
const char *fnpdb = "pdb/1L2Y.pdb";
double kb = 200.0;
double ka = 40.0;
double kd1 = 1.0;
double kd3 = 0.5;
double nbe = 1.0;
double nbc = 4.0;
double rc = 6.0;

//double tp = 1.1;
double tp = 1.0;
double amp = 0.2;
double nsteps = 1e10;
double nstrep = 1000000;

double rmsdmin = 1.0;
//double rmsdmax = 11.0;
//double rmsddel = 0.05;
double rmsdmax = 6.0;
double rmsddel = 0.1; /* should be small enough */

const char *fnpos = "go.pos";
const char *fnvrmsd = "vrmsd.dat";

double wl_lnf0 = 4e-4;
double wl_flatness = 0.3;
double wl_frac = 0.5;
double invt_c = 1.0;

int nstblk = 10;
int nsthmc = 1;
int nvswaps = 1;

double tot = 0, acc = 0;



int main(void)
{
  cago_t *go;
  cagormsd_t *r;
  wl_t *wl;
  double rmsd = 0;
  int it;
  double t;

  if ( (go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, tp);
  go->dof = go->n * D;

  /* open a Wang-Landau object */
  wl = wl_openf(rmsdmin, rmsdmax, rmsddel,
      wl_lnf0, wl_flatness, wl_frac, invt_c, 0);
  r = cagormsd_open(go, wl);

  /* warm up the system such that RMSD > rmsdmin */
  for ( t = nstblk; rmsd < rmsdmin; t += nstblk ) {
    for ( it = 0; it < nstblk; it++ ) {
      cago_metro(go, amp, 1/tp);
    }
    rmsd = cago_rmsd(go, go->x, NULL);
  }
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, cago_ncontacts(go, go->x, 1.2, NULL, NULL), go->ncont);

  for ( t = nstblk; t <= nsteps; t += nstblk ) {
    for ( it = 0; it < nstblk; it++ ) {
      /* a step of Metropolis alogrithm */
      tot += 1;
      acc += cagormsd_metro(r, amp, 1/tp, &rmsd);

      wl_addf(wl, rmsd);
    }

    if ( wl_updatelnf(wl) ) {
      /* update the MC amplitude */
      update_mcamp(&amp, 0.5, &acc, &tot);
    }

    if ( fmod(t, nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      int nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
      cago_writepos(go, go->x, go->v, fnpos);
      wl_save(wl, fnvrmsd);
      fprintf(stderr, "%g: acc %.2f%%, ep %g, rmsd %g, nc %d/%d, flatness %g%%, lnf %g\n",
          t, 100.0 * acc / tot, go->epot,
          rmsd, nc, go->ncont, 100.0 * flatness, wl->lnf);
      cago_rmcom(go, go->x, go->v);
    }
  }

  cago_close(go);
  cagormsd_close(r);
  wl_close(wl);
  return 0;
}
