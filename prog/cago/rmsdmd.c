/* sample a flat histogram along the RMSD from the native structure */
#include "cagomodel.h"
#include "cagormsd.h"



int main(int argc, char **argv)
{
  cagomodel_t m[1];
  cago_t *go;
  cagormsd_t *r;
  wl_t *wl;
  hmc_t *hmc;
  double t;
  double fdat[2], rmsd = 0;
  double hmcacc = 0, hmctot = DBL_MIN;

  cagomodel_default(m);
  m->invt_c = 5.0;
  m->nstrep = 100000;
  cagomodel_doargs(m, argc, argv);

  /* open a Go model */
  if ( (go = cago_open(m->fnpdb, m->kb, m->ka, m->kd1, m->kd3, m->nbe, m->nbc, m->rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", m->fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, m->temp);

  /* open a Wang-Landau object */
  wl = wl_openf(m->rmsdmin, m->rmsdmax, m->rmsddel,
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);
  r = cagormsd_open(go, wl);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( m->nvswaps > 0 ) {
    go->dof = (go->n - 1) * D;
  }

  /* warm up the system such that RMSD > rmsdmin */
  for ( t = 1; rmsd < m->rmsdmin; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);
    rmsd = cago_rmsd(go, go->x, NULL);
  }
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, cago_ncontacts(go, go->x, 1.2, NULL, NULL), go->ncont);

  /* make a hybrid Monte-Carlo object
   * and two extra floating-point numbers:
   * RMSD and potential energy */
  hmc = hmc_open(go->n, 0, 2);
  fdat[0] = rmsd;
  fdat[1] = go->epot;
  hmc_push(hmc, go->x, go->v, go->f, NULL, fdat);

  for ( t = 1; t <= m->nsteps; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( fmod(t, m->nsthmc) < 0.1 ) {
      hmctot += 1;
      hmcacc += cagormsd_hmc(r, hmc, &rmsd);
      go->ekin = cago_vscramble(go, go->v, m->nvswaps);
    } else {
      rmsd = cago_rmsd(go, go->x, NULL);
    }

    /* Wang-Landau updating */
    wl_addf(wl, rmsd);

    if ( fmod(t, m->nstblk) < 0.1 ) {
      wl_updatelnf(wl);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvrmsd);
      printf("%g: ep %g, rmsd %g, hmcacc %.2f%%, flatness %.2f%%, lnf %g\n",
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
