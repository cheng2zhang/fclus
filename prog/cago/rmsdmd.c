/* sample a flat histogram along the RMSD from the native structure */
#include "cagomodel.h"
#include "cagormsd.h"




/* use molecular dynamics to warm-up the system
 * till the RMSD reaches m->rmsdmin */
static void warmup_md_rmsd(cago_t *go, cagomodel_t *m)
{
  double t, rmsd = 0;

  /* warm up the system such that RMSD > rmsdmin */
  for ( t = 1; rmsd < m->rmsdmin; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);
    rmsd = cago_rmsd(go, go->x, NULL);
  }
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, cago_ncontacts(go, go->x, 1.2, NULL, NULL), go->ncont);
}



/* explicit HMC */
static int run_hmc_rmsd(cago_t *go, wl_t *wl, cagomodel_t *m)
{
  hmc_t *hmc;
  double t, rmsd;
  double hmcacc = 0, hmctot = DBL_MIN;

  /* make a hybrid Monte-Carlo object */
  hmc = cago_hmc_rmsd_init(go);

  for ( t = 1; t <= m->nsteps; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( fmod(t, m->nsthmc) < 0.1 ) {
      hmctot += 1;
      hmcacc += cago_hmc_rmsd(go, wl, hmc, &rmsd);
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

  hmc_close(hmc);

  return 0;
}



/* implicit HMC */
static int run_ihmc_rmsd(cago_t *go, wl_t *wl, cagomodel_t *m)
{
  hmc_t *hmc;
  double t, rmsd = 0;
  double hmcacc = 0, hmctot = DBL_MIN;
  double *fdat;

  /* make a hybrid Monte-Carlo object */
  hmc = cago_ihmc_rmsd_init(go, &fdat);

  for ( t = 1; t <= m->nsteps; t++ ) {
    hmctot += 1;
    hmcacc += cago_vv_rmsd(go, 1.0, m->mddt, go->x1,
        wl, m->mflmin, m->mflmax, m->mfhmin, m->mfhmax,
        m->temp, hmc, fdat);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvrmsd);
      rmsd = hmc->fdat[0];
      printf("%g: ep %g, rmsd %g, hmcacc %.2f%%, flatness %.2f%%, lnf %g\n",
          t, go->epot, rmsd,
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
    }
  }

  hmc_close(hmc);

  return 0;
}



int main(int argc, char **argv)
{
  cagomodel_t m[1];
  cago_t *go;
  wl_t *wl;

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

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved
   * NOTE: if the mass is not uniform, linear momenta are not
   * preserved as well, and `go->dof` should be `go->n * D` */
  if ( m->nvswaps > 0 ) {
    go->dof = (go->n - 1) * D;
  }

  /* warm-up the system */
  warmup_md_rmsd(go, m);

  if ( m->implicithmc ) {
    run_ihmc_rmsd(go, wl, m);
  } else {
    run_hmc_rmsd(go, wl, m);
  }

  cago_close(go);
  wl_close(wl);
  return 0;
}

