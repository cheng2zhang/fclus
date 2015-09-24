/* sample a flat histogram along the number of contacts */
#include "cagomodel.h"
#include "cagonc.h"



/* use molecular dynamics to warm-up the system
 * till the NC is lowered below m->ncmax */
static void warmup_md_nc(cago_t *go, cagomodel_t *m,
    int ncmin, int ncmax)
{
  int nc;
  double t, rmsd = 0;

  /* warm up the system such that nc < ncmax */
  nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  for ( t = 1; nc <= ncmin || nc >= ncmax; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);
    nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  }
  rmsd = cago_rmsd(go, go->x, NULL);
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, nc, go->ncont);
}



/* explicit HMC */
static int run_hmc_nc(cago_t *go, wl_t *wl, cagomodel_t *m)
{
  hmc_t *hmc;
  int nc;
  double t, rmsd;
  double hmcacc = 0, hmctot = DBL_MIN;

  /* make a hybrid Monte-Carlo object
   * with one extra integer: number of contacts
   * and one extra floating-point: potential energy */
  hmc = hmc_open(go->n, 1, 1);
  hmc_push(hmc, go->x, go->v, go->f, &nc, &go->epot);

  for ( t = 1; t <= m->nsteps; t++ ) {
    cago_vv(go, 1.0, m->mddt);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);

    /* use hybrid MC to sample a flat histogram along the cluster size */
    if ( fmod(t, m->nsthmc) < 0.1 ) {
      hmctot += 1;
      hmcacc += cago_hmc_nc(go, wl, hmc, &nc);
      go->ekin = cago_vscramble(go, go->v, m->nvswaps);
    } else {
      nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
    }

    /* Wang-Landau updating */
    wl_addi(wl, nc);

    if ( fmod(t, m->nstblk) < 0.1 ) {
      wl_updatelnf(wl);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvnc);
      rmsd = cago_rmsd(go, go->x, NULL);
      printf("%g: ep %g, nc %d, rmsd %g, hmcacc %.2f%%, flatness %.2f%%, lnf %g\n",
          t, go->epot, nc, rmsd,
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
    }
  }

  hmc_close(hmc);

  return 0;
}



/* implicit HMC */
static int run_ihmc_nc(cago_t *go, wl_t *wl, cagomodel_t *m)
{
  hmc_t *hmc;
  double t;
  int idat = 0;
  double hmcacc = 0, hmctot = DBL_MIN;
  double *fdat;
  xymap_t *xy;

  /* make a hybrid Monte-Carlo object */
  hmc = cago_ihmc_nc_init(go, &idat, &fdat);
  xy = xymap_open(-go->ncont, 0, 1.0);

  for ( t = 1; t <= m->nsteps; t++ ) {
    hmctot += 1;
    hmcacc += cago_vv_nc(go, 1.0, m->mddt, go->x1,
        wl, xy, m->mflmin, m->mflmax, m->mfhmin, m->mfhmax,
        m->temp, hmc, &idat, fdat);
    go->ekin = cago_vrescale(go, go->v, m->temp, m->thdt);

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvrmsd);
      printf("%g: ep %g, nc %d, hmcacc %.2f%%, flatness %.2f%%, lnf %g\n",
          t, go->epot, hmc->idat[0],
          100.0 * hmcacc / hmctot, 100.0 * flatness, wl->lnf);
    }
  }

  hmc_close(hmc);
  xymap_close(xy);

  return 0;
}



int main(int argc, char **argv)
{
  cagomodel_t m[1];
  cago_t *go;
  wl_t *wl;
  int ncmin, ncmax;

  cagomodel_default(m);
  m->fncmin = 0.2;
  m->fncmax = 0.9;
  m->invt_c = 10.0;
  m->nstrep = 100000;
  cagomodel_doargs(m, argc, argv);

  /* open a Go model */
  if ( (go = cago_open(m->fnpdb, m->kb, m->ka, m->kd1, m->kd3, m->nbe, m->nbc, m->rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", m->fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, m->temp);

  ncmin = (int) (m->fncmin * go->ncont + 0.5);
  ncmax = (int) (m->fncmax * go->ncont + 0.5);

  /* open a Wang-Landau object */
  wl = wl_openi(ncmin, ncmax,
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);

  /* change the degrees of freedom, with velocity swaps
   * the angular momenta are no longer conserved */
  if ( m->nvswaps > 0 ) {
    go->dof = (go->n - 1) * D;
  }

  /* warm up the system such that nc < ncmax */
  warmup_md_nc(go, m, ncmin, ncmax);

  if ( m->implicithmc ) {
    run_ihmc_nc(go, wl, m);
  } else {
    run_hmc_nc(go, wl, m);
  }

  cago_close(go);
  wl_close(wl);
  return 0;
}

