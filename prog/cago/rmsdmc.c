/* sample a flat histogram along the RMSD from the native structure */
#include "cagomodel.h"
#include "cagormsd.h"



static void warmup_mc_rmsd(cago_t *go, cagomodel_t *m)
{
  int it;
  double t, rmsd = 0;

  /* warm up the system such that RMSD > rmsdmin */
  for ( t = m->nstblk; rmsd < m->rmsdmin; t += m->nstblk ) {
    for ( it = 0; it < m->nstblk; it++ ) {
      cago_metro(go, m->mcamp, 1/m->temp);
    }
    rmsd = cago_rmsd(go, go->x, NULL);
  }
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, cago_ncontacts(go, go->x, 1.2, NULL, NULL), go->ncont);
}



static int run_metro_rmsd(cago_t *go, wl_t *wl, cagomodel_t *m)
{
  int it;
  double t, rmsd = 0, tot = 0, acc = 0;

  for ( t = m->nstblk; t <= m->nsteps; t += m->nstblk ) {
    for ( it = 0; it < m->nstblk; it++ ) {
      /* a step of Metropolis alogrithm */
      tot += 1;
      acc += cago_metro_rmsd(go, wl, m->mcamp, 1/m->temp, &rmsd);

      wl_addf(wl, rmsd);
    }

    if ( wl_updatelnf(wl) ) {
      /* update the MC amplitude */
      update_mcamp(&m->mcamp, 0.5, &acc, &tot);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      int nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvrmsd);
      fprintf(stderr, "%g: acc %.2f%%, ep %g, rmsd %g, nc %d/%d, flatness %g%%, lnf %g\n",
          t, 100.0 * acc / tot, go->epot,
          rmsd, nc, go->ncont, 100.0 * flatness, wl->lnf);
      cago_rmcom(go, go->x, go->v);
    }
  }

  return 0;
}



int main(int argc, char **argv)
{
  cagomodel_t m[1];
  cago_t *go;
  wl_t *wl;

  cagomodel_default(m);
  cagomodel_doargs(m, argc, argv);

  /* open a Go model */
  if ( (go = cago_open(m->fnpdb, m->kb, m->ka, m->kd1, m->kd3, m->nbe, m->nbc, m->rc,
                       PDB_CONTACT_HEAVY, 4)) == NULL ) {
    fprintf(stderr, "cannot initialize Go model from %s\n", m->fnpdb);
    return -1;
  }
  cago_initmd(go, 0, 0.01, m->temp);
  go->dof = go->n * D;

  /* open a Wang-Landau object */
  wl = wl_openf(m->rmsdmin, m->rmsdmax, m->rmsddel,
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, NULL, 0);

  warmup_mc_rmsd(go, m);

  run_metro_rmsd(go, wl, m);

  cago_close(go);
  wl_close(wl);
  return 0;
}
