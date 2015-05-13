/* sample a flat histogram along the number of native contacts */
#include "cagomodel.h"
#include "cagonc.h"





int main(int argc, char **argv)
{
  cagomodel_t m[1];
  cago_t *go;
  cagonc_t *r;
  wl_t *wl;
  int it, nc, ncmin, ncmax;
  double t, rmsd = 0, tot = 0, acc = 0;

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

  ncmin = (int) (m->fncmin * go->ncont + 0.5);
  ncmax = (int) (m->fncmax * go->ncont + 0.5);

  /* open a Wang-Landau object */
  wl = wl_openi(ncmin, ncmax,
      m->wl_lnf0, m->wl_flatness, m->wl_frac, m->invt_c, 0);
  r = cagonc_open(go, wl);

  /* warm up the system such that nc < ncmax */
  nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  for ( t = m->nstblk; nc >= ncmax; t += m->nstblk ) {
    for ( it = 0; it < m->nstblk; it++ ) {
      cago_metro(go, m->mcamp, 1/m->temp);
    }
    nc = cago_ncontacts(go, go->x, -1, NULL, NULL);
  }
  rmsd = cago_rmsd(go, go->x, NULL);
  fprintf(stderr, "warm-up t %g, RMSD %g, ncont %d/%d\n",
      t, rmsd, nc, go->ncont);

  for ( t = m->nstblk; t <= m->nsteps; t += m->nstblk ) {
    for ( it = 0; it < m->nstblk; it++ ) {
      /* a step of Metropolis alogrithm */
      tot += 1;
      acc += cagonc_metro(r, m->mcamp, 1/m->temp, &nc);

      wl_addi(wl, nc);
    }

    if ( wl_updatelnf(wl) ) {
      /* update the MC amplitude */
      update_mcamp(&m->mcamp, 0.5, &acc, &tot);
    }

    if ( fmod(t, m->nstrep) < 0.1 ) {
      double flatness = wl_getflatness(wl);
      rmsd = cago_rmsd(go, go->x, NULL);
      cago_writepos(go, go->x, go->v, m->fnpos);
      wl_save(wl, m->fnvnc);
      fprintf(stderr, "%g: acc %.2f%%, ep %g, rmsd %g, nc %d/%d, flatness %g%%, lnf %g\n",
          t, 100.0 * acc / tot, go->epot,
          rmsd, nc, go->ncont, 100.0 * flatness, wl->lnf);
      cago_rmcom(go, go->x, go->v);
    }
  }

  cago_close(go);
  cagonc_close(r);
  wl_close(wl);
  return 0;
}
