#ifndef GMXGOMODEL_H__
#define GMXGOMODEL_H__


typedef struct {
  double mflmin;
  double mflmax;
  double mfhmin;
  double mfhmax;
  char *fnvrmsd;

  double wl_lnf0;
  double wl_flatness;
  double wl_frac;
  double invt_c;

  int nstrep;
  char *fnrhis;
} gmxgomodel_t;



static gmxgomodel_t *gmxgomodel_open(void)
{
  gmxgomodel_t *m;

  snew(m, 1);

  m->mflmin = -1.0;
  m->mflmax = 1.0;
  m->mfhmin = -10.0;
  m->mfhmax = 0.0;
  m->fnvrmsd = "vrmsd.dat";

  m->wl_lnf0 = 1e-4;
  m->wl_flatness = 0.3;
  m->wl_frac = 0.5;
  m->invt_c = 1e4;

  m->nstrep = 10000;
  m->fnrhis = "rhis.dat";
  return m;
}



static void gmxgomodel_close(gmxgomodel_t *m)
{
  sfree(m);
}


#endif /* GMXGOMODEL_H__ */

