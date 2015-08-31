#ifndef GMXGOMODEL_H__
#define GMXGOMODEL_H__



#include "util.h"



typedef struct {
  const char *fnpdb; /* reference PDB structure */

  double mfmin;
  double mfmax;

  double mflmin;
  double mflmax;
  double mfhmin;
  double mfhmax;
  char *fnvrmsd;

  double rmsdmin; /* minimal RMSD in angstroms */
  double rmsdmax;
  double rmsddel;

  double wl_lnf0;
  double wl_flatness;
  double wl_frac;
  double invt_c;

  double rhis_dx;
  double rhis_max;
  char *fnrhis;

  int nstrep;
} gmxgomodel_t;



static void gmxgomodel_default(gmxgomodel_t *m)
{
  m->fnpdb = "1VII.pdb";
  m->mfmin = -4.0;
  m->mfmax = 4.0;

  m->mflmin = -1.0;
  m->mflmax = 0.0;
  m->mfhmin = 1.0;
  m->mfhmax = 3.0;
  m->fnvrmsd = "vrmsd.dat";

  m->rmsdmin = 0.08;
  m->rmsdmax = 0.22;
  m->rmsddel = 0.01;

  m->wl_lnf0 = 1e-6;
  m->wl_flatness = 0.3;
  m->wl_frac = 0.5;
  m->invt_c = 1;

  m->rhis_dx = 0.002;
  m->rhis_max = 4.0;
  m->fnrhis = "rhis.dat";

  m->nstrep = 10000;
}



/* load settings from the configuration file `fn` */
static int gmxgomodel_load(gmxgomodel_t *m, const char *fn)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int inpar;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    strstrip(buf); /* remove trailing spaces */
    if ( buf[0] == '\0' || buf[0] == '#'
      || buf[0] == ';' )
      continue;

    /* parse the line to a key and a value */
    /* find the end of the key
     * note the key may contain an index, e.g., in
     *   arr(1) = val;
     * the key is `arr(1)`, the value is `val`. */
    inpar = 0; /* within paretheses (...) or brackets [...] */
    for ( p = buf; *p; p++ ) {
      if ( (inpar == 0 && isspace(*p)) || *p == '=' ) {
        *p = '\0';
        break;
      }
      if ( *p == '(' || *p == '[' )  {
        inpar++; /* entering a paretheses region */
      } else if ( inpar && ( *p == ')' || *p == ']' ) ) {
        if ( --inpar <= 0 ) inpar = 0;
      }

      /* convert the key to lower case */
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; ) p++;
    val = p;
    for ( ; *p; p++ ) *p = (char) tolower(*p);

    if ( strcmpfuzzy(key, "pdb") == 0
      || strcmpfuzzy(key, "fnpdb") == 0 ) {
      m->fnpdb = val;
    } else if ( strcmpfuzzy(key, "mfmin") == 0 ) {
      m->mfmin = atof(val);
    } else if ( strcmpfuzzy(key, "mfmax") == 0 ) {
      m->mfmax = atof(val);
    } else if ( strcmpfuzzy(key, "mflmin") == 0 ) {
      m->mflmin = atof(val);
    } else if ( strcmpfuzzy(key, "mflmax") == 0 ) {
      m->mflmax = atof(val);
    } else if ( strcmpfuzzy(key, "mfhmin") == 0 ) {
      m->mfhmin = atof(val);
    } else if ( strcmpfuzzy(key, "mfhmax") == 0 ) {
      m->mfhmax = atof(val);
    } else if ( strcmpfuzzy(key, "fnvrmsd") == 0 ) {
      m->fnvrmsd = val;
    } else if ( strcmpfuzzy(key, "rmsdmin") == 0 ) {
      m->rmsdmin = atof(val);
    } else if ( strcmpfuzzy(key, "rmsdmax") == 0 ) {
      m->rmsdmax = atof(val);
    } else if ( strcmpfuzzy(key, "rmsddel") == 0 ) {
      m->rmsddel = atof(val);
    } else if ( strcmpfuzzy(key, "wl_lnf0") == 0 ) {
      m->wl_lnf0 = atof(val);
    } else if ( strcmpfuzzy(key, "wl_flatness") == 0 ) {
      m->wl_flatness = atof(val);
    } else if ( strcmpfuzzy(key, "wl_frac") == 0 ) {
      m->wl_frac = atof(val);
    } else if ( strcmpfuzzy(key, "invt_c") == 0 ) {
      m->invt_c = atof(val);
    } else if ( strcmpfuzzy(key, "rhis_dx") == 0 ) {
      m->rhis_dx = atof(val);
    } else if ( strcmpfuzzy(key, "rhis_max") == 0 ) {
      m->rhis_max = atof(val);
    } else if ( strcmpfuzzy(key, "nstrep") == 0 ) {
      m->nstrep = atoi(val);
    } else if ( strcmpfuzzy(key, "fnrhis") == 0 ) {
      m->fnrhis = val;
    } else {
      fprintf(stderr, "Warning: unknown option %s = %s\n", key, val);
      getchar();
    }
  }

  fclose(fp);

  return 0;
}

#endif /* GMXGOMODEL_H__ */

