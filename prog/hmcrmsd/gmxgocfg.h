#ifndef GMXGOMODEL_H__
#define GMXGOMODEL_H__



/* input parameters */



#include "util.h"


enum {
  SEL_CA = 0,
  SEL_HEAVY,
  SEL_ALL,
  SEL_CAENDTOEND,
  SEL_COUNT
};

const char *seltype_names[] = {
  "CA",
  "heavy",
  "all",
  "CA-end-to-end",
  "count"
};


typedef struct {
  const char *fnpdb; /* reference PDB structure */

  int seltype; /* SEL_CA, SEL_HEAVY, SEL_ALL */

  int passive; /* only observe the RMSD distribution
                  do not change it */

  int dohmc; /* use HMC to reject unwanted configurations */

  int bias_mf; /* use mean force to estimate the bias potential */

  double mfmin;
  double mfmax;

  double mflmin;
  double mflmax;
  double mfhmin;
  double mfhmax;

  double minh;

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

  char *fnlog;

  int nstlog;
  int nstchat;
  int nstrep;

  int debug;
} gmxgocfg_t;



static void gmxgocfg_default(gmxgocfg_t *cfg)
{
  cfg->fnpdb = "1VII.pdb";

  cfg->seltype = SEL_ALL;

  cfg->passive = 0;

  cfg->dohmc = 1;

  cfg->bias_mf = 0;

  cfg->mfmin = -1000.0;
  cfg->mfmax = +1000.0;

  /* for rmsd < rmsdmin, a bias with dvdx < 0 helps
   * the RMSD increase and go back to the range */
  cfg->mflmin = -1000.0;
  cfg->mflmax = 0.0;

  /* for rmsd > rmsdmax, a bias with dvdx > 0 helps
   * the RMSD decrease and go back to the range */
  cfg->mfhmin = 10.0;
  cfg->mfhmax = 1000.0;

  cfg->minh = 0;

  cfg->fnvrmsd = "vrmsd.dat";

  /* RMSD range */
  cfg->rmsdmin = 0.10;
  cfg->rmsdmax = 1.00;
  cfg->rmsddel = 0.01;

  cfg->wl_lnf0 = 1e-4;
  cfg->wl_flatness = 0.3;
  cfg->wl_frac = 0.5;
  cfg->invt_c = 1;

  cfg->rhis_dx = 0.002;
  cfg->rhis_max = 4.0;
  cfg->fnrhis = "rhis.dat";

  cfg->fnlog = "rmsd.log";

  cfg->nstlog = 100;
  cfg->nstchat = 1000;
  cfg->nstrep = 10000;

  cfg->debug = 0;
}



static int array_select(const char *val, const char **names, int cnt)
{
  int i;

  for ( i = 0; i < cnt; i++ ) {
    if ( strcmpfuzzy(val, names[i]) == 0 ) {
      return i;
    }
  }

  return -1;
}


/* load settings from the configuration file `fn` */
static int gmxgocfg_load(gmxgocfg_t *cfg, const char *fn)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int inpar;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load configuration file %s\n", fn);
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
    for ( p++; *p && (isspace(*p) || *p == '=') ; )
      p++;
    val = p;

    if ( strcmpfuzzy(key, "pdb") == 0
      || strcmpfuzzy(key, "fnpdb") == 0 ) {
      cfg->fnpdb = strclone(val);
    } else if ( strcmpfuzzy(key, "seltype") == 0 ) {
      cfg->seltype = array_select(val, seltype_names, SEL_COUNT);
      //printf("cfg->seltype %d\n", cfg->seltype); getchar();
    } else if ( strcmpfuzzy(key, "passive") == 0 ) {
      if ( *val ) {
        cfg->passive = atoi(val);
      } else {
        cfg->passive = 1;
      }
    } else if ( strcmpfuzzy(key, "mfmin") == 0 ) {
      cfg->mfmin = atof(val);
    } else if ( strcmpfuzzy(key, "mfmax") == 0 ) {
      cfg->mfmax = atof(val);
    } else if ( strcmpfuzzy(key, "mflmin") == 0 ) {
      cfg->mflmin = atof(val);
    } else if ( strcmpfuzzy(key, "mflmax") == 0 ) {
      cfg->mflmax = atof(val);
    } else if ( strcmpfuzzy(key, "mfhmin") == 0 ) {
      cfg->mfhmin = atof(val);
    } else if ( strcmpfuzzy(key, "mfhmax") == 0 ) {
      cfg->mfhmax = atof(val);
    } else if ( strcmpfuzzy(key, "minh") == 0 ) {
      cfg->minh = atof(val);
    } else if ( strcmpfuzzy(key, "fnvrmsd") == 0 ) {
      cfg->fnvrmsd = strclone(val);
    } else if ( strcmpfuzzy(key, "rmsdmin") == 0 ) {
      cfg->rmsdmin = atof(val);
    } else if ( strcmpfuzzy(key, "rmsdmax") == 0 ) {
      cfg->rmsdmax = atof(val);
    } else if ( strcmpfuzzy(key, "rmsddel") == 0 ) {
      cfg->rmsddel = atof(val);
    } else if ( strcmpfuzzy(key, "wl_lnf0") == 0 ) {
      cfg->wl_lnf0 = atof(val);
    } else if ( strcmpfuzzy(key, "wl_flatness") == 0 ) {
      cfg->wl_flatness = atof(val);
    } else if ( strcmpfuzzy(key, "wl_frac") == 0 ) {
      cfg->wl_frac = atof(val);
    } else if ( strcmpfuzzy(key, "invt_c") == 0 ) {
      cfg->invt_c = atof(val);
    } else if ( strcmpfuzzy(key, "rhis_dx") == 0 ) {
      cfg->rhis_dx = atof(val);
    } else if ( strcmpfuzzy(key, "rhis_max") == 0 ) {
      cfg->rhis_max = atof(val);
    } else if ( strcmpfuzzy(key, "fnrhis") == 0 ) {
      cfg->fnrhis = strclone(val);
    } else if ( strcmpfuzzy(key, "fnlog") == 0 ) {
      cfg->fnlog = strclone(val);
    } else if ( strcmpfuzzy(key, "nstlog") == 0 ) {
      cfg->nstlog = atoi(val);
    } else if ( strcmpfuzzy(key, "nstchat") == 0 ) {
      cfg->nstchat = atoi(val);
    } else if ( strcmpfuzzy(key, "nstrep") == 0 ) {
      cfg->nstrep = atoi(val);
    } else if ( strcmpfuzzy(key, "hmc") == 0
             || strcmpfuzzy(key, "dohmc") == 0 ) {
      cfg->dohmc = atoi(val);
    } else if ( strcmpfuzzy(key, "biasmf") == 0
             || strcmpfuzzy(key, "bias_mf") == 0 ) {
      cfg->bias_mf = atoi(val);
    } else if ( strcmpfuzzy(key, "debug") == 0
             || strcmpfuzzy(key, "dbg") == 0 ) {
      cfg->debug = (*val) ? atoi(val) : 1;
    } else {
      fprintf(stderr, "Warning: unknown option %s = %s\n", key, val);
      getchar();
    }
  }

  fclose(fp);

  return 0;
}



#endif /* GMXGOMODEL_H__ */

