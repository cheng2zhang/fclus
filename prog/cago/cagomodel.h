#ifndef CAGOMODEL_H__
#define CAGOMODEL_H__





/* parameters and basic data types */





#include "util.h"





typedef struct {
  const char *fnpdb;
  double kb;
  double ka;
  double kd1;
  double kd3;
  double nbe;
  double nbc;
  double temp;
  double beta;
  double rc; /* contact distance cutoff */
  double mcamp; /* MC amplitude */
  double mddt; /* MD time step */
  double thdt; /* thermostat step size */
  int nstblk;
  int nsthmc;
  int implicithmc; /* implicit HMC */
  double nsteps;
  double nstrep;
  int verbose;
  const char *prog;
  const char *fncfg;
  const char *fnpos;
  const char *fnvnc;
  const char *fnvrmsd;
  const char *fnrep;
  double fncmin; /* minimal fraction of number of contacts */
  double fncmax; /* maximal fraction of number of contacts */
  double rmsdmin;
  double rmsdmax;
  double rmsddel;
  double wl_lnf0;
  double wl_flatness;
  double wl_frac;
  double invt_c;
  double mflmin; /* minimal mean force for r < rmin */
  double mflmax; /* maximal mean force for r < rmin */
  double mfhmin; /* minimal mean force for r > rmax */
  double mfhmax; /* maximal mean force for r > rmax */
  int changeseed;
  int nvswaps;
} cagomodel_t;



/* set default values of the parameters */
__inline static void cagomodel_default(cagomodel_t *m)
{
  memset(m, 0, sizeof(*m));
  m->fnpdb = "pdb/1L2Y.pdb";
  m->kb = 200.0;
  m->ka = 40.0;
  m->kd1 = 1.0;
  m->kd3 = 0.5;
  m->nbe = 1.0;
  m->nbc = 4.0;
  m->rc = 6.0;
  m->temp = 1.02;
  m->beta = 1/m->temp;
  m->mcamp = 0.2;
  m->mddt = 0.002;
  m->thdt = 0.1;
  m->nstblk = 10;
  m->nsthmc = 1;
  m->implicithmc = 0;
  m->nsteps = 1e10;
  m->nstrep = 1000000;
  m->verbose = 0;
  m->prog = NULL;
  m->fncfg = NULL;
  m->fnpos = "lj.pos";
  m->fnvnc = "vnc.dat";
  m->fnvrmsd = "vrmsd.dat";
  m->fnrep = NULL;
  m->fncmin = 0.0;
  m->fncmax = 1.0;
  m->rmsdmin = 1.0;
  m->rmsdmax = 6.0;
  m->rmsddel = 0.05;
  m->wl_lnf0 = 1e-4;
  m->wl_flatness = 0.3;
  m->wl_frac = 0.5;
  m->invt_c = 1.0;
  m->mflmin = 0.1;
  m->mflmax = 1e10;
  m->mfhmin = -1e10;
  m->mfhmax = -0.1;
  m->changeseed = 0;
  m->nvswaps = 1;
}



/* print help message and die */
__inline static void cagomodel_help(const cagomodel_t *m)
{
  fprintf(stderr, "Alpha-carbon Go model\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] [input.cfg]\n\n", m->prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i, --pdb=     set the input protein databank PDB, default %s\n", m->fnpdb);
  fprintf(stderr, "  -T:            set the temperautre, default %g\n", m->temp);
  fprintf(stderr, "  --rc=:         set the contact distance cutoff, default %g\n", m->rc);
  fprintf(stderr, "  -A, --mcamp=:  set the MC amplitude, default %g\n", m->mcamp);
  fprintf(stderr, "  --dt=:         set the MD time step, default %g\n", m->mddt);
  fprintf(stderr, "  --thdt=:       set the thermostat time step, default %g\n", m->thdt);
  fprintf(stderr, "  --nstblk=:     set the number of steps in each block, default %d\n", m->nstblk);
  fprintf(stderr, "  --nsthmc=:     set the number of steps in each hybrid MC step, default %d\n", m->nsthmc);
  fprintf(stderr, "  --ihmc:        use implicit HMC, default %d\n", m->implicithmc);
  fprintf(stderr, "  --nsteps=:     set the number of steps, default %g\n", m->nsteps);
  fprintf(stderr, "  --nstrep=:     set the number of steps for reporting, default %g\n", m->nstrep);
  fprintf(stderr, "  --cfg=:        set the configuration file, default: %s\n", m->fncfg);
  fprintf(stderr, "  --pos=:        set the output coordinates file, default: %s\n", m->fnpos);
  fprintf(stderr, "  --vnc=:        set the output file for the NC bias potential, default: %s\n", m->fnvnc);
  fprintf(stderr, "  --vrmsd=:      set the output file for the RMSD bias potential, default: %s\n", m->fnvrmsd);
  fprintf(stderr, "  --rep=:        set the output report file, default: %s\n", m->fnrep);
  fprintf(stderr, "  --rmin=:       set the minimal RMSD, default %g\n", m->rmsdmin);
  fprintf(stderr, "  --rmax=:       set the maximal RMSD, default %g\n", m->rmsdmax);
  fprintf(stderr, "  --lnf0=:       set the initial lnf, default: %g\n", m->wl_lnf0);
  fprintf(stderr, "  --flatness=:   set the histogram flatness for lnf, default: %g\n", m->wl_flatness);
  fprintf(stderr, "  --frac=:       set the reduction factor for lnf, default: %g\n", m->wl_frac);
  fprintf(stderr, "  --invt_c=:     set the constant c in lnf = c/t, default: %g\n", m->invt_c);
  fprintf(stderr, "  --mflmin=      set the minimal mean-force for r < rmin, default %g\n", m->mflmin);
  fprintf(stderr, "  --mflmax=      set the maximal mean-force for r < rmin, default %g\n", m->mflmax);
  fprintf(stderr, "  --mfhmin=      set the minimal mean-force for r > rmax, default %g\n", m->mfhmin);
  fprintf(stderr, "  --mfhmax=      set the maximal mean-force for r > rmax, default %g\n", m->mfhmax);
  fprintf(stderr, "  --chseed:      randomly change the seed particle, default: %d\n", m->changeseed);
  fprintf(stderr, "  --nvswaps:     set the number of velocity swaps, default: %d\n", m->nvswaps);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* compute dependent variables */
__inline static void cagomodel_compute(cagomodel_t *m)
{
  m->beta = 1 / m->temp;
}



/* load settings from the configuration file `fn' */
__inline static int cagomodel_load(cagomodel_t *m, const char *fn)
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
    if ( buf[0] == '\0' || buf[0] == '#' ) continue;

    /* parse the line to a key and a value */
    /* find the end of the key */
    inpar = 0; /* within (...) */
    for ( p = buf; *p; p++ ) {
      if ( (!inpar && isspace(*p)) || *p == '=' ) {
        *p = '\0'; /* end the key part */
        break;
      }
      if ( !inpar && (*p == '(' || *p == '[') )
        inpar = 1; /* enter a parentheses block */
      else if ( inpar && (*p == ')' || *p == ']') )
        inpar = 0; /* leave a parentheses block */

      /* convert the key to lower case */
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; ) p++;
    val = p;
    for ( ; *p; p++ ) *p = (char) tolower(*p);

    if ( strcmpfuzzy(key, "pdb") == 0 ) {
      m->fnpdb = val;
    } else if ( strcmpfuzzy(key, "kb") == 0 ) {
      m->kb = atof(val);
    } else if ( strcmpfuzzy(key, "ka") == 0 ) {
      m->ka = atof(val);
    } else if ( strcmpfuzzy(key, "kd1") == 0 ) {
      m->kd1 = atof(val);
    } else if ( strcmpfuzzy(key, "kd3") == 0 ) {
      m->kd3 = atof(val);
    } else if ( strcmpfuzzy(key, "nbe") == 0 ) {
      m->nbe = atof(val);
    } else if ( strcmpfuzzy(key, "nbc") == 0 ) {
      m->nbc = atof(val);
    } else if ( strcmpfuzzy(key, "ihmc") == 0 ) {
      if ( strcmpfuzzy(key, "true") == 0 ) {
        m->implicithmc = 1;
      } else if ( strcmpfuzzy(key, "false") == 0 ) {
        m->implicithmc = 0;
      } else { /* numerical value, 0 or 1 */
        m->implicithmc = atoi(key);
      }
    } else if ( strcmpfuzzy(key, "rc") == 0 ) {
      m->rc = atof(val);
    } else if ( strcmpfuzzy(key, "T") == 0
             || strcmpfuzzy(key, "temp") == 0 ) {
      m->temp = atof(val);
    } else if ( strcmpfuzzy(key, "rmin") == 0 ) {
      m->rmsdmin = atof(val);
    } else if ( strcmpfuzzy(key, "rmax") == 0 ) {
      m->rmsdmax = atof(val);
    } else if ( strcmpfuzzy(key, "mcamp") == 0 ) {
      m->mcamp = atof(val);
    } else if ( strcmpfuzzy(key, "mddt") == 0 ) {
      m->mddt = atof(val);
    } else if ( strcmpfuzzy(key, "thdt") == 0
             || strcmpfuzzy(key, "thermdt") == 0 ) {
      m->thdt = atof(val);
    } else if ( strcmpfuzzy(key, "nstblk") == 0 ) {
      m->nstblk = atoi(val);
    } else if ( strcmpfuzzy(key, "nsthmc") == 0 ) {
      m->nsthmc = atoi(val);
    } else if ( strcmpfuzzy(key, "nsteps") == 0 ) {
      m->nsteps = atof(val);
    } else if ( strcmpfuzzy(key, "nstrep") == 0 ) {
      m->nstrep = atof(val);
    } else if ( strcmpfuzzy(key, "fnpos") == 0 ) {
      m->fnpos = val;
    } else if ( strcmpfuzzy(key, "fnvnc") == 0 ) {
      m->fnvnc = val;
    } else if ( strcmpfuzzy(key, "fnvrmsd") == 0 ) {
      m->fnvrmsd = val;
    } else if ( strcmpfuzzy(key, "fnrep") == 0 ) {
      m->fnrep = val;
    } else if ( strcmpfuzzy(key, "wl_lnf0") == 0 ) {
      m->wl_lnf0 = atof(val);
    } else if ( strcmpfuzzy(key, "wl_flatness") == 0 ) {
      m->wl_flatness = atof(val);
    } else if ( strcmpfuzzy(key, "wl_frac") == 0 ) {
      m->wl_frac = atof(val);
    } else if ( strcmpfuzzy(key, "invt_c") == 0 ) {
      m->invt_c = atof(val);
    } else if ( strncmpfuzzy(key, "mflmin", 6) == 0 ) {
      m->mflmin = atof(val);
    } else if ( strncmpfuzzy(key, "mflmax", 6) == 0 ) {
      m->mflmax = atof(val);
    } else if ( strncmpfuzzy(key, "mfhmin", 6) == 0 ) {
      m->mfhmin = atof(val);
    } else if ( strncmpfuzzy(key, "mfhmax", 6) == 0 ) {
      m->mfhmax = atof(val);
    } else if ( strcmpfuzzy(key, "changeseed") == 0
            ||  strcmpfuzzy(key, "chseed") == 0 ) {
      m->changeseed = ( strcmpfuzzy(val, "true") == 0 );
    } else if ( strcmpfuzzy(key, "nvswaps") == 0 ) {
      m->nvswaps = atoi(val);
    } else {
      fprintf(stderr, "Warning: unknown option %s = %s\n", key, val);
      getchar();
    }
  }

  fclose(fp);
  cagomodel_compute(m);

  return 0;
}



/* handle command line arguments */
__inline static void cagomodel_doargs(cagomodel_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  /* reset */
  cagomodel_default(m);

  /* set the program name */
  m->prog = argv[0];

  for ( i = 1; i < argc; i++ ) {
    /* it's an argument */
    if ( argv[i][0] != '-' ) {
      m->fncfg = argv[i];
      if ( cagomodel_load(m, m->fncfg) != 0 ) {
        fprintf(stderr, "failed to load %s\n", m->fncfg);
        cagomodel_help(m);
      }
      continue;
    }

    /* long argument, like --help */
    if ( argv[i][1] == '-' ) {
      /* try to parse the argment
         e.g., `--prog=aaa' is parsed to `--prog' and `aaa' */
      p = argv[i] + 2;
      /* let q point to the argument of the option */
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }

      if ( strcmpfuzzy(p, "pdb") == 0 ) {
        m->fnpdb = q;
      } else if ( strcmpfuzzy(p, "ihmc") == 0 ) {
        m->implicithmc = 1;
      } else if ( strncmp(p, "temp", 4) == 0 ) {
        m->temp = atof(q);
      } else if ( strcmp(p, "rc") == 0 ) {
        m->rc = atof(q);
      } else if ( strcmpfuzzy(p, "rmin") == 0 ) {
        m->rmsdmin = atof(q);
      } else if ( strcmpfuzzy(p, "rmax") == 0 ) {
        m->rmsdmax = atof(q);
      } else if ( strcmp(p, "mcamp") == 0 ) {
        m->mcamp = atof(q);
      } else if ( strcmp(p, "mddt") == 0 ) {
        m->mddt = atof(q);
      } else if ( strcmp(p, "thdt") == 0 ) {
        m->thdt = atof(q);
      } else if ( strcmp(p, "nstblk") == 0 ) {
        m->nstblk = atoi(q);
      } else if ( strcmp(p, "nsthmc") == 0 ) {
        m->nsthmc = atoi(q);
      } else if ( strcmp(p, "nsteps") == 0 ) {
        m->nsteps = atof(q);
      } else if ( strcmp(p, "nstrep") == 0 ) {
        m->nstrep = atof(q);
      } else if ( strcmp(p, "cfg") == 0 ) {
        m->fncfg = q;
      } else if ( strcmp(p, "pos") == 0 ) {
        m->fnpos = q;
      } else if ( strcmp(p, "vnc") == 0 ) {
        m->fnvnc = q;
      } else if ( strcmp(p, "vrmsd") == 0 ) {
        m->fnvrmsd = q;
      } else if ( strcmp(p, "rep") == 0 ) {
        m->fnrep = q;
      } else if ( strcmp(p, "lnf0") == 0 ) {
        m->wl_lnf0 = atof(q);
      } else if ( strcmp(p, "flatness") == 0 ) {
        m->wl_flatness = atof(q);
      } else if ( strcmp(p, "frac") == 0 ) {
        m->wl_frac = atof(q);
      } else if ( strcmp(p, "invt_c") == 0 ) {
        m->invt_c = atof(q);
      } else if ( strncmp(p, "mflmin", 6) == 0 ) {
        m->mflmin = atof(q);
      } else if ( strncmp(p, "mflmax", 6) == 0 ) {
        m->mflmax = atof(q);
      } else if ( strncmp(p, "mfhmin", 6) == 0 ) {
        m->mfhmin = atof(q);
      } else if ( strncmp(p, "mfhmax", 6) == 0 ) {
        m->mfhmax = atof(q);
      } else if ( strcmp(p, "chseed") == 0 ) {
        m->changeseed = 1;
      } else if ( strcmp(p, "nvswaps") == 0 ) {
        m->nvswaps = atoi(q);
      } else if ( strcmp(p, "help") == 0 ) {
        cagomodel_help(m);
      } else {
        fprintf(stderr, "unknown option %s\n", argv[i]);
        cagomodel_help(m);
      }

      continue;
    }

    /* it is an option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( strchr("iTA", ch) != NULL ) {
        /* handle options that require an argument */
        q = p = argv[i] + j + 1;
        if ( *p != '\0' ) {
          /* the argument follows immediately after the option
           * e.g., -oa.dat */
          q = p;
        } else if ( ++i < argc ) {
          /* the option and argument are separated by a space
           * then the argument belongs to the next argv[] element,
           * hence ++i
           * e.g., -o a.dat */
          q = argv[i];
        } else {
          fprintf(stderr, "-%c requires an argument!\n", ch);
          cagomodel_help(m);
        }

        if ( ch == 'i' ) { /* override the input PDB */
          m->fnpdb = q;
        } else if ( ch == 'T' ) { /* override the temperature */
          m->temp = atof(q);
        } else if ( ch == 'A' ) { /* override the MC amplitude */
          m->mcamp = atof(q);
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == 'v' ) {
        m->verbose++;
      } else if ( ch == 'h' ) {
        cagomodel_help(m);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        cagomodel_help(m);
      }
    }
  }

  cagomodel_compute(m);
}





#endif /* CAGOMODEL_H__ */
