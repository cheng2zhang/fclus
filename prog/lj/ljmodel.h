#ifndef LJMODEL_H__
#define LJMODEL_H__





/* parameters and basic data types */





#include "util.h"





typedef struct {
  int n; /* number of particles */
  double rho;
  double temp;
  double beta;
  double rcdef; /* preferred radius cutoff */
  double rcls; /* clustering distance cutoff */
  double mcamp; /* MC amplitude */
  double mddt; /* MD time step */
  double thdt; /* thermostat step size */
  int nstblk;
  int nsthmc;
  double nsteps;
  double nstrep;
  int verbose;
  const char *prog;
  const char *fncfg;
  const char *fnpos;
  const char *fnvcls;
  const char *fnrep;
  double wl_lnf0;
  double wl_flatness;
  double wl_frac;
  double invt_c;
  int changeseed;
  int nvswaps;
} ljmodel_t;



/* set default values of the parameters */
__inline static void ljmodel_default(ljmodel_t *m)
{
  memset(m, 0, sizeof(*m));
  m->n = 108;
#if D == 2
  m->rho = 0.35;
#else
  m->rho = 0.25;
#endif
  m->temp = 2;
  m->beta = 0.5;
  m->rcdef = 1e9;
  m->rcls = 1.6;
  m->mcamp = 0.2;
  m->mddt = 0.002;
  m->thdt = 0.1;
  m->nstblk = 10;
  m->nsthmc = 1;
  m->nsteps = 1e10;
  m->nstrep = 10000;
  m->verbose = 0;
  m->prog = NULL;
  m->fncfg = NULL;
  m->fnpos = "lj.pos";
  m->fnvcls = "vcls.dat";
  m->fnrep = NULL;
  m->wl_lnf0 = 4e-4;
  m->wl_flatness = 0.3;
  m->wl_frac = 0.5;
  m->invt_c = 1.0;
  m->changeseed = 0;
  m->nvswaps = 1;
}



/* print help message and die */
__inline static void ljmodel_help(const ljmodel_t *m)
{
  fprintf(stderr, "Lennard-Jones fluid\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] [input.cfg]\n\n", m->prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -n:            set the number of particles, default %d\n", m->n);
  fprintf(stderr, "  -r, --rho=:    set the (maximal) density, default %g\n", m->rho);
  fprintf(stderr, "  -T:            set the temperautre, default %g\n", m->temp);
  fprintf(stderr, "  --rc=:         set the default cutoff, default %g\n", m->rcdef);
  fprintf(stderr, "  --rcls=:       set the clustering distance cutoff, default %g\n", m->rcls);
  fprintf(stderr, "  -A, --mcamp=:  set the MC amplitude, default %g\n", m->mcamp);
  fprintf(stderr, "  --dt=:         set the MD time step, default %g\n", m->mddt);
  fprintf(stderr, "  --thdt=:       set the thermostat time step, default %g\n", m->thdt);
  fprintf(stderr, "  --nstblk=:     set the number of steps in each block, default %d\n", m->nstblk);
  fprintf(stderr, "  --nsthmc=:     set the number of steps in each hybrid MC step, default %d\n", m->nsthmc);
  fprintf(stderr, "  --nsteps=:     set the number of steps, default %g\n", m->nsteps);
  fprintf(stderr, "  --nstrep=:     set the number of steps for reporting, default %g\n", m->nstrep);
  fprintf(stderr, "  --cfg=:        set the configuration file, default: %s\n", m->fncfg);
  fprintf(stderr, "  --pos=:        set the output coordinates file, default: %s\n", m->fnpos);
  fprintf(stderr, "  --vcls=:       set the output file for the bias potential, default: %s\n", m->fnvcls);
  fprintf(stderr, "  --rep=:        set the output report file, default: %s\n", m->fnrep);
  fprintf(stderr, "  --lnf0=:       set the initial lnf, default: %g\n", m->wl_lnf0);
  fprintf(stderr, "  --flatness=:   set the histogram flatness for lnf, default: %g\n", m->wl_flatness);
  fprintf(stderr, "  --frac=:       set the reduction factor for lnf, default: %g\n", m->wl_frac);
  fprintf(stderr, "  --invt_c=:     set the constant c in lnf = c/t, default: %g\n", m->invt_c);
  fprintf(stderr, "  --chseed:      randomly change the seed particle, default: %d\n", m->changeseed);
  fprintf(stderr, "  --nvswaps:     set the number of velocity swaps, default: %d\n", m->nvswaps);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* compute dependent variables */
__inline static void ljmodel_compute(ljmodel_t *m)
{
  m->beta = 1 / m->temp;
}



/* load settings from the configuration file `fn' */
__inline static int ljmodel_load(ljmodel_t *m, const char *fn)
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
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; ) p++;
    val = p;
    for ( ; *p; p++ ) *p = (char) tolower(*p);

    if ( strcmpfuzzy(key, "n") ) {
      m->n = atoi(val);
    } else if ( strcmpfuzzy(key, "rho") ) {
      m->rho = atof(val);
    } else if ( strcmpfuzzy(key, "T") == 0
             || strcmpfuzzy(key, "temp") == 0 ) {
      m->temp = atof(val);
    } else if ( strcmpfuzzy(key, "rc") == 0 ) {
      m->rcdef = atof(val);
    } else if ( strcmpfuzzy(key, "rcls") == 0 ) {
      m->rcls = atof(val);
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
      m->fnpos = strclone(val);
    } else if ( strcmpfuzzy(key, "fnvcls") == 0 ) {
      m->fnvcls = strclone(val);
    } else if ( strcmpfuzzy(key, "fnrep") == 0 ) {
      m->fnrep = strclone(val);
    } else if ( strcmpfuzzy(key, "wl_lnf0") == 0 ) {
      m->wl_lnf0 = atof(val);
    } else if ( strcmpfuzzy(key, "wl_flatness") == 0 ) {
      m->wl_flatness = atof(val);
    } else if ( strcmpfuzzy(key, "wl_frac") == 0 ) {
      m->wl_frac = atof(val);
    } else if ( strcmpfuzzy(key, "invt_c") == 0 ) {
      m->invt_c = atof(val);
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
  ljmodel_compute(m);

  return 0;
}



/* handle command line arguments */
__inline static void ljmodel_doargs(ljmodel_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  /* reset */
  ljmodel_default(m);

  /* set the program name */
  m->prog = argv[0];

  for ( i = 1; i < argc; i++ ) {
    /* it's an argument */
    if ( argv[i][0] != '-' ) {
      m->fncfg = argv[i];
      if ( ljmodel_load(m, m->fncfg) != 0 ) {
        fprintf(stderr, "failed to load %s\n", m->fncfg);
        ljmodel_help(m);
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

      if ( strcmp(p, "rho") == 0 ) {
        m->rho = atof(q);
      } else if ( strncmp(p, "temp", 4) == 0 ) {
        m->temp = atof(q);
      } else if ( strcmp(p, "rc") == 0 ) {
        m->rcdef = atof(q);
      } else if ( strcmp(p, "rcls") == 0 ) {
        m->rcls = atof(q);
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
      } else if ( strcmp(p, "vcls") == 0 ) {
        m->fnvcls = q;
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
      } else if ( strcmp(p, "chseed") == 0 ) {
        m->changeseed = 1;
      } else if ( strcmp(p, "nvswaps") == 0 ) {
        m->nvswaps = atoi(q);
      } else if ( strcmp(p, "help") == 0 ) {
        ljmodel_help(m);
      } else {
        fprintf(stderr, "unknown option %s\n", argv[i]);
        ljmodel_help(m);
      }

      continue;
    }

    /* it is an option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( strchr("nrTA", ch) != NULL ) {
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
          ljmodel_help(m);
        }

        if ( ch == 'n' ) { /* override the number of particles */
          m->n = atoi(q);
        } else if ( ch == 'r' ) { /* override the density */
          m->rho = atof(q);
        } else if ( ch == 'T' ) { /* override the temperature */
          m->temp = atof(q);
        } else if ( ch == 'A' ) { /* override the MC amplitude */
          m->mcamp = atof(q);
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == 'v' ) {
        m->verbose++;
      } else if ( ch == 'h' ) {
        ljmodel_help(m);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        ljmodel_help(m);
      }
    }
  }

  ljmodel_compute(m);
}





#endif /* LJMODEL_H__ */
