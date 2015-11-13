#include "../util.h"
#include "../mtrand.h"



#define ALIAS_MAX 8

/* methods of computing the acceptance probability */
enum {
  HMCACC_FOUR_POTENTIAL = 0,
  HMCACC_TOTAL_ENERGY,
  HMCACC_WORK,
  HMCACC_COUNT
};

const char *hmcacc_names[][ALIAS_MAX] = {
  {"4", "4-pot", "4-potential", "four", "four-potential", ""},
  {"e", "tote", "totene", "total-energy", ""},
  {"w", "work", ""}
};



/* MD integrators */
enum {
  INTEGRATOR_VV = 0,
  INTEGRATOR_LEAPFROG,
  INTEGRATOR_COUNT
};

const char *integrator_names[][ALIAS_MAX] = {
  {"v", "vv", "Verlet", "velocity-Verlet", ""},
  {"l", "leapfrog", ""},
};


typedef struct {
  int dohmc;
  double nsteps;
  double tp;
  double dt;
  double thdt;
  double xpart;   /* the location of the discontinuity of the double well */
  double xmax;
  double dx;
  int nstblk;     /* block of MD steps between HMC trials */
  int hmcacc;     /* method of computing the acceptance probability */
  int integrator; /* MD integrator */
  int vflip;
  const char *prog;
  const char *fnhis;
} model_t;


static void model_default(model_t *m)
{
  m->dohmc = 0;
  m->nsteps = 10000000;
  m->tp = 1.0;
  m->dt = 0.005;
  m->thdt = 0.1;
  m->xpart = -0.5;
  m->xmax = 10;
  m->dx = 0.05;
  m->nstblk = 1;
  m->hmcacc = HMCACC_FOUR_POTENTIAL;
  m->integrator = INTEGRATOR_VV;
  m->vflip = 1;
  m->fnhis = "hist.dat";
}



static void model_help(const model_t *m)
{
  fprintf(stderr, "%s [OPTIONS]\n", m->prog);
  fprintf(stderr, "  --hmc:               do hmc\n");
  fprintf(stderr, "  --nsteps=%-8.0f:   set the number of steps\n", m->nsteps);
  fprintf(stderr, "  --tp=%-6.3f:         set the temperature\n", m->tp);
  fprintf(stderr, "  --dt=%-6.3f:         set the time step\n", m->dt);
  fprintf(stderr, "  --thdt=%-6.3f:       set the thermostat time step\n", m->thdt);
  fprintf(stderr, "  --xc=%-7.3f:        set the partition position\n", m->xpart);
  fprintf(stderr, "  --block=%-4d:        set the number of MD steps between HMC trials\n", m->nstblk);
  fprintf(stderr, "  --hmcacc=%s:          set the method of computing the acceptance probability, `4` (four-potential),`e` (total energy), `w` (work)\n", hmcacc_names[0][0]);
  fprintf(stderr, "  --novi=%d:            do not flip velocity\n", m->vflip);
  exit(1);
}



/* match the name of option */
static int select_opt(const char *val, const char *names[][ALIAS_MAX], int cnt)
{
  int i, j;

  for ( i = 0; i < cnt; i++ ) {
    for ( j = 0; j < ALIAS_MAX; j++ ) {
      if ( names[i][j][0] == '\0' ) {
        break;
      }

      /* compare `val` with the `j`th alias of name `i` */
      if ( strcmpfuzzy(val, names[i][j]) == 0 ) {
        return i;
      }
    }
  }

  fprintf(stderr, "cannot match the option for [%s]\n", val);

  return 0; /* default */
}



static int model_doargs(model_t *m, int argc, char **argv)
{
  int i;
  char *p, *q;

  m->prog = argv[0];
  for ( i = 1; i < argc; i++ ) {
    if ( argv[i][0] != '-' ) {
      continue;
    }

    if ( argv[i][1] == '-' ) { /* long option */
      p = argv[i] + 2;
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }

      if ( strcmp(p, "hmc") == 0 ) {
        m->dohmc = 1;
      } else if ( strcmp(p, "nsteps") == 0 ) {
        m->nsteps = atof(q);
      } else if ( strcmp(p, "tp") == 0 ) {
        m->tp = atof(q);
      } else if ( strcmp(p, "dt") == 0 ) {
        m->dt = atof(q);
      } else if ( strcmp(p, "thdt") == 0 ) {
        m->thdt = atof(q);
      } else if ( strcmp(p, "xpart") == 0 ) {
        m->xpart = atof(q);
      } else if ( strcmp(p, "nstblk") == 0
               || strcmp(p, "block") == 0 ) {
        m->nstblk = atoi(q);
      } else if ( strcmp(p, "hmcacc") == 0 ) {
        m->hmcacc = select_opt(q, hmcacc_names, HMCACC_COUNT);
      } else if ( strcmp(p, "int") == 0
               || strcmp(p, "integrator") == 0 ) {
        m->integrator = select_opt(q, integrator_names, INTEGRATOR_COUNT);
      } else if ( strcmp(p, "novi") == 0
               || strcmp(p, "novflip") == 0 ) {
        m->vflip = 0;
      } else if ( strcmp(p, "help") == 0 ) {
        model_help(m);
      }
      continue;
    } else {  /* short options */
      for ( p = argv[i] + 1; *p; p++ ) {
        if ( *p == 'h' ) {
          model_help(m);
        }
      }
    }
  }
  return 0;
}



static double force(model_t *m, double x,
    double *f, double *xi)
{
  double dx;

  *xi = ( x > m->xpart ? 1.0 : -1.0 );
  dx = x - *xi;
  *f = -dx;
  return 0.5 * dx * dx;
}



/* velocity rescaling thermostat */
static double vrescale(double thdt, double *v, double tp)
{
  double ek1, ek2, s, c, r, r2;

  /* normal velocity rescaling */
  ek1 = (*v) * (*v) / 2;

  /* approximate algorithm of computing ek2
   * only valid for a small thdt <= 0.001 */
  c = ( thdt < 700 ) ? exp(-thdt) : 0;
  r = randgaus();
  r2 = 0; /* randchisqr(dof - 1); */
  ek2 = c * ek1 + (1 - c) * (r2 + r*r) * .5 * tp
      + r * sqrt(c * (1 - c) * 2 * tp * ek1);

  if ( ek2 < 1e-30 ) ek2 = 1e-30;
  s = sqrt(ek2/ek1);
  *v *= s;

  return ek2;
}



static int savehis(model_t *m, double *his, int cnt)
{
  int i;
  FILE *fp;

  if ( (fp = fopen(m->fnhis, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", m->fnhis);
    return -1;
  }

  for ( i = 0; i < cnt; i++ ) {
    fprintf(fp, "%g %g\n", -m->xmax + (i + 0.5) * m->dx, his[i]);
  }

  fclose(fp);
  return 0;
}



static int domd(model_t *m)
{
  int t, nhis, i;
  double x = 0.01, v = 0, f = 0, dt = m->dt;
  double ep, xi;
  double *hist;

  nhis = (int) ( m->xmax * 2 / m->dx ) + 1;
  xnew(hist, nhis);

  v = sqrt(m->tp);
  ep = force(m, x, &f, &xi);
  for ( t = 0; t < m->nsteps; t++ ) {
    double x0, xi0, v0, f0, ep0;
    double dee = 0; /* change of the total energy */
    double ep1, ek1;

    x0 = x; /* save the old x */
    xi0 = xi; /* save the old reference point */
    v0 = v;
    f0 = f;
    ep0 = ep;

    if ( m->integrator == INTEGRATOR_VV || 1 ) {
      /* velocity-verlet */
      for ( i = 0; i < m->nstblk; i++ ) {
        vrescale(m->thdt * 0.5, &v, m->tp);
        ep1 = ep;
        ek1 = 0.5 * v * v;
        v += f * 0.5 * dt / m->nstblk;
        x += v * dt / m->nstblk;
        ep = force(m, x, &f, &xi);
        v += f * 0.5 * dt / m->nstblk;
        dee += ep - ep1 + 0.5 * v * v - ek1;
        vrescale(m->thdt * 0.5, &v, m->tp);
      }
    } else {
      /* leap-frog */
      for ( i = 0; i < m->nstblk; i++ ) {
        ep1 = ep;
        x += v * dt / m->nstblk;
        ep = force(m, x, &f, &xi);
        dee += ep - ep1;
        vrescale(m->thdt, &v, m->tp);
        ek1 = 0.5 * v * v;
        v += f * dt / m->nstblk;
        dee += 0.5 * v * v - ek1;
      }
    }

    /* use HMC to effect the change
     * due to the reference coordinate */
    if ( m->dohmc ) {
      double de, r;
      double ep2, ep3;
      int acc = 1;

      if ( m->hmcacc == HMCACC_WORK ) {
        /* use the work formula to compute
         * the acceptance probability */
        de = ep - ep0 + (x - x0) * (f + f0) * 0.5;
      } else if ( m->hmcacc == HMCACC_TOTAL_ENERGY ) {
        /* use the total energy change to compute
         * the acceptance probability */
        de = dee;
      } else {
        /* use the 4-potential formula to compute
         * the acceptance probability */
        ep2 = 0.5 * (x - xi0) * (x - xi0);
        ep3 = 0.5 * (x0 - xi) * (x0 - xi);
        de = 0.5 * (ep - ep0 - ep2 + ep3);
      }

      if ( de > 0 ) {
        r = rand01();
        if ( r > exp(-de/m->tp) ) {
          acc = 0;
        }
      }

      if ( !acc ) {
        if ( m->vflip ) {
          v = -v0;
        } else {
          v = v0;
        }
        x = x0;
        xi = xi0;
        f = f0;
      }
    }
#if 0
    /* velocity rescaling to adjust the velocity */
    if ( m->vflip ) {
      vrescale(m->thdt, &v, m->tp);
    } else {
      v = randgaus();
    }
#endif
    if ( x >= -m->xmax && x < m->xmax ) {
      hist[ (int) ((x + m->xmax) / m->dx) ] += 1;
    }
  }

  savehis(m, hist, nhis);
  free(hist);
  return 0;
}



int main(int argc, char **argv)
{
  model_t m[1];

  model_default(m);
  model_doargs(m, argc, argv);
  domd(m);
  return 0;
}

