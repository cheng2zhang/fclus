/* implicit Gibbs-ensemble simulation */
#define D 3
#include "ljdiv.h"



int n = 108;
double rho = 0.3;
double tp = 1.2;
double rcdef = 1e9;
double amp = 0.2;

double nsteps = 1e10;
double nstrep = 100000;

const char *fnpos = "lj.pos";
const char *fnhist = "hist.dat";

double mctot = 0, mcacc = 0;
double divtot = 0, divacc = 0;



double *hist;

static int writehist(double *h, int nh, const char *fn)
{
  int i;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  for ( i = 1; i <= nh; i++ ) {
    fprintf(fp, "%d %g\n", i, h[i]);
  }
  fclose(fp);
  return 0;
}



int main(void)
{
  double t;
  lj_t *lj;
  ljdiv_t *d;
  int divsize;
  double epsm = 0;

  /* make a Lennard-Jones object */
  lj = lj_open(n, rho, rcdef);
  d = ljdiv_open(lj);
  /* turn on the division energy */
  d->lamdiv = 0.0;

  xnew(hist, n + 1);

  /* compute the potential energy */
  lj_energy(lj);
  d->ediv = ljdiv_ediv(d, d->idiv, d->divseed, d->divr);

  for ( t = 1; t <= nsteps; t++ ) {
    /* a step of Metropolis algorithm
     * graph is implicitly computed */
    mctot += 1;
    mcacc += ljdiv_metro(d, amp, 1/tp);
    /* try to change the division */
    divtot += 1;
    divacc += ljdiv_changediv(d, lj->l * 0.1, 1/tp);
    //printf("%d, %g, %g\n", lj_getdivsize(lj, lj->idiv), lj->ediv, lj_ediv(lj, lj->idiv)); getchar();

    /* add the division size to histogram */
    divsize = ljdiv_getdivsize(d, d->idiv);

    hist[divsize] += 1;

    if ( fmod(t, nstrep) < 0.1 ) {
      lj_writepos(lj, lj->x, lj->v, fnpos);
      writehist(hist, n, fnhist);
      fprintf(stderr, "t %g, mcacc %.2f%% ep %g, divacc %.2f%%, ediv %g, "
          "divsize %d, divr %g, divrho %g\n",
          t, 100*mcacc/mctot, lj->epot,
          100*divacc/divtot, d->ediv,
          divsize, d->divr, divsize/lj_sphrvol(d->divr));
    }
    epsm += lj->epot;
  }
  //printf("epot: %g, %g\n", lj->epot, lj_energy_low(lj, lj->x, lj->r2ij, NULL, NULL, NULL));
  lj_close(lj);
  ljdiv_close(d);
  free(hist);
  fprintf(stderr, "rho %g, tp %g, ep %g\n", rho, tp, epsm/nsteps/n);
  return 0;
}

