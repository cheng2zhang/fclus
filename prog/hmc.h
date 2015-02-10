#ifndef HMC_H__
#define HMC_H__



/* Hybrid Monte Carlo */


typedef struct {
  int n;
  double (*x)[D];
  double (*v)[D];
  double (*f)[D];
  int ni; /* number of extra integers */
  int nf; /* number of extra floating points */
  int *idat;
  double *fdat;
} hmc_t;



__inline static hmc_t *hmc_open(int n, int ni, int nf)
{
  hmc_t *h;

  xnew(h, 1);
  h->n = n;
  xnew(h->x, n);
  xnew(h->v, n);
  xnew(h->f, n);
  h->ni = ni;
  if ( ni > 0 ) {
    xnew(h->idat, ni);
  }
  h->nf = nf;
  if ( nf > 0 ) {
    xnew(h->fdat, nf);
  }
  return h;
}



__inline static void hmc_close(hmc_t *h)
{
  free(h->x);
  free(h->v);
  free(h->f);
  free(h->idat);
  free(h->fdat);
  free(h);
}



static void hmc_push(hmc_t *h, double (*x)[D], double (*v)[D], double (*f)[D],
    const int *idat, const double *fdat)
{
  memcpy(h->x, x, sizeof(x[0]) * h->n);
  memcpy(h->v, v, sizeof(v[0]) * h->n);
  memcpy(h->f, f, sizeof(f[0]) * h->n);
  if ( h->ni > 0 ) {
    memcpy(h->idat, idat, sizeof(idat[0]) * h->ni);
  }
  if ( h->nf > 0 ) {
    memcpy(h->fdat, fdat, sizeof(fdat[0]) * h->nf);
  }
}



static void hmc_pop(hmc_t *h, double (*x)[D], double (*v)[D], double (*f)[D],
    int *idat, double *fdat, int reversev)
{
  int i, n = h->n;

  if ( reversev ) {
    for ( i = 0; i < n; i++ ) {
      vneg( h->v[i] );
    }
  }
  memcpy(x, h->x, sizeof(x[0]) * n);
  memcpy(v, h->v, sizeof(v[0]) * n);
  memcpy(f, h->f, sizeof(f[0]) * n);
  if ( h->ni > 0 ) {
    memcpy(idat, h->idat, sizeof(idat[0]) * h->ni);
  }
  if ( h->nf > 0 ) {
    memcpy(fdat, h->fdat, sizeof(fdat[0]) * h->nf);
  }
}




#endif /* HMC_H__ */

