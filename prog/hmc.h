#ifndef HMC_H__
#define HMC_H__



/* Hybrid Monte Carlo */


typedef struct {
  int n;
  double (*x)[D];
  double (*v)[D];
  double (*f)[D];
} hmc_t;



__inline static hmc_t *hmc_open(int n)
{
  hmc_t *h;

  xnew(h, 1);
  h->n = n;
  xnew(h->x, n);
  xnew(h->v, n);
  xnew(h->f, n);
  return h;
}



__inline static void hmc_close(hmc_t *h)
{
  free(h->x);
  free(h->v);
  free(h->f);
  free(h);
}



static void hmc_push(hmc_t *h, double (*x)[D], double (*v)[D], double (*f)[D])
{
  memcpy(h->x, x, sizeof(x[0]) * h->n);
  memcpy(h->v, v, sizeof(v[0]) * h->n);
  memcpy(h->f, f, sizeof(f[0]) * h->n);
}



static void hmc_pop(hmc_t *h, double (*x)[D], double (*v)[D], double (*f)[D],
    int reversev)
{
  int i, n = h->n;

  if ( reversev ) {
    for ( i = 0; i < n; i++ )
      vinv( h->v[i] );
  }
  memcpy(x, h->x, sizeof(x[0]) * n);
  memcpy(v, h->v, sizeof(v[0]) * n);
  memcpy(f, h->f, sizeof(f[0]) * n);
}




#endif /* HMC_H__ */
