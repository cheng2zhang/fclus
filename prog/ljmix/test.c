#include "ljmixmodel.h"
#include "ljmixcore.h"



/* test if force matches the energy */
static void ljmix_testforce(ljmix_t *lj)
{
  int i;
  double ep, f2 = 0, delta = 0.0001;

  ep = ljmix_force(lj);
  for ( i = 0; i < lj->n; i++ ) {
    f2 += vsqr(lj->f[i]);
  }
  /* move along the force */
  for ( i = 0; i < lj->n; i++ ) {
    vsinc(lj->x[i], lj->f[i], delta/f2);
  }
  ljmix_energy(lj);
  printf("energy %g, delta %g\n", lj->epot, (ep - lj->epot)/delta);
}



int main(int argc, char **argv)
{
  ljmixmodel_t m[1];
  ljmix_t *lj;

  ljmixmodel_default(m);
  ljmixmodel_doargs(m, argc, argv);
  printf("ns %d, np %d %d, sig %g %g\n", m->ns, m->np[0], m->np[1], m->sig[0], m->sig[1]);

  lj = ljmix_open(m->ns, m->np, m->sig, m->eps, m->rho, m->rcdef);
  ljmix_testforce(lj);
  ljmix_close(lj);
  return 0;
}
