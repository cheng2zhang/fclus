#ifndef MDTIL_H__
#define MDUTIL_H__



#include "mat.h"



/* remove the center of mass motion */
static void rmcom(double (*x)[D], const double *m, int n)
{
  int i;
  double xc[D] = {0}, mtot = 0;

  for ( i = 0; i < n; i++ ) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double xc[D] = {0, 0, 0}, xi[D], ang[D], am[D] = {0, 0, 0};
  double dv[D], mat[D][D], inv[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  double mtot = 0;

  for (i = 0; i < n; i++) {
    vsinc(xc, x[i], m[i]);
    mtot += m[i];
  }
  vsmul(xc, 1.0 / mtot);

  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(ang, xi, v[i]);
    vsinc(am, ang, m[i]);

    xx += m[i] * xi[0] * xi[0];
    yy += m[i] * xi[1] * xi[1];
    zz += m[i] * xi[2] * xi[2];
    xy += m[i] * xi[0] * xi[1];
    yz += m[i] * xi[1] * xi[2];
    zx += m[i] * xi[2] * xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  minv(inv, mat);
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  /* ang is the solution of M^(-1) * L */
  for (i = 0; i < n; i++) {
    vdiff(xi, x[i], xc);
    vcross(dv, ang, xi);
    vinc(v[i], dv);
  }
}



#endif /* MDUTIL_H__ */

