#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Let's make the problem HS55
 * min  x1 + 2x2 + 4x5 + exp(x1x4)
 * s.t. x1 + 2x2 +           5x5      = 6
 *      x1 +  x2 +  x3                = 3
 *                       x4 + x5 + x6 = 2
 *      x1             + x4           = 1
 *            x2            + x5      = 2
 *                  x3           + x6 = 2
 *      x >= 0
 *      x1 <= 1, x4 <= 1
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 6)
    return;
  Real x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
  *f = x1 + 2*x2 + 4*x5 + exp(x1*x4);
  if (*grad == dciTrue) {
    g[0] = 1 + x4*exp(x1*x4);
    g[1] = 2;
    g[2] = 0;
    g[3] = x1*exp(x1*x4);
    g[4] = 4;
    g[5] = 0;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*n != 6) || (*m != 6) || (*mmax < *m) )
    return;
  Real x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
  Real y1 = y[0], y2 = y[1], y3 = y[2], y4 = y[3], y5 = y[4], y6 = y[5];
  q[0] = x4*x4*exp(x1*x4) * p[0];
  q[1] = 0;
  q[2] = 0;
  q[3] = x1*x1*exp(x1*x4) * p[3];
  q[4] = 0;
  q[5] = 0;
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
  *f = x1 + 2*x2 + 4*x5 + exp(x1*x4);
  c[0] = x1 + 2*x2 + 5*x5 - 6;
  c[1] = x1 + x2 + x3 - 3;
  c[2] = x4 + x5 + x6 - 2;
  c[3] = x1 + x4 - 1;
  c[4] = x2 + x5 - 2;
  c[5] = x3 + x6 - 2;
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3], x5 = x[4], x6 = x[5];
  c[0] = x1 + 2*x2 + 5*x5 - 6;
  c[1] = x1 + x2 + x3 - 3;
  c[2] = x4 + x5 + x6 - 2;
  c[3] = x1 + x4 - 1;
  c[4] = x2 + x5 - 2;
  c[5] = x3 + x6 - 2;
  if (*Grad == dciFalse)
    return;
  /*
   * A = |1 2 0 0 5 0|
   *     |1 1 1 0 0 0|
   *     |0 0 0 1 1 1|
   *     |1 0 0 1 0 0|
   *     |0 1 0 0 1 0|
   *     |0 0 1 0 0 1|
   */
  indvar[0] = 0;
  indfun[0] = 0;
  J[0] = 1;
  indvar[1] = 1;
  indfun[1] = 0;
  J[1] = 2;
  indvar[2] = 4;
  indfun[2] = 0;
  J[2] = 5;

  indvar[3] = 0;
  indfun[3] = 1;
  J[3] = 1;
  indvar[4] = 1;
  indfun[4] = 1;
  J[4] = 1;
  indvar[5] = 2;
  indfun[5] = 1;
  J[5] = 1;

  indvar[6] = 3;
  indfun[6] = 2;
  J[6] = 1;
  indvar[7] = 4;
  indfun[7] = 2;
  J[7] = 1;
  indvar[8] = 5;
  indfun[8] = 2;
  J[8] = 1;

  indvar[9] = 0;
  indfun[9] = 3;
  J[9] = 1;
  indvar[10] = 3;
  indfun[10] = 3;
  J[10] = 1;

  indvar[11] = 1;
  indfun[11] = 4;
  J[11] = 1;
  indvar[12] = 4;
  indfun[12] = 4;
  J[12] = 1;

  indvar[13] = 2;
  indfun[13] = 5;
  J[13] = 1;
  indvar[14] = 5;
  indfun[14] = 5;
  J[14] = 1;

  *nnzJ = 15;
}

int main () {
  Int n = 6, m = 6;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m], linear[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  x[0] = 1;
  x[1] = 2;
  x[2] = 0;
  x[3] = 0;
  x[4] = 0;
  x[5] = 2;

  for (Int i = 0; i < 6; i++) {
    bl[i] = 0;
    bu[i] = dciInf;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = dciTrue;
    linear[i] = dciTrue;
  }

  bu[0] = 1;
  bu[3] = 1;

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
//  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_equatn (m, equatn);
  dci.set_linear (m, linear);

  dci.start ();
  dci.solve ();
  dci.show();

}
