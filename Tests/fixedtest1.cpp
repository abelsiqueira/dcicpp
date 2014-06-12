#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/*
 * min x_1^2 + x_2^2 + x_3^2
 * s.t. x_1 + x_2 + x_3 = 1
 *      1 <= x_2 <= 1
 *
 * sol = [0, 1, 0]
 * fsol = 1
 *
 *
 */

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = x1 * x1 + x2 * x2 + x3 * x3;
  if (*grad == dciTrue) {
    g[0] = 2 * x1;
    g[1] = 2 * x2;
    g[2] = 2 * x3;
  }
}

void CPROD (pInt, Int *, Int *, Bool * , Real * , Real *, Real * p, Real * q) {
  q[0] =  2 * p[0];
  q[1] =  2 * p[1];
  q[2] =  2 * p[2];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = x1 * x1 + x2 * x2 + x3 * x3;
  c[0] = x1 + x2 + x3 - 1;
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *,
      Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  c[0] = x1 + x2 + x3 - 1;
  if (*Grad == dciFalse)
    return;
  indvar[0] = 0;
  indfun[0] = 0;
  J[0] = 1;
  indvar[1] = 1;
  indfun[1] = 0;
  J[1] = 1;
  indvar[2] = 2;
  indfun[2] = 0;
  J[2] = 1;
  *nnzJ = 3;
}

int main () {
  Int n = 3, m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  x[0] = 2;
  x[1] = 2;
  y[0] = 2;
  bl[0] = -dciInf;
  bl[1] = 1;
  bl[2] = -dciInf;
  bu[0] = dciInf;
  bu[1] = 1;
  bu[2] = dciInf;

  cl[0] = 0;
  cu[0] = 0;
  equatn[0] = dciTrue;

  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
