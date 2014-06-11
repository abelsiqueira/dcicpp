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

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 3)
    return;
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = x1 * x1 + x2 * x2 + x3 * x3;
  if (*grad == dciTrue) {
    g[0] = 2 * x1;
    g[1] = 2 * x2;
    g[2] = 2 * x3;
  }
}

void CPROD (Int * n, Int * m, Bool * , Real * , Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*n != 3) || (*m != 1) || (*mmax < *m) )
    return;
  Real y1 = y[0];
  q[0] =  2 * p[0];
  q[1] =  2 * p[1];
  q[2] =  2 * p[2];
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  if ( (*n != 3) || (*m != 1) || (*mmax < *m) )
    throw("Error");
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = x1 * x1 + x2 * x2 + x3 * x3;
  c[0] = x1 + x2 + x3 - 1;
}

void CCFSG (Int * , Int * , Real * x, Int * , Real * c, Int * nnzJ,
    Int * , Real * J, Int * indvar, Int * indfun, Bool * Grad) {
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

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_equatn (m, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
