#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Let's make the problem HS34
 * min -x1
 * s.t. cI(x) => 0
 *      l <= x <= u
 *
 * where
 *
 * f(x) = - x(1)
 *
 * and
 *
 * cI(x) = [x(2) - exp(x(1));
 *          x(3) - exp(x(2))]
 *
 * l = [0; -inf; -inf]
 * u = [inf; inf; 10]
 *
 * Expected solution: 
 *
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 3)
    return;
  Real x1 = x[0];
  *f = -x1;
  if (*grad == dciTrue) {
    g[0] = -x1;
    g[1] = 0;
    g[2] = 0;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*n != 3) || (*m != 2) || (*mmax < *m) )
    return;
  Real x1 = x[0], x2 = x[1];
  Real y1 = y[0], y2 = y[2];
  q[0] = - y1 * exp(x1) * p[0];
  q[1] = - y2 * exp(x2) * p[1];
  q[2] = 0;
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = -x1;
  c[0] = x2 - exp(x1);
  c[1] = x3 - exp(x2);
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  c[0] = x2 - exp(x1);
  c[1] = x3 - exp(x2);
  if (*Grad == dciFalse)
    return;
  indvar[0] = 0;
  indfun[0] = 0;
  J[0] = -exp(x1);
  indvar[1] = 1;
  indfun[1] = 1;
  J[1] = -exp(x2);
  indvar[2] = 1;
  indfun[2] = 0;
  J[2] = 1;
  indvar[3] = 2;
  indfun[3] = 1;
  J[3] = 1;
  *nnzJ = 4;
}

int main () {
  Int n = 3, m = 2;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  x[0] = 0;
  x[1] = 1.05;
  x[2] = 2.9;
  y[0] = 0;
  y[1] = 0;

  bl[0] = 0;
  bl[1] = -dciInf;
//  bl[1] = 0;
  bl[2] = -dciInf;
//  bl[2] = 0;
//  bu[0] = 100;
//  bu[1] = 100;
  bu[0] = dciInf;
  bu[1] = dciInf;
  bu[2] = 10;
  
  cl[0] = 0;
  cl[1] = 0;
  cl[2] = 0;
  cu[0] = dciInf;
  cu[1] = dciInf;
  cu[2] = dciInf;

  equatn[0] = dciFalse;
  equatn[1] = dciFalse;
  equatn[2] = dciFalse;

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_equatn (m, equatn);

  dci.start ();
  dci.solve ();
  dci.show (1);

}
