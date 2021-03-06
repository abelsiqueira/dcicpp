#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Loosened HS34
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
 * sol = [log(log(10)), log(10), 10]
 *     ~ [0.83403, 2.3026, 10]
 * fsol = -log(log(10))
 *
 */

void COFG (pInt, const Int *, const Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0];
  *f = -x1;
  if (*grad == dciTrue) {
    g[0] = -x1;
    g[1] = 0;
    g[2] = 0;
  }
}

void CPROD (pInt, const Int *, const Int *, const Bool *, const Real * x, const Real * y, Real * p, Real * q) {
  Real x1 = x[0], x2 = x[1];
  Real y1 = y[0], y2 = y[1];
  q[0] = - y1 * exp(x1) * p[0];
  q[1] = - y2 * exp(x2) * p[1];
  q[2] = 0;
}

void CFN (pInt, const Int *, const Int *, const Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1], x3 = x[2];
  *f = -x1;
  c[0] = x2 - exp(x1);
  c[1] = x3 - exp(x2);
}

void CCFSG (pInt, const Int *, const Int *, const Real * x, Real * c, Int * nnzJ, const Int *, Real * J,
    Int * indvar, Int * indfun, const Bool * Grad) {
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

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
