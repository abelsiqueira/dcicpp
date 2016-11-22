#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem HS18
 * min f(x)
 * s.t. cI(x) >= 0
 *      2 <= x(1) <= 50
 *      0 <= x(2) <= 50
 *
 * where
 *
 * f(x) = 0.01*x(1)^2 + x(2)^2
 *
 * and
 *
 * cI(x) = [x(1) * x(2) - 25; x(1)^2 + x(2)^2 - 25]
 *
 * Expected solution:
 *
 * x = [sqrt(10)*5, sqrt(10)/2] ~ [15.811, 1.5811]
 * f(x) = 5
 *
 */

void COFG (pInt, const Int *, const Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  if (*grad == dciTrue) {
    g[0] = 0.02 * x1;
    g[1] = 2 * x2;
  }
}

void CPROD (pInt, const Int *, const Int *, const Bool *, const Real *, const Real * y, Real * p, Real * q) {
  Real y1 = y[0], y2 = y[1];
  q[0] = (0.02 + 2*y2) * p[0] + y1 * p[1];
  q[1] = y1 * p[0] + (2 + 2*y2) * p[1];
}

void CFN (pInt, const Int *, const Int *, const Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  c[0] = x1 * x2 - 25;
  c[1] = x1 * x1 + x2 * x2 - 25;
}

void CCFSG (pInt, const Int *, const Int *, const Real * x, Real * c, Int * nnzJ, const Int *, Real * J,
    Int * indvar, Int * indfun, const Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x1 * x2 - 25;
  c[1] = x1 * x1 + x2 * x2 - 25;
  if (*Grad == dciFalse)
    return;
  if (x1 == 0) {
    if (x2 == 0) {
      *nnzJ = 0;
    } else {
      *nnzJ = 2;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
      indvar[1] = 1;
      indfun[1] = 1;
      J[1] = 2*x2;
    }
  } else {
    if (x2 == 0) {
      *nnzJ = 2;
      indvar[0] = 0;
      indfun[0] = 1;
      J[0] = 2*x1;
      indvar[1] = 1;
      indfun[1] = 0;
      J[1] = x1;
    } else {
      *nnzJ = 4;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
      indvar[1] = 0;
      indfun[1] = 1;
      J[1] = 2*x1;
      indvar[2] = 1;
      indfun[2] = 0;
      J[2] = x1;
      indvar[3] = 1;
      indfun[3] = 1;
      J[3] = 2*x2;
    }
  }
}

int main () {
  Int n = 2, m = 2;
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
  bl[0] = 2;
  bl[1] = 0;
  bu[0] = 50;
  bu[1] = 50;

  y[0] = 0;
  y[1] = 0;
  cl[0] = 0;
  cl[1] = 0;
  cu[0] = dciInf;
  cu[1] = dciInf;
  equatn[0] = dciFalse;
  equatn[1] = dciFalse;

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
