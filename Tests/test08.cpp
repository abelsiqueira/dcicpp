#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cE(x) = 0
 *      x1 >= 2
 *      x2 >= 0
 *
 * where
 *
 * f(x) = 0.01*x(1)^2 + x(2)^2
 *
 * and
 *
 * cI(x) = x(1) * x(2) - 25
 *
 * Expected solution:
 *
 * x = [sqrt(10)*5, sqrt(10)/2] ~ [15.811, 1.5811]
 * f(x) = 5
 *
 */

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  if (*grad == dciTrue) {
    g[0] = 0.02 * x1;
    g[1] = 2 * x2;
  }
}

void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  Real y1 = y[0];
  q[0] = 0.02 * p[0] + y1 * p[1];
  q[1] = y1 * p[0] + 2 * p[1];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  c[0] = x1 * x2 - 25;
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J,
    Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x1 * x2 - 25;
  if (*Grad == dciFalse)
    return;
  if (x1 == 0) {
    if (x2 == 0) {
      *nnzJ = 0;
    } else {
      *nnzJ = 1;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
    }
  } else {
    if (x2 == 0) {
      *nnzJ = 1;
      indvar[0] = 1;
      indfun[0] = 0;
      J[0] = x1;
    } else {
      *nnzJ = 2;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
      indvar[1] = 1;
      indfun[1] = 0;
      J[1] = x1;
    }
  }
}

int main () {
  Int n = 2, m = 1;
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
  bu[0] = dciInf;
  bu[1] = dciInf;
  bl[0] = 2;
  bl[1] = 0;
  
  y[0] = 1;
  cl[0] = 0;
  cu[0] = 0;
  equatn[0] = dciTrue;

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
