#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) => 0
 *
 * where
 *
 * f(x) = 0.01*x(1)^2 + x(2)^2
 *
 * and
 *
 * cI(x) = x(1) * x(2) - 25
 *
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 2)
    return;
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  if (*grad == dciTrue) {
    g[0] = 0.02 * x1;
    g[1] = 2 * x2;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*n != 2) || (*m != 1) || (*mmax < *m) )
    return;
  Real y1 = y[0];
  q[0] = 0.02 * p[0] + y1 * p[1];
  q[1] = y1 * p[0] + 2 * p[1];
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  if ( (*n < *m) || (*m != 1) || (*mmax < *m) )
    throw("Error");
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  c[0] = x1 * x2 - 25;
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
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
  y[0] = 1;
  bu[1] = dciInf;
  bl[0] = 2;
  bl[1] = 0;
  bu[0] = 50;
  
  cl[0] = 0;
  cu[0] = dciInf;
  equatn[0] = dciFalse;

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
