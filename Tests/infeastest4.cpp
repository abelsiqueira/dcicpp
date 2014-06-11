#include <iostream>
#include "dci.h"

using namespace DCI;

/* Infeasible problem
 *
 * min f(x)
 * s.t. cI(x) >= 0
 *
 * where
 *
 * f(x) = x(1) + x(2)
 *
 * and
 *
 * cE(x) = [x(1)^2 + x(2)^2 - 1;    <= 0
 *          x(1)^2 + x(2)^2 - 4]    >= 0
 *
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 2)
    return;
  Real x1 = x[0], x2 = x[1];
  *f = 7*x1 + 11*x2;
  if (*grad == dciTrue) {
    g[0] = 7;
    g[1] = 11;
  }
}

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real y1 = y[0], y2 = y[1];
  q[0] = 2*(y1 + y2)*p[0];
  q[1] = 2*(y1 + y2)*p[1];
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = 7*x1 + 11*x2;
  c[0] = x1*x1 + x2*x2 - 1;
  c[1] = c[0] - 3;
  c[2] = 3*x1 + 4*x2 - 1;
}

void CCIFG (Int *, Int * i, Real * x, Real * ci, Real * gci, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  if (*i <= 2) {
    *ci = x1*x1 + x2*x2 - 1;
    if (*Grad == dciTrue) {
      gci[0] = 2*x1;
      gci[1] = 2*x2;
    }
    if (*i == 2) {
      *ci -= 3;
    }
  } else {
    *ci = 3*x1 + 4*x2 - 1;
    if (*Grad == dciTrue) {
      gci[0] = 3;
      gci[1] = 4;
    }
  }
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x1*x1 + x2*x2 - 1;
  c[1] = c[0] - 3;
  c[2] = 3*x1 + 4*x2 - 1;
  if (*Grad == dciFalse)
    return;

  int k = 2;
  indvar[0] = 0;
  indfun[0] = 2;
  J[0] = 3;
  indvar[1] = 1;
  indfun[1] = 2;
  J[1] = 4;
  if (x1 != 0) {
    indvar[k] = 0;
    indfun[k] = 0;
    J[k] = 2*x1;
    k++;
    indvar[k] = 0;
    indfun[k] = 1;
    J[k] = 2*x1;
    k++;
  }
  if (x2 != 0) {
    indvar[k] = 1;
    indfun[k] = 0;
    J[k] = 2*x2;
    k++;
    indvar[k] = 1;
    indfun[k] = 1;
    J[k] = 2*x2;
    k++;
  }
  *nnzJ = k;
}

int main () {
  Int n = 2, m = 3;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);
  dci.set_ccifg (CCIFG);

  for (Int i = 0; i < 2; i++) {
    x[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
    y[i] = 0;
    equatn[i] = dciFalse;
  }
  cl[0] = -dciInf;
  cu[0] = 0;
  cl[1] = 0;
  cu[1] = dciInf;
  cl[2] = 0;
  cu[2] = 0;
  equatn[2] = dciTrue;

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
