#include <iostream>
#include "dci.h"

using namespace DCI;

/* Infeasible problem
 *
 * min f(x)
 * s.t. cE(x) = 0
 *
 * where
 *
 * f(x) = x(1)^2 + x(2)^2
 *
 * and
 *
 * cE(x) = [x(2) - x(1)^2 - 1;
 *          x(1) - x(2)^2 - 1]
 *
 *
 */

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 2)
    return;
  Real x1 = x[0], x2 = x[1];
  *f = x1 * x1 + x2 * x2;
  if (*grad == dciTrue) {
    g[0] = 2 * x1;
    g[1] = 2 * x2;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  Real y1 = y[0], y2 = y[1];
  q[0] = 2*(1 - y1)*p[0];
  q[1] = 2*(1 - y2)*p[1];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = x1 * x1 + x2 * x2;
  c[0] = x2 - x1*x1 - 1;
  c[1] = x1 - x2*x2 - 1;
}

void CCIFG (pInt, Int *, Int * i, Real * x, Real * ci, Real * gci, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  if (*i == 1) {
    *ci = x2 - x1*x1 - 1;
    if (*Grad == dciTrue) {
      gci[0] = -2*x1;
      gci[1] = 1;
    }
  } else if (*i == 2) {
    *ci = x1 - x2*x2 - 1;
    if (*Grad == dciTrue) {
      gci[0] = 1;
      gci[1] = -2*x2;
    }
  } else {
    throw ("CCIFG out of bound");
  }
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x2 - x1*x1 - 1;
  c[1] = x1 - x2*x2 - 1;
  if (*Grad == dciFalse)
    return;

  int k = 0;
  if (x1 != 0) {
    indvar[k] = 0;
    indfun[k] = 0;
    J[0] = -2*x1;
    k++;
  }
  if (x2 != 0) {
    indvar[k] = 1;
    indfun[k] = 1;
    J[1] = -2*x2;
    k++;
  }
  indvar[k] = 1;
  indfun[k] = 0;
  J[k] = 1;
  k++;
  indvar[k] = 0;
  indfun[k] = 1;
  J[k] = 1;
  *nnzJ = k + 1;
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
  dci.set_ccifg (CCIFG);

  for (Int i = 0; i < 2; i++) {
    x[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
    y[i] = 0;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = dciTrue;
  }


  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
