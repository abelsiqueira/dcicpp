#include <iostream>
#include "dci.h"
#include <cmath>

using namespace DCI;

/* Infeasible problem
 *
 * min f(x)
 * s.t. cE(x) = 0
 *
 * where
 *
 * f(x) = x(1)
 *
 * and
 *
 * cE(x) = [ x(1)^2 - 1 ]^2 + 1
 *
 *
 */

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0];
  *f = x1;
  if (*grad == dciTrue) {
    g[0] = 1;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int *, Int *, Bool *, Real * x, Real * y, Real * p, Real * q) {
  Real x1 = x[0];
  q[0] = (12*x1*x1 - 4)*y[0]*p[0];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0];
  *f = x1;
  c[0] = pow(x1*x1 - 1, 2) + 1;
}

void CCIFG (pInt, Int *, Int * i, Real * x, Real * ci, Real * gci, Bool * Grad) {
  Real x1 = x[0];
  if (*i == 1) {
    *ci = pow(x1*x1 - 1, 2) + 1;
    if (*Grad == dciTrue) {
      gci[0] = 4*x1*(x1*x1 - 1);
    }
  } else {
    throw ("CCIFG out of bound");
  }
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0];
  c[0] = pow(x1*x1 - 1, 2) + 1;
  if (*Grad == dciFalse)
    return;

  int k = 0;
  if (x1 != 0 && x1 != 1 && x1 != -1) {
    indvar[0] = 0;
    indfun[0] = 0;
    J[0] = 4*x1*(x1*x1 - 1);
    k++;
  }
  *nnzJ = k;
}

int main () {
  Int n = 1, m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);
  dci.set_ccifg (CCIFG);

  for (Int i = 0; i < 1; i++) {
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
