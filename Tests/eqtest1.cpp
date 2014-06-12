#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. h(x) = 0
 *
 * where
 *
 * f(x) = x1^2 + 2 * x2^2
 *
 * and
 *
 * h(x) = x1^2 + (x2 - 1)^2/4 - 1
 *
 * Expected solution:
 *
 * x = [sqrt(33)/7, -1/7];
 * lambda = -1;
 *
 */

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += (i + 1)*xi*xi;
    if (*grad == 1)
      g[i] = 2*xi*(i + 1);
  }
}

void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  q[0] = 2*(1 + y[0])*p[0];
  q[1] = 0.5*(8 + y[0])*p[1];
}

void CFN (pInt, Int * n, Int *, Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += (i + 1)*xi*xi;
  }
  c[0] = x[0]*x[0] + (x[1] - 1)*(x[1] - 1)/4 - 1;
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = x[0]*x[0] + (x[1] - 1)*(x[1] - 1)/4 - 1;
  if (*Grad == dciFalse)
    return;
  Int i = 0;
  if (x[0] != 0) {
    J[i] = 2*x[0];
    indvar[i] = 0;
    indfun[i] = 0;
    i++;
  }
  if (x[1] != 1) {
    J[i] = (x[1] - 1)/2;
    indvar[i] = 1;
    indfun[i] = 0;
    i++;
  }
  *nnzJ = i;
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

  for (Int i = 0; i < n; i++) {
    x[i] = (i + 1)*(i + 2);
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (Int i = 0; i < m; i++) {
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
