#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. h(x) = 0
 *
 * where
 *
 * f(x) = sum_{i = 1}^n x_i^2
 *
 * and
 *
 * h_i(x) = x_i - x_1^2 - 1
 *
 * sol = [0; 1; 1; ...; 1; 1];
 *       lambda = [-2; -2; -2; ...; -2];
 */

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
    if (*grad == 1)
      g[i] = 2*xi;
  }
}

void CPROD (pInt, Int * n, Int * m, Bool *, Real *, Real * y, Real * p, Real * q) {
  q[0] = 0;
  for (Int j = 0; j < *m; j++)
    q[0] += y[j];
  q[0] = (2 - 2*q[0])*p[0];

  for (Int i = 1; i < *n; i++)
    q[i] = 2*p[i];
}

void CFN (pInt, Int * n, Int * m, Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
  }
  for (Int i = 0; i < *m; i++)
    c[i] = x[i+1] - x[0]*x[0] - 1;
}

void CCFSG (pInt, Int *, Int * m, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  for (Int i = 0; i < *m; i++)
    c[i] = x[i+1] - x[0]*x[0] - 1;

  if (*Grad == dciFalse)
    return;

  Int k = 0;
  Real x1 = x[0];
  if (x1 != 0) {
    for (Int i = 0; i < *m; i++) {
      J[k] = -2*x1;
      indvar[k] = 0;
      indfun[k] = i;
      k++;
    }
  }
  for (Int i = 0; i < *m; i++) {
    J[k] = 1;
    indvar[k] = i+1;
    indfun[k] = i;
    k++;
  }

  *nnzJ = k;
}

int main () {
  Int n = 30;
  Int m = n - 1;
  Int amax = 2*n;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = 2;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = dciTrue;
  }

  dci.set_amax (amax);

  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
