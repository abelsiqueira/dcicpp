#include <iostream>
#include "dci.h"
#include <cmath>

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. h(x) = 0
 *
 * where
 *
 * f(x) = (1 - x_1)^2 + 100 * (x_2 - x_1^2)^2
 *
 * and
 *
 * h(x) = x_2 - x_1^2 - sin(pi * x_1)
 *
 * sol = [1; 1]
 *       lambda = 0
 */

Real rosencst = 100;

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1];
  *f = pow (1 - x1, 2) + rosencst * pow (x2 - x1*x1, 2);
  if (*grad) {
    g[0] = 2*(x1 - 1) - 4 * rosencst * x1 * (x2 - x1*x1);
    g[1] = 2 * rosencst * (x2 - x1*x1);
  }
}

void CPROD (pInt, Int *, Int *, Bool *, Real * x, Real * y, Real * p, Real * q) {
  q[0] = (2 - 4 * rosencst*(x[1] - x[0]*x[0]) + 8 * rosencst *x[0]*x[0] + y[0]*(-2 + M_PI*M_PI*sin(M_PI*x[0])) )*p[0] - 4 * rosencst * x[0] * p[1];
  q[1] = -4 * rosencst*x[0]*p[0] + 2 * rosencst *p[1];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = pow (1 - x1, 2) + rosencst * pow (x2 - x1*x1, 2);
  c[0] = x2 - x1*x1 - sin(M_PI*x1);
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {

  c[0] = x[1] - x[0]*x[0] - sin(M_PI*x[0]);

  if (*Grad == dciFalse)
    return;

  Int k = 0;
  Real aux = 2*x[0] + M_PI*cos(M_PI*x[0]);
  if (aux != 0) {
    J[k] = -aux;
    indvar[k] = 0;
    indfun[k] = 0;
    k++;
  }
  J[k] = 1;
  indvar[k] = 1;
  indfun[k] = 0;
  k++;

  *nnzJ = k;
}

int main () {
  Int n = 2;
  Int m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = 1.8;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf;
    cu[i] = dciInf;
    equatn[i] = dciTrue;
  }


  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
