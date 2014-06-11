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

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 2)
    return;
  Real x1 = x[0], x2 = x[1];
  *f = pow (1 - x1, 2) + rosencst * pow (x2 - x1*x1, 2);
  g[0] = 2*(x1 - 1) - 4 * rosencst * x1 * (x2 - x1*x1);
  g[1] = 2 * rosencst * (x2 - x1*x1);
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if (*n != 2)
    return;
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != 1) || (*mmax < *m) )
      return;
    q[0] = (2 - 4 * rosencst*(x[1] - x[0]*x[0]) + 8 * rosencst *x[0]*x[0] + y[0]*(-2 + M_PI*M_PI*sin(M_PI*x[0])) )*p[0] - 4 * rosencst * x[0] * p[1];
    q[1] = -4 * rosencst*x[0]*p[0] + 2 * rosencst *p[1];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = pow (1 - x1, 2) + rosencst * pow (x2 - x1*x1, 2);
  if (*m != 1)
    return;
  if (*mmax < 1)
    return;
  c[0] = x2 - x1*x1 - sin(M_PI*x1);
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  if (*n != 2)
    return;
  if (*m != 1)
    return;

  c[0] = x[1] - x[0]*x[0] - sin(M_PI*x[0]);

  if (*Grad == dciFalse)
    return;
  if (*mmax < 1)
    return;
  if (*jmax < 0)
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
  }

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);

  dci.start ();
  dci.solve ();
  dci.show();

}
