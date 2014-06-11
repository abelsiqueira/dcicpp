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

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
    if (*grad == 1)
      g[i] = 2*xi;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real unused = x[0];
  unused = x[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != *n - 1) || (*mmax < *m) )
      return;
    q[0] = 0;
    for (Int j = 0; j < *m; j++)
      q[0] += y[j];
    q[0] = (2 - 2*q[0])*p[0];

    for (Int i = 1; i < *n; i++)
      q[i] = 2*p[i];

  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
  }
  if (*m != *n - 1)
    return;
  if (*mmax < *m)
    return;
  for (Int i = 0; i < *m; i++)
    c[i] = x[i+1] - x[0]*x[0] - 1;
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  if (*m != *n - 1)
    return;

  for (Int i = 0; i < *m; i++)
    c[i] = x[i+1] - x[0]*x[0] - 1;

  if (*Grad == dciFalse)
    return;
  if (*mmax < *m)
    return;
  if (*jmax < 2*(*m))
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
    cl[i] = -dciInf;
    cu[i] = dciInf;
  }

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_amax (amax);

  dci.start ();
  dci.solve ();
  dci.show();

}
