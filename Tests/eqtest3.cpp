#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. h(x) = 0
 *
 * where
 *
 * f(x) = - sum_{i = 1}^n x_i^2
 *
 * and
 *
 * h(x) = sum_{i = 1}^n (x_i - 1)^2 - 1
 *
 * sol = (1 + sqrt(n))/sqrt(n)
 *       lambda = 1 - sqrt(n)
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f -= xi*xi;
    if (*grad == 1)
      g[i] = -2*xi;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real unused = x[0];
  unused = x[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != 1) || (*mmax < *m) )
      return;
    for (Int i = 0; i < *n; i++) {
      q[i] = (-2 + 2 * y[0]) * p[i];
    }
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f -= xi*xi;
  }
  if (*m != 1)
    return;
  if (*mmax < 1)
    return;
  c[0] = -1;
  for (Int i = 0; i < *n; i++)
    c[0] += (x[i] - 1) * (x[i] - 1);
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  if (*m != 1)
    return;

  c[0] = -1;
  for (Int i = 0; i < *n; i++)
    c[0] += (x[i] - 1) * (x[i] - 1);

  if (*Grad == dciFalse)
    return;
  if (*mmax < 1)
    return;
  if (*jmax < 0)
    return;

  Int k = 0;
  Real xi = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    if (xi != 1) {
      J[k] = 2*(xi - 1);
      indvar[k] = i;
      indfun[k] = 0;
      k++;
    }
  }

  *nnzJ = k;
}

int main () {
  Int n = 20000;
  Int m = 1;
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

  dci.start ();
  dci.solve ();
  dci.show (1);

}
