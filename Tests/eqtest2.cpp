#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. h(x) = 0
 *
 * where
 *
 * f(x) = -x1^4 -x2^4
 *
 * and
 *
 * h(x) = x1^2 + (x2 - 1)^2/4 - 1
 *
 * sol = [0; 3], f = -81
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    xi2 = xi*xi;
    *f -= xi2*xi2;
    if (*grad == 1)
      g[i] = -4*xi*xi2;
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real unused = x[0];
  unused = x[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*n != 2) || (*m != 1) || (*mmax < *m) )
      return;
    q[0] = ( -12*x[0]*x[0] + 2 * y[0]) * p[0];
    q[1] = ( -12*x[1]*x[1] + 0.5 * y[0]) * p[1];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    xi2 = xi*xi;
    *f -= xi2*xi2;
  }
  if (*m != 1)
    return;
  if (*mmax < 1)
    return;
  c[0] = x[0]*x[0] + (x[1] - 1)*(x[1] - 1)/4 - 1;
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = x[0]*x[0] + (x[1] - 1)*(x[1] - 1)/4 - 1;
  if (*Grad == dciFalse)
    return;
  Int i = 0;
  if (*n != 2)
    return;
  if (*m != 1)
    return;
  if (*mmax != 1)
    return;
  if (*jmax < 0)
    return;
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

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = (i + 1)*(i + 2);
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
