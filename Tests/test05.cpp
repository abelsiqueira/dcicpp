#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) <= 0
 *
 * where
 *
 * f(x) = dot(x,x)
 *
 * and
 *
 * cI(x) = -x.^3
 *
 * Expected solution:
 *
 * sol = [0, 0, 0, ..., 0]
 * fsol = 0
 * lambda = any
 *
 */

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
    if (*grad == dciTrue)
      g[i] = 2*xi;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int * n, Int *, Bool *, Real * x, Real * y, Real * p,
    Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = (1 - 6*x[i]*y[i])*p[i];
}

void CFN (pInt, Int * n, Int *, Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
  }
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[i] = - xi * xi * xi;
  }
}

void CCFSG (pInt, Int * n, Int *, Real * x, Real * c, Int * nnzJ, Int *,
    Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[i] = - xi * xi * xi;
  }
  if (*Grad == dciFalse)
    return;
  Int k = 0;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    if (xi != 0) {
      J[k] = -3*xi*xi;
      indvar[i] = i;
      indfun[i] = i;
      k++;
    }
  }
  *nnzJ = k;
}

int main () {
  Int n = 10, m = n;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = 1;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  for (Int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf;
    cu[i] = 0;
    equatn[i] = dciFalse;
  }

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
