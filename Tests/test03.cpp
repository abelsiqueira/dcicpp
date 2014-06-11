#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) <= 0
 *
 * where
 *
 * f(x) = 0.5*dot(x,x)
 *
 * and
 *
 * cI(x) = r^2 - dot(x,Dx),
 *
 * where D = diag(1,2,3,...,n).
 *
 * Global solution:
 *
 * f = 0.5*r^2/n
 * x = [0, 0, 0, ..., r/sqrt(n)]
 * lambda = r/(n*sqrt(n))
 *
 * All local solutions:
 *
 * x = +- r * k^{-1/2} * e_k, k = 1, ..., n
 * f = r^2/(2*k)
 *
 */

Real r = 3;

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += 0.5*xi*xi;
    if (*grad == dciTrue)
      g[i] = xi;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int * n, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = (1 - 2*y[0]*(i+1))*p[i];
}

void CFN (pInt, Int * n, Int *, Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += 0.5*xi*xi;
  }
  c[0] = r*r;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[0] -= xi*xi*(i+1);
  }
}

void CCFSG (pInt, Int * n, Int *, Real * x, Real * c, Int * nnzJ, Int *,
    Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = r*r;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[0] -= xi*xi*(i+1);
  }
  if (*Grad == dciFalse)
    return;
  Int k = 0;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    if (xi != 0) {
      J[k] = -2*(i+1)*xi;
      indvar[k] = i;
      indfun[k] = 0;
      k++;
    }
  }
  *nnzJ = k;
}

int main () {
  Int n = 10, m = 1;
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
