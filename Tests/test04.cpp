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
 * f(x) = - 0.5*dot(x,x)
 *
 * and
 *
 * cI_1(x) = x(1)^2 + x(2)^2/r^2 - 1
 * cI_2(x) = x(1)^2/r^2 + x(2)^2 - 1
 *
 * Expected solution:
 *
 * s = r/sqrt(r^2 + 1)
 * x = [+- s; +- s] (4 solutions)
 * f = -s^2 = -r^2/(r^2 + 1)
 *
 */

Real r = 2;

void COFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f -= 0.5*xi*xi;
    if (*grad == dciTrue)
      g[i] = -xi;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  q[0] = (-1 + 2 * y[0] + 2 * y[1] / (r*r) ) * p[0];
  q[1] = (-1 + 2 * y[0] / (r*r) + 2 * y[1] ) * p[1];
}

void CFN (pInt, Int *n, Int *, Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f -= 0.5*xi*xi;
  }
  Real v1 = x[0]*x[0], v2 = x[1]*x[1], rr = r*r;
  c[0] = v1 + v2/rr - 1;
  c[1] = v1/rr + v2 - 1;
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J,
    Int * indvar, Int * indfun, Bool * Grad) {
  Real v1 = x[0]*x[0], v2 = x[1]*x[1], rr = r*r;
  c[0] = v1 + v2/rr - 1;
  c[1] = v1/rr + v2 - 1;
  if (*Grad == dciFalse)
    return;
  Int k = 0;
  v1 = x[0];
  v2 = x[1];
  if (v1 != 0) {
    J[k] = 2*v1;
    indfun[k] = 0;
    indvar[k] = 0;
    k++;
    J[k] = 2*v1/rr;
    indfun[k] = 1;
    indvar[k] = 0;
    k++;
  }
  if (v2 != 0) {
    J[k] = 2*v2;
    indfun[k] = 1;
    indvar[k] = 1;
    k++;
    J[k] = 2*v2/rr;
    indfun[k] = 0;
    indvar[k] = 1;
    k++;
  }
  *nnzJ = k;
}

int main () {
  Int n = 2, m = 2;
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
  
  for (int i = 0; i < m; i++) {
    y[i] = 1;
    cl[i] = -dciInf;
    cu[i] = 0;
    equatn[i] = dciFalse;
  }

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
