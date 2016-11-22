#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) <= 0
 *
 * where
 *
 * f(x) = 0.5*dot(x-1,x-1)
 *
 * and
 *
 * cI(x) = x(2) - x(1)^3 + x(1)
 *
 * One Solution:
 * sol = [0;0]
 * fsol = 1
 *
 * Global Solution:
 * sol ~= [1.30697; 0.92556]
 * fsol ~= 0.049989
 *
 */

void COFG (pInt, const Int * n, const Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += 0.5*xi*xi;
    if (*grad == dciTrue)
      g[i] = xi;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, const Int *, const Int *, const Bool *, const Real * x, const Real * y, Real * p,
    Real * q) {
  q[0] = (1 - 6*y[0]*x[0])*p[0];
  q[1] = p[1];
}

void CFN (pInt, const Int * n, const Int *, const Real * x, Real * f, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += 0.5*xi*xi;
  }
  c[0] = x[1] - x[0]*x[0]*x[0] + x[0];
}

void CCFSG (pInt, const Int *, const Int *, const Real * x, Real * c, Int * nnzJ, const Int *,
    Real * J, Int * indvar, Int * indfun, const Bool * Grad) {
  c[0] = x[1] - x[0]*x[0]*x[0] + x[0];
  if (*Grad == dciFalse)
    return;
  Real val = 1 - 3*x[0]*x[0];
  Int k = 0;
  if (val != 0) {
    J[0] = val;
    indvar[0] = 0;
    indfun[0] = 0;
    k++;
  }
  J[k] = 1;
  indvar[k] = 1;
  indfun[k] = 0;

  *nnzJ = k + 1;
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
    x[i] = -i;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (int i = 0; i < m; i++) {
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
