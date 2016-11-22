#include "dci.h"
#include <cmath>

using namespace DCI;

/* min  f(x)
 * s.t. cI(x) <= 0
 *
 * f(x) = dot(x-e, x-e)
 * cI(x) = sum(x) - 1
 *
 * e = (1, 1, ..., 1)
 * Expected solution:
 *
 * x = (1/n, ... , 1/n)
 * lambda = 2*(n-1)/n;
 * f = (n-1)^2/n
 *
 */

//g(x) = 2*(x - e)
void COFG (pInt, const Int * n, const Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += xi*xi;
    if (*grad == dciTrue)
      g[i] = 2*xi;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, const Int * n, const Int *, const Bool *, const Real *, const Real * , Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = 2*p[i];
}

void CFN (pInt, const Int * n, const Int *, const Real * x, Real * f, Real * c) {
  Real xi = 0;

  *f = 0;
  c[0] = -1;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += xi*xi;
    c[0] += x[i];
  }
}

//J(x) = e'
void CCFSG (pInt, const Int * n, const Int *, const Real * x, Real * c, Int * nnzJ, const Int *,
    Real * J, Int * indvar, Int * indfun, const Bool * Grad) {
  c[0] = -1;
  for (Int i = 0; i < *n; i++)
    c[0] += x[i];
  if (*Grad == dciFalse)
    return;
  for (Int i = 0; i < *n; i++) {
    J[i] = 1;
    indvar[i] = i;
    indfun[i] = 0;
  }
  *nnzJ = *n;
}

int main () {
  Int n = 5, m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = -1;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (Int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf; // No lower bound
    cu[i] = 0; // c(x) <= 0
    equatn[i] = dciFalse; // Inequality
  }

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
