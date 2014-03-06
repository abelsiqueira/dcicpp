#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) >= 0
 *      cE(x) = 0
 *
 * where
 *
 * f(x) = 0.5*( x(1)^2  + (x(2)+1)^2 )
 *
 * and
 *
 * cI(x) = x(2) - x(1)^2
 * cE(x) = x(1) + x(2) - 1
 *
 */

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1]+1;
  *f = 0.5*(x1*x1 + x2*x2);
  if (*grad == dciTrue) {
    g[0] = 2*x1;
    g[1] = x2;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y,
    Real * p, Real * q) {
  q[0] = (1 - 2*y[0])*p[0];
  q[1] = p[1];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1]+1;
  *f = 0.5*(x1*x1 + x2*x2);
  x2 = x[1];
  c[0] = x2 - x1*x1;
  c[1] = x1 + x2 - 1;
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J,
    Int * indvar, Int * indfun, Bool *) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x2 - x1*x1;
  c[1] = x1 + x2 - 1;
  J[0] = -2*x1;
  indvar[0] = 0;
  indfun[0] = 0;
  J[1] = 1;
  indvar[1] = 1;
  indfun[1] = 0;
  J[2] = 1;
  indvar[2] = 0;
  indfun[2] = 1;
  J[3] = 1;
  indvar[3] = 1;
  indfun[3] = 1;
  *nnzJ = 4;
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

  x[0] = 1;
  x[1] = -1;
  for (Int i = 0; i < n; i++) {
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  for (int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = 0;
  }
  cu[0] = dciInf;
  cu[1] = 0;
  equatn[0] = dciFalse;
  equatn[1] = dciTrue;

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
