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
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
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
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*n != 2) || (*m != 1) || (*mmax < *m) )
      return;
    q[0] = (1 - 6*y[0]*x[0])*p[0];
    q[1] = p[1];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += 0.5*xi*xi;
  }
  if (*m != 1)
    return;
  if (*mmax < 1)
    return;
  c[0] = x[1] - x[0]*x[0]*x[0] + x[0];
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = x[1] - x[0]*x[0]*x[0] + x[0];
  if (*Grad == dciFalse)
    return;
  if (*n != 2)
    return;
  if (*m != 1)
    return;
  if (*mmax != 1)
    return;
  if (*jmax < 0)
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
  Real x[n], bl[n], bu[n], sol[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = -i;
    sol[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  for (int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf;
    cu[i] = 0;
    equatn[i] = dciFalse;
  }

  dci.set_x (n, x);
  dci.set_sol (n, sol);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);
  dci.set_lambda (m, y);
  dci.set_cl (m, cl);
  dci.set_cu (m, cu);
  dci.set_equatn (m, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
