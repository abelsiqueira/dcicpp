#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) <= 0
 *
 * where
 *
 * f(x) = dot(x-1,x-1)
 *
 * and
 *
 * cI(x) = sum(x) - 1
 *
 * Expected solution:
 *
 * x = [1/n, ... , 1/n]
 * lambda = 2*(n-1)/n;
 *
 */

void UFN (Int * n, Real * x, Real * f) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += xi*xi;
  }
}

//g(x) = 2*(x - e);
void UGR (Int * n, Real * x, Real * g) {
  Real xi = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    g[i] = 2*xi;
  }
}

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  UFN (n, x, f);
  if (*grad == 1)
    UGR (n, x, g);
}

/*void UPROD (Int * n, Bool * getder, Real * x, Real * p, Real * q) {
  Real unused = x[0];
  unused = x[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    for (Int i = 0; i < *n; i++) {
      q[i] = 2 * (i + 1) * p[i];
    }
  }
}*/

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real unusedx = x[0], unusedy = y[0];
  unusedx = x[0];
  unusedy = y[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != 1) || (*mmax < *m) )
      return;
    for (Int i = 0; i < *n; i++)
      q[i] = 2*p[i];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  UFN (n, x, f);
  if (*m != 1)
    return;
  if (*mmax < 1)
    return;
  c[0] = -1;
  for (Int i = 0; i < *n; i++)
    c[0] += x[i];
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = -1;
  for (Int i = 0; i < *n; i++)
    c[0] += x[i];
  if (*Grad == dciFalse)
    return;
  if (*m != 1)
    return;
  if (*mmax != 1)
    return;
  if (*jmax < 0)
    return;
  for (Int i = 0; i < *n; i++) {
    J[i] = 1;
    indvar[i] = i;
    indfun[i] = 0;
  }
  *nnzJ = *n;
}

int main () {
  Int n = 50, m = 1;
//  Int nmax = 1000, mmax = 1000;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n], sol[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

//  dci.set_ufn (UFN);
//  dci.set_ugr (UGR);
  dci.set_cofg (COFG);
//  dci.set_uprod (UPROD);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = -1;
    sol[i] = 1.0/n;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  for (Int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf;
    cu[i] = 0; // c(x) <= 0
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
