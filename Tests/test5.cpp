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

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
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
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != *n) || (*mmax < *m) )
      throw("Error");
    for (Int i = 0; i < *n; i++)
      q[i] = (1 - 6*x[i]*y[i])*p[i];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    *f += xi*xi;
  }
  if (*m != *n)
    throw("Error");
  if (*mmax < 1)
    throw("Error");
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[i] = - xi * xi * xi;
  }
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[i] = - xi * xi * xi;
  }
  if (*Grad == dciFalse)
    return;
  if (*m != *n)
    throw("Error");
  if (*mmax != *n)
    throw("Error");
  if (*jmax < 0)
    throw("Error");
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
  Real x[n], bl[n], bu[n], sol[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = 1;
    sol[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  for (Int i = 0; i < m; i++) {
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
