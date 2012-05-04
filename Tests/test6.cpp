#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem
 * min f(x)
 * s.t. cI(x) <= 0
 *
 * where
 *
 * f(x) = sum(x.^4)
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
    xi = x[i];
    xi *= xi;
    *f += xi*xi;
  }
}

//g(x) = 2*(x - e);
void UGR (Int * n, Real * x, Real * g) {
  Real xi = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    g[i] = 4*xi*xi*xi;
  }
}

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i];
    xi2 = xi*xi;
    *f += xi2*xi2;
    if (*grad == dciTrue)
      g[i] = 4*xi*xi2;
  }
}

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real unusedx = x[0], unusedy = y[0];
  unusedx = x[0];
  unusedy = y[0];
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*n < *m) || (*m != 2) || (*mmax < *m) )
      throw("Error");
    for (Int i = 0; i < *n - 1; i++) {
      Real xi = x[i];
      q[i] = (12*xi*xi + 2*(y[0] + y[1]))*p[i];
    }
    Real xn = x[*n - 1];
    q[*n - 1] = (12*xn*xn + 2*y[0])*p[*n - 1];
  }
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  UFN (n, x, f);
  if ( (*n < *m) || (*m != 2) || (*mmax < *m) )
    throw("Error");
  c[0] = x[*n-1]*x[*n-1] - 1;
  c[1] = -x[*n-1];
  for (Int i = 0; i < *n - 1; i++) {
    Real xi = x[i];
    xi *= xi;
    c[0] += xi;
    c[1] += xi;
  }
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = x[*n-1]*x[*n-1] - 1;
  c[1] = -x[*n-1];
  for (Int i = 0; i < *n - 1; i++) {
    Real xi = x[i];
    xi *= xi;
    c[0] += xi;
    c[1] += xi;
  }
  if ( (*n < *m) || (*m != 2) || (*mmax < *m) || (*jmax < 0) )
    throw("Error");

  if (*Grad == dciFalse)
    return;
  Int k = 0;
  for (Int i = 0; i < *n - 1; i++) {
    Real xi = x[i];
    if (xi != 0) {
      J[k] = 2*xi;
      indvar[k] = i;
      indfun[k] = 0;
      k++;
      J[k] = 2*xi;
      indvar[k] = i;
      indfun[k] = 1;
      k++;
    }
  }
  Real xn = x[*n - 1];
  if (xn != 0) {
    J[k] = 2*xn;
    indvar[k] = *n - 1;
    indfun[k] = 0;
    k++;
  }
  J[k] = -1;
  indvar[k] = *n - 1;
  indfun[k] = 1;
  k++;
  *nnzJ = k;
}

int main () {
  Int n = 3, m = 2;
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
  
  y[0] = 0;
  y[1] = 0;
  cl[0] = 0;
  cl[1] = -dciInf;
  cu[0] = 0;
  cu[1] = 0;
  equatn[0] = dciTrue;
  equatn[1] = dciFalse;

  dci.set_x (n, x);
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
