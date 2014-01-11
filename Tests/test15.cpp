#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* 
 * min  
 * s.t.
 *    
 *   
 *  
 * 
 *   
 *  
 * 
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  *f = -(*n);
  for (Int i = 0; i < *n; i++) {
    Real xi = (x[i] - 1);
    xi = xi*xi;
    xi = xi*xi;
    *f += xi;
  }
  if (*grad == dciTrue) {
    for (Int i = 0; i < *n; i++) {
      Real xi = (x[i] - 1);
      xi = xi*xi*xi;
      g[i] = 4*xi;
    }
  }
}

void CPROD (Int * n, Int *, Bool * , Real * x, Int * , Real * y, Real * p, Real * q) {
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i] - 1;
    xi = xi*xi;
    q[i] = 12*p[i]*xi + 2*y[i];
  }
}

void CFN (Int * n, Int * , Real * x, Real * f, Int * , Real * c) {
  *f = -(*n);
  c[0] = -(*n);
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i] - 1;
    xi = xi*xi;
    *f += xi;
    c[0] += (x[i] + 1)*(x[i] + 1);
  }
}

void CCFSG (Int * n, Int * , Real * x, Int * , Real * c, Int * nnzJ, Int * , Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = -(*n);
  for (Int i = 0; i < *n; i++)
    c[0] += (x[i] + 1)*(x[i] + 1);
  if (*Grad == dciFalse)
    return;
  for (Int i = 0; i < *n; i++) {
    indvar[i] = i;
    indfun[i] = 0;
    J[i] = 2*(x[i] + 1);
  }
  *nnzJ = *n;
}

int main () {
  Int n = 1000, m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m], linear[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  Real a = 100;
  for (Int i = 0; i < n; i++)
    x[i] = -1;

  for (Int i = 0; i < n; i++) {
    bl[i] = -a;
    bu[i] = a;
  }

  y[0] = 0;
  cl[0] = 0;
  cu[0] = 0;
  equatn[0] = dciTrue;
  linear[0] = dciTrue;

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);
  dci.set_linear (m, linear);

  dci.start ();
  dci.solve ();
  dci.show();

}
