#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/* Let's make the problem ARGLINB
 * min sum_{i = 1}^M 
 *      ( sum_{j = 1}^N   x_j * i * j - 1.0 )^2
 *     
 * with N = 10 and M = 20
 *
 */

int m;

void UOFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  *f = 0;
  if (*grad == dciTrue) {
    for (Int i = 0; i < *n; i++)
      g[i] = 0.0;
  }
  for (Int i = 0; i < m; i++) {
    Real xjj = 0.0;
    for (Int j = 0; j < *n; j++) {
      xjj += x[j]*(j+1);
    }
    xjj = (i+1)*xjj - 1;
    *f += xjj*xjj;
    if (*grad == dciTrue) {
      for (Int k = 0; k < *n; k++)
        g[k] += 2*xjj*(i+1)*(k+1);
    }
  }
}

void UPROD (Int * n, Bool *, Real *, Real * p, Real * q) {
  Real sum_pjj = 0.0;
  for (Int j = 0; j < *n; j++) {
    sum_pjj += p[j]*(j+1);
  }
  Real aux = m*(m+1)*(2*m+1)/3;
  for (Int i = 0; i < *n; i++) {
    q[i] = aux*(i+1)*sum_pjj;
  }
}

void UFN (Int * n, Real * x, Real * f) {
  *f = 0;
  for (Int i = 0; i < m; i++) {
    Real xjj = 0.0;
    for (Int j = 0; j < *n; j++) {
      xjj += x[j]*(j+1);
    }
    xjj = (i+1)*xjj - 1;
    *f += xjj*xjj;
  }
}

int main () {
  m = 4;
  Int n = 2;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];

  dci.set_ufn (UFN);
  dci.set_uofg (UOFG);
  dci.set_uprod (UPROD);

  for (Int i = 0; i < n; i++) {
    x[i] = 1.0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }

  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);

  dci.start ();
  dci.solve ();
  dci.show();

}
