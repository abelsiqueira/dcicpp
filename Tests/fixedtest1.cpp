#include <iostream>
#include <cmath>
#include "dci.h"

using namespace DCI;

/*
 * min  sum_i x_i^2
 * s.to sum_i x_i = 1
 *      x_1 = 1
 *
 * x0 = (-1, -1, ..., -1)
 * sol = [0, 1, 0]
 * fsol = 1
 *
 *
 */

void COFG (pInt, Int *n, Real * x, Real * f, Real * g, Bool * grad) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    *f += xi*xi;
    if (*grad == dciTrue)
      g[i] = 2*xi;
  }
}

void CPROD (pInt, Int *n, Int *, Bool * , Real * , Real *, Real * p, Real * q) {
  for  (Int i = 0; i < *n; i++)
    q[i] =  2 * p[i];
}

void CFN (pInt, Int *n, Int *, Real * x, Real * f, Real * c) {
  *f = 0.0;
  c[0] = -1;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    *f += xi*xi;
    c[0] += xi;
  }
}

void CCFSG (pInt, Int *n, Int *, Real * x, Real * c, Int * nnzJ, Int *,
      Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = -1;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    c[0] += xi;
    if (*Grad) {
      indvar[i] = i;
      indfun[i] = 0;
      J[i] = 1;
    }
  }
  *nnzJ = *n;
}

int main () {
  Int n = 100000, m = 1;
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
  bl[1] = 1;
  bu[1] = 1;

  cl[0] = 0;
  cu[0] = 0;
  equatn[0] = dciTrue;

  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
