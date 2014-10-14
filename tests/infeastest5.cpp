#include <iostream>
#include "dci.h"

using namespace DCI;

/* Infeasible problem
 *
 * min f(x)
 * s.t. cI(x) >= 0
 *
 * where
 *
 * f(x) = x(1) + x(2)
 *
 * and
 *
 * cI(x) = [x(1)^2 + x(2)^2 - 1;    <= 0
 *          x(1)^2 + x(2)^2 - 4;  >= 0
 *          x(1)^2 + x(2)^2 - 1/4]    >= 0
 *
 *
 */

void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  Real x1 = x[0], x2 = x[1];
  *f = 7*x1 + 11*x2;
  if (*grad == dciTrue) {
    g[0] = 7;
    g[1] = 11;
  }
}

//H(x,y) = 2*I
void CPROD (pInt, Int *, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  Real y1 = y[0], y2 = y[1], y3 = y[2];
  q[0] = 2*(y1 + y2 + y3)*p[0];
  q[1] = 2*(y1 + y2 + y3)*p[1];
}

void CFN (pInt, Int *, Int *, Real * x, Real * f, Real * c) {
  Real x1 = x[0], x2 = x[1];
  *f = 7*x1 + 11*x2;
  c[0] = x1*x1 + x2*x2 - 1;
  c[1] = c[0] - 3;
  c[2] = c[0] + 3.0/4.0;
}

void CCIFG (pInt, Int *, Int * i, Real * x, Real * ci, Real * gci, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  *ci = x1*x1 + x2*x2 - 1;
  if (*Grad == dciTrue) {
    gci[0] = 2*x1;
    gci[1] = 2*x2;
  }
  if (*i == 2) {
    *ci -= 3;
  } else if (*i == 3) {
    *ci += 3.0/4.0;
  }
}

void CCFSG (pInt, Int *, Int *, Real * x, Real * c, Int * nnzJ, Int *, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x1*x1 + x2*x2 - 1;
  c[1] = c[0] - 3;
  c[2] = c[0] + 3.0/4.0;
  if (*Grad == dciFalse)
    return;

  int k = 0;
  if (x1 != 0) {
    indvar[k] = 0;
    indfun[k] = 0;
    J[k] = 2*x1;
    k++;
    indvar[k] = 0;
    indfun[k] = 1;
    J[k] = 2*x1;
    k++;
    indvar[k] = 0;
    indfun[k] = 2;
    J[k] = 2*x1;
    k++;
  }
  if (x2 != 0) {
    indvar[k] = 1;
    indfun[k] = 0;
    J[k] = 2*x2;
    k++;
    indvar[k] = 1;
    indfun[k] = 1;
    J[k] = 2*x2;
    k++;
    indvar[k] = 1;
    indfun[k] = 2;
    J[k] = 2*x2;
    k++;
  }
  *nnzJ = k;
}

int main () {
  Int n = 2, m = 3;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);
  dci.set_ccifg (CCIFG);

  for (Int i = 0; i < 2; i++) {
    x[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
    y[i] = 0;
    equatn[i] = dciFalse;
  }
  cl[0] = -dciInf;
  cu[0] = 0;
  cl[1] = 0;
  cu[1] = dciInf;
  cl[2] = 0;
  cu[2] = dciInf;
  equatn[2] = dciFalse;


  dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
