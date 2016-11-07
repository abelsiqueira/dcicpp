#include "dci.h"
#include <cmath>

using namespace DCI;

const Real a = -1.0;
const Real b = 0.5;

//f(x) = x1
//g(x) = [1.0; 0.0; 0.0]
void COFG (pInt, Int *, Real * x, Real * f, Real * g, Bool * grad) {
  *f = x[0];
  if (*grad == dciTrue) {
    g[0] = 1.0;
    g[1] = 0.0;
    g[2] = 0.0;
  }
}

//H(x,y) = [2y_1 0 0;0 0 0;0 0 0] 
void CPROD (pInt, Int * n, Int *, Bool *, Real *, Real * y, Real * p, Real * q) {
  q[0] = 2*y[0]*p[0];
  q[1] = q[2] = 0.0;
}

// c(x) = [x_1^2 - x_2 + a;
//         x_1 - x_3 - b]
void CFN (pInt, Int * n, Int *, Real * x, Real * f, Real * c) {
  Real xi = 0;

  *f = x[0];
  c[0] = x[0]*x[0] - x[1] + a;
  c[1] = x[0] - x[2] - b;
}

//J(x) = [2x_1  -1.0   0.0;
//         1.0   0.0  -1.0]
void CCFSG (pInt, Int * n, Int *, Real * x, Real * c, Int * nnzJ, Int *,
    Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = x[0]*x[0] - x[1] + a;
  c[1] = x[0] - x[2] - b;
  if (*Grad == dciFalse)
    return;

  J[0] = 2*x[0];
  J[1] = -1.0;
  J[2] = 1.0;
  J[3] = -1.0;
  indvar[0] = 0;
  indvar[1] = 1;
  indvar[2] = 0;
  indvar[3] = 2;
  indfun[0] = 0;
  indfun[1] = 0;
  indfun[2] = 1;
  indfun[3] = 1;
  *nnzJ = 4;
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

  x[0] = -3.0;
  x[1] = 1.0;
  x[2] = 1.0;
  bl[0] = -dciInf;
  bl[1] = 0.0;
  bl[2] = 0.0;
  bu[0] = dciInf;
  bu[1] = dciInf;
  bu[2] = dciInf;

  for (Int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = dciTrue;
  }

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
