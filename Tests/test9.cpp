#include <iostream>
#include "dci.h"

using namespace DCI;

/* Let's make the problem HS18
 * min f(x)
 * s.t. l <= cI(x) <= u
 *
 * where
 *
 * f(x) = 0.01*x(1)^2 + x(2)^2
 *
 * and
 *
 * cI(x) = [x(1) * x(2) - 25; 
 *          x(1)^2 + x(2)^2 - 25;
 *          x(1);
 *          x(2)]
 *
 * l = [  0,   0,  2,  0]
 * u = [inf, inf, 50, 50]
 *
 * x = [sqrt(10)*5, sqrt(10)/2] ~ [15.811, 1.5811]
 * f(x) = 5
 *
 */

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  if (*n != 2)
    return;
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  if (*grad == dciTrue) {
    g[0] = 0.02 * x1;
    g[1] = 2 * x2;
  }
}

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * , Real * , Int * mmax, Real * y, Real * p, Real * q) {
  if ( (*n != 2) || (*m != 4) || (*mmax < *m) )
    return;
  Real y1 = y[0], y2 = y[1], y3 = y[2], y4 = y[3];
  q[0] = (0.02 + 2*y2 + y3) * p[0] + y1 * p[1];
  q[1] = y1 * p[0] + (2 + 2*y2 + y4) * p[1];
}

void CFN (Int * , Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  if ( (*m != 4) || (*mmax < *m) )
    throw("Error");
  Real x1 = x[0], x2 = x[1];
  *f = 0.01 * x1 * x1 + x2 * x2;
  c[0] = x1 * x2 - 25;
  c[1] = x1 * x1 + x2 * x2 - 25;
  c[2] = x1;
  c[3] = x2;
}

void CCFSG (Int * , Int * , Real * x, Int * , Real * c, Int * nnzJ, Int * , Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  Real x1 = x[0], x2 = x[1];
  c[0] = x1 * x2 - 25;
  c[1] = x1 * x1 + x2 * x2 - 25;
  c[2] = x1;
  c[3] = x2;
  if (*Grad == dciFalse)
    return;
  /*
   * J = [  x(2),   x(1);
   *      2*x(1), 2*x(2);
   *           1,      0;
   *           0,      1]
   */
  if (x1 == 0) {
    if (x2 == 0) {
      *nnzJ = 0;
    } else {
      *nnzJ = 2;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
      indvar[1] = 1;
      indfun[1] = 1;
      J[1] = 2*x2;
    }
  } else {
    if (x2 == 0) {
      *nnzJ = 2;
      indvar[0] = 0;
      indfun[0] = 1;
      J[0] = 2*x1;
      indvar[1] = 1;
      indfun[1] = 0;
      J[1] = x1;
    } else {
      *nnzJ = 4;
      indvar[0] = 0;
      indfun[0] = 0;
      J[0] = x2;
      indvar[1] = 0;
      indfun[1] = 1;
      J[1] = 2*x1;
      indvar[2] = 1;
      indfun[2] = 0;
      J[2] = x1;
      indvar[3] = 1;
      indfun[3] = 1;
      J[3] = 2*x2;
    }
  }
  int k = *nnzJ;
  indvar[k] = 0;
  indfun[k] = 2;
  J[k] = 1;
  k++;
  indvar[k] = 1;
  indfun[k] = 3;
  J[k] = 1;
  *nnzJ += 2;
}

int main () {
  Int n = 2, m = 4;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  x[0] = 2;
  x[1] = 2;
  for (Int i = 0; i < 2; i++) {
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  y[0] = 0;
  y[1] = 0;
  y[2] = 0;
  y[3] = 0;
  cl[0] = 0;
  cl[1] = 0;
  cl[2] = 2;
  cl[3] = 0;
  cu[0] = dciInf;
  cu[1] = dciInf;
  cu[2] = 50;
  cu[3] = 50;
  equatn[0] = dciFalse;
  equatn[1] = dciFalse;
  equatn[2] = dciFalse;
  equatn[3] = dciFalse;

  dci.con_setup (n, x, bl, bu, m, y, cl, cu, equatn);

  dci.start ();
  dci.solve ();
  dci.show();

}
