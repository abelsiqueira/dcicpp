#include "dci.h"
#include <cmath>

using namespace DCI;

/* min  f(x)
 * s.t. cI(x) <= 0
 *
 * f(x) = dot(x-e, x-e)
 * cI(x) = sum(x) - 1
 *
 * e = (1, 1, ..., 1)
 * Expected solution:
 *
 * x = (1/n, ... , 1/n)
 * lambda = 2*(n-1)/n;
 * f = (n-1)^2/n
 *
 */

//g(x) = 2*(x - e)
void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += xi*xi;
    if (*grad == dciTrue)
      g[i] = 2*xi;
  }
}

//H(x,y) = 2*I
void CPROD (Int * n, Int * m, Bool * getder, Real * , Int * mmax,
            Real * , Real * p, Real * q) {
  if ( (*getder == 0) || (*getder == 1) ) {
    if ( (*m != 1) || (*mmax < *m) )
      return;
    for (Int i = 0; i < *n; i++)
      q[i] = 2*p[i];
  }
}

void CFN (Int * n, Int * , Real * x, Real * f, Int * , Real * c) {
  Real xi = 0;
  
  *f = 0;
  c[0] = -1;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    *f += xi*xi;
    c[0] += x[i];
  }
}

//J(x) = e'
void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, 
            Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
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
  Int n = 5, m = 1;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  Real sol[n], fsol;

  for (Int i = 0; i < n; i++)
    sol[i] = 1.0/n;
  fsol = (n-1)*(n-1)/(1.0*n);

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
    x[i] = -1;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  
  for (Int i = 0; i < m; i++) {
    y[i] = 0;
    cl[i] = -dciInf; // No lower bound
    cu[i] = 0; // c(x) <= 0
    equatn[i] = dciFalse; // Inequality
  }

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

  pReal px = dci.get_x();
  Real difnorm = 0.0;
  for (Int i = 0; i < n; i++)
    difnorm += pow(px[i] - sol[i], 2);
  std::cout << "|x* - sol|^2 = " << difnorm << std::endl;
  std::cout << "|f* - f| = " << DCI::AbsValue(dci.get_f() - fsol)
    << std::endl;

}
