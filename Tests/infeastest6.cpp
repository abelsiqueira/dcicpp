#include <iostream>
#include "dci.h"
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace DCI;

/* Infeasible problem
 *
 * min f(x) = dot(c,x)
 * s.t. cI(x) >= 0
 *
 * where
 *
 * f(x) = c(1)*x(1) + c(2)*x(2)
 *
 * and
 *
 * cE(x) = [dot(x,x) - 1;        <= 0
 *          dot(x,x) - 4]        >= 0
 *
 * sol = -sqrt(5/2) * c/norm(c)
 *
 */

Int DIM = 10;
Real v[100];

Real rand_between (Real a, Real b) {
  Real x = (rand()%1001/1000.0);
  return (b - a)*x + a;
}

void COFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += v[i] * x[i];
  if (*grad == dciTrue) {
    for (Int i = 0; i < *n; i++)
      g[i] = v[i];
  }
}

void CPROD (Int * n, Int * m, Bool * getder, Real * x, Int * mmax, Real * y, Real * p, Real * q) {
  Real y1 = y[0], y2 = y[1];
  for (Int i = 0; i < *n; i++)
    q[i] = 2*(y1 + y2)*p[i];
}

void CFN (Int * n, Int * m, Real * x, Real * f, Int * mmax, Real * c) {
  *f = 0.0;
  c[0] = -1.0;
  for (Int i = 0; i < *n; i++) {
    *f += v[i] * x[i];
    c[0] += x[i] * x[i];
  }
  c[1] = c[0] - 3;
}

void CCFSG (Int * n, Int * m, Real * x, Int * mmax, Real * c, Int * nnzJ, Int * jmax, Real * J, Int * indvar, Int * indfun, Bool * Grad) {
  c[0] = -1.0;
  for (Int i = 0; i < *n; i++) {
    c[0] += x[i] * x[i];
  }
  c[1] = c[0] - 3;
  if (*Grad == dciFalse)
    return;

  int k = 0;
  for (Int i = 0; i < *n; i++) {
    Real xi = x[i];
    if (xi != 0) {
      indvar[k] = i;
      indfun[k] = 0;
      J[k] = 2*xi;
      k++;
      indvar[k] = i;
      indfun[k] = 1;
      J[k] = 2*xi;
      k++;
    }
  }
  *nnzJ = k;
}

int main () {
  Int n = DIM, m = 2;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];
  Real y[m], cl[m], cu[m];
  Bool equatn[m];

  srand(time(0));
  for (Int i = 0; i < n; i++)
    v[i] = rand_between(-1, 1);

  dci.set_cofg (COFG);
  dci.set_cprod (CPROD);
  dci.set_cfn (CFN);
  dci.set_ccfsg (CCFSG);

  for (Int i = 0; i < n; i++) {
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
  Real vnorm = 0.0;
  Real residue = 0.0;

  for (Int i = 0; i < n; i++)
    vnorm += v[i]*v[i];
  vnorm = sqrt(vnorm);

  Real alpha = -sqrt(5.0/2.0)/vnorm;

  for (Int i = 0; i < n; i++)
    residue += pow( px[i] - alpha*v[i], 2 );

/*   if (residue > DIM*1e-6) {
 *     std::cout << "x = " << std::endl;
 *     for (Int i = 0; i < n; i++)
 *       std::cout << px[i] << ' ';
 *     std::cout << std::endl;
 *     std::cout << "c = " << std::endl;
 *     for (Int i = 0; i < n; i++)
 *       std::cout << v[i] << ' ';
 *     std::cout << std::endl;
 *   }
 */

  Real dot_sol = -sqrt(5.0/2.0)*vnorm, dot_x = 0.0;

  for (Int i = 0; i < n; i++) {
    dot_x += v[i] * px[i];
  }

  std::cout << "dot(sol,c) = " << dot_sol << std::endl;
  std::cout << "dot(x,c) = " << dot_x << std::endl;
  std::cout << "|sol - x| = " << residue << std::endl;

}
