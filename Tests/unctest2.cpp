#include "dci.h"
#include <cmath>

using namespace DCI;

/* min  f(x)
 *
 * f(x) = x1 * exp(x2*ti) - yi
 *
 *
 */

const Int m = 2;
Real t[m], y[m];

void generate_t_and_y () {
  for (Int i = 0; i < m; i++) {
    t[i] = i + 5;
    y[i] = 0.001 * exp (2 * t[i]);
  }
}

Real hi (Real * x, Int i) {
  return x[0] * exp( x[1] * t[i] ) - y[i];
}

void Ji (Real * x, Int i, Real * g) {
  Real etix2 = exp(x[1] * t[i]);
  g[0] = etix2;
  g[1] = etix2 * x[0] * t[i];
}

void Hiprod (Real * x, Int i, Real * p, Real * q) {
  Real etix2 = exp(x[1] * t[i]);
  q[0] = etix2 * t[i] * p[1];
  q[1] = etix2 * t[i] * ( p[0] + t[i] * x[0] * p[1] );
}

void UOFG (Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < m; i++)
    *f += 0.5*pow( hi(x, i), 2);
    
  if (*grad == dciTrue) {
    Real aux[2], auxhi;
    g[0] = 0.0;
    g[1] = 0.0;
    for (Int i = 0; i < m; i++) {
      Ji(x, i, aux);
      auxhi = hi(x, i);
      g[0] += auxhi * aux[0];
      g[1] += auxhi * aux[1];
    }
  }
}

void UPROD (Int * n, Bool * getder, Real * x, Real * p, Real * q) {
  Real xi = 0;
  q[0] = 0.0;
  q[1] = 0.0;
  Real aux[2], auxdot, auxhi;
  for (Int i = 0; i < m; i++) {
    Ji(x, i, aux);
    auxdot = aux[0] * p[0] + aux[1] * p[1];
    q[0] += auxdot * aux[0];
    q[1] += auxdot * aux[1];
    Hiprod(x, i, p, aux);
    auxhi = hi(x, i);
    q[0] += auxhi * aux[0];
    q[1] += auxhi * aux[1];
  }
}

void UFN (Int * n, Real * x, Real * f) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < m; i++)
    *f += 0.5 * pow( hi(x, i), 2);
}

int main () {
  Int n = 2;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];

  generate_t_and_y();

  dci.set_uofg (UOFG);
  dci.set_uprod (UPROD);
  dci.set_ufn (UFN);

  for (Int i = 0; i < n; i++) {
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  x[0] = 0.001;
  x[1] = 2.1;
  
  dci.set_x (n, x);
  dci.set_bl (n, bl);
  dci.set_bu (n, bu);

  dci.start ();
  dci.solve ();
  dci.show();

}
