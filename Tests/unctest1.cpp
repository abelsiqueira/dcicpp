#include "dci.h"

using namespace DCI;

/* min  f(x)
 *
 * f(x) = (x[i] - 1)^4
 *
 *
 */

//g(x) = 2*(x - e)
void UOFG (pInt, Int * n, Real * x, Real * f, Real * g, Bool * grad) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    xi2 = xi*xi;
    *f += xi2*xi2;
    if (*grad == dciTrue)
      g[i] = 4*xi*xi2;
  }
}

//H(x,y) = 2*I
void UPROD (pInt, Int * n, Bool *, Real * x, Real * p, Real * q) {
  Real xi = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    q[i] = 12*xi*xi*p[i];
  }
}

void UFN (pInt, Int * n, Real * x, Real * f) {
  Real xi = 0, xi2 = 0;
  *f = 0;
  for (Int i = 0; i < *n; i++) {
    xi = x[i] - 1;
    xi2 = xi*xi;
    *f += xi2*xi2;
  }
}

int main () {
  Int n = 5;
  DCI::Interface dci;
  Real x[n], bl[n], bu[n];

  dci.set_uofg (UOFG);
  dci.set_uprod (UPROD);
  dci.set_ufn (UFN);

  for (Int i = 0; i < n; i++) {
    x[i] = -1;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }


  dci.unc_setup(n, x, bl, bu);

  dci.start ();
  dci.solve ();
  dci.show();

}
