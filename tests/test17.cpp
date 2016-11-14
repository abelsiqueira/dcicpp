#include <iostream>
#include "dci.h"
#include <assert.h>
#include <cmath>
extern "C" {
#include "nope.h"
}

using namespace DCI;

/* Let's make the problem
 * min  f(x) = 0.5*norm(x)^2
 * s.t. x_{i+1} = x_i + 1
 *      and x_1 = 1
 *
 * Solution: x = [1, 2, 3, ..., n], f = 0.5*(1 + 2² + 3² + ... + n²) = n(n+1)(2n+1)/12
 *
 */

Nope *nope;
#include "nope_interface.c"

Int nvar = 10;
Int ncon = nvar - 1;

void core_uofg (pInt, const Int *n, const Real * x, Real * f, Real * g, const Bool * grad) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  if (!(*grad))
    return;
  for (Int i = 0; i < *n; i++)
    g[i] = x[i];
}

void core_cofg (pInt, const Int *n, const Real * x, Real * f, Real * g, Bool * grad) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  if (!(*grad))
    return;
  for (Int i = 0; i < *n; i++)
    g[i] = x[i];
}

//H(x,y) = 2*I
void core_uhprod (pInt, const Int *n, const Bool *, const Real *, const Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = p[i];
}

void core_chprod (pInt, const Int *n, const Int *, const Bool *, const Real *, const Real *,
    Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = p[i];
}

void core_ufn (pInt, const Int *n, const Real * x, Real * f) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
}

void core_cfn (pInt, const Int *n, const Int *, const Real * x, Real * f, Real * c) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  for (Int i = 0; i < *n-1; i++)
    c[i] = x[i+1] - x[i] - 1;
}

void core_ccfsg (pInt, const Int *n, const Int *, const Real * x, Real * c, Int * nnzJ, const Int *,
    Real *Jval, Int *Jvar, Int *Jfun, const Bool *grad) {
  for (Int i = 0; i < *n-1; i++) {
    c[i] = x[i+1] - x[i] - 1;
    if (*grad) {
      Jval[2*i] = -1;
      Jval[2*i+1] = 1;
      Jvar[2*i] = i+1;
      Jvar[2*i+1] = i+2;
      Jfun[2*i] = i+1;
      Jfun[2*i+1] = i+1;
    }
  }
  *nnzJ = 2*(*n-1);
}

void core_cdimen(Int *, const Int *, Int *n, Int *m) {
  *n = nvar;
  *m = ncon;
}

void core_csetup(Int *, const Int *, const Int *, const Int *, Int *n, Int *m, Real *x,
    Real *bl, Real *bu, Real *y, Real *cl, Real *cu, Bool *equatn,
    Bool *linear, const Int *, const Int *, const Int *) {
  for (Int i = 0; i < *n; i++) {
    x[i] = 0;
    bl[i] = -dciInf;
    bu[i] = dciInf;
  }
  bl[0] = 1;
  bu[0] = 1;

  for (Int i = 0; i < *m; i++) {
    y[i] = 0;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = dciTrue;
    linear[i] = dciTrue;
  }
}

void core_cdimsj(Int *, Int *nnzj) {
  *nnzj = 2*ncon;
}

void print_array(size_t n, double * x, const char * name) {
  printf("%s = \n", name);
  for (size_t i = 0; i < n; i++) {
    printf(" %8.1e\n", x[i]);
  }
}

int main () {
  Int n = 0, m = 0;
  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, core_ufn, core_uofg, core_uhprod,
      core_csetup, 0, core_cfn, core_cofg, core_chprod, core_ccfsg,
      core_cdimsj);
  runNope(nope);
  ppDIMEN(nope, &n, &m);

  assert(n == 0);
  for (int i = 0; i < nvar; i++)
    assert(nope->x[i] == i+1);
  int st = 0;
  double f = 0.0;
  core_cfn(&st, &nvar, &ncon, nope->x, &f, nope->c);
  assert(fabs(f - nvar*(nvar+1)*(2*nvar+1)/12.0) < 1e-12);
  for (int i = 0; i < ncon; i++)
    assert(nope->c[i] == 0.0);

  destroyNope(nope);

}
