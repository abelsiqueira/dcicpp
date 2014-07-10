#include <iostream>
#include "dci.h"
extern "C" {
#include "preprocessor.h"
}

using namespace DCI;

/* Let's make the problem
 * min  f(x) = 0.5*norm(x)^2
 * s.t. x_{i+1} = x_i
 *      and x_1 = 1
 *
 * Solution: x = ones(n,1), f = n/2
 *
 */

Preprocessor *prep;
#include "prep_interface.c"

Int nvar = 10;
Int ncon = nvar - 1;

void core_uofg (pInt, Int *n, Real * x, Real * f, Real * g, Bool * grad) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  if (!(*grad))
    return;
  for (Int i = 0; i < *n; i++)
    g[i] = x[i];
}

void core_cofg (pInt, Int *n, Real * x, Real * f, Real * g, Bool * grad) {
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
void core_uhprod (pInt, Int *n, Bool *, Real *, Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = p[i];
}

void core_chprod (pInt, Int *n, Int *, Bool *, Real *, Real *,
    Real * p, Real * q) {
  for (Int i = 0; i < *n; i++)
    q[i] = p[i];
}

void core_ufn (pInt, Int *n, Real * x, Real * f) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
}

void core_cfn (pInt, Int *n, Int *, Real * x, Real * f, Real * c) {
  *f = 0.0;
  for (Int i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  for (Int i = 0; i < *n-1; i++)
    c[i] = x[i+1] - x[i];
}

void core_ccfsg (pInt, Int *n, Int *, Real * x, Real * c, Int * nnzJ, Int *,
    Real *Jval, Int *Jvar, Int *Jfun, Bool *grad) {
  for (Int i = 0; i < *n-1; i++) {
    c[i] = x[i+1] - x[i];
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

void core_cdimen(Int *, Int *, Int *n, Int *m) {
  *n = nvar;
  *m = ncon;
}

void core_csetup(Int *, Int *, Int *, Int *, Int *n, Int *m, Real *x,
    Real *bl, Real *bu, Real *y, Real *cl, Real *cu, Bool *equatn,
    Bool *linear, Int *, Int *, Int *) {
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

int main () {
  Int n = 0, m = 0;
  DCI::Interface dci;

  dci.cuter();
  prep = initializePreprocessor();

  setFuncs(prep, core_cdimen, 0, core_ufn, core_uofg, core_uhprod,
      core_csetup, core_cfn, core_cofg, core_chprod, core_ccfsg,
      core_cdimsj);
  runPreprocessor(prep);
  ppDIMEN(prep, &n, &m);

  Real x[n], bl[n], bu[n];

  if (m == 0) {
    dci.set_ufn(ufn);
    dci.set_uofg(uofg);
    dci.set_uprod(uhprod);

    runUncSetup(prep, &n, x, bl, bu);
    dci.unc_setup(n, x, bl, bu);
  } else {
    Real y[m], cl[m], cu[m];
    Bool equatn[m], linear[m];
    Int jmax;

    dci.set_cfn(cfn);
    dci.set_cofg(cofg);
    dci.set_cprod(chprod);
    dci.set_ccfsg(ccfsg);

    runConSetup(prep, &n, x, bl, bu, &m, y, cl, cu, equatn, linear, &jmax);
    dci.set_linear(m, linear);
    dci.set_amax(jmax);
    dci.con_setup(n, x, bl, bu, m, y, cl, cu, equatn);
  }

  dci.start ();
  dci.solve ();
  dci.show();

  destroyPreprocessor(prep);

}
