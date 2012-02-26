#include <iostream>
#include "dci.h"
extern "C" {
#include "cuter.h"
}

using namespace DCI;

int MAINENTRY () {
  DCI::Interface dci;
  char fname[10] = "OUTSDIF.d";
  Int nvar = 0, ncon = 0, amax = 0;
  Int nmax, mmax;
  Int funit = 42, ierr = 0, fout = 6;

  FORTRAN_OPEN ((&funit), fname, (&ierr));
  CDIMEN ((&funit), (&nvar), (&ncon));

  dci.Cuter ();

  Real x[nvar], bl[nvar], bu[nvar];
  nmax = nvar;

  if (ncon == 0) {
    USETUP ((&funit), (&fout), (&nvar), x, bl, bu, (&nmax));

    dci.set_uofg (UOFG);
    dci.set_uprod (UPROD);
    dci.set_ufn (UFN);
    dci.set_unames (UNAMES);

  } else {
    mmax = ncon;
    Real y[ncon], cl[ncon], cu[ncon];
    Bool equatn[ncon], linear[ncon];
    Bool efirst = 0, lfirst = 0, nvfirst = 0;

    CSETUP ((&funit), (&fout), (&nvar), (&ncon), x, bl, bu, (&nmax), equatn, linear, y, cl, cu, (&mmax), (&efirst), (&lfirst), (&nvfirst));
    CDIMSJ ((&amax));
    amax *= 2;
  
    dci.set_cofg (COFG);
    dci.set_cprod (CPROD);
    dci.set_cfn (CFN);
    dci.set_ccfsg (CCFSG);
    dci.set_ccifg (CCIFG);
    dci.set_cnames (CNAMES);

    dci.set_lambda (ncon, y);
    dci.set_cl (ncon, cl);
    dci.set_cu (ncon, cu);
    dci.set_equatn (ncon, equatn);
    dci.set_amax (amax);
  }

  dci.set_x (nvar, x);
  dci.set_bl (nvar, bl);
  dci.set_bu (nvar, bu);

  dci.start ();
  try {
    dci.solve ();
  } catch (...) {
  }
  dci.show (1);
  dci.printLatex ();

  FORTRAN_CLOSE ((&funit), (&ierr));

  return 0;

}
