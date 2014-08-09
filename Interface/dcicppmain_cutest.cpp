#include <iostream>
#include "dci.h"
extern "C" {
#include "nope.h"
#include "cutest.h"
}

using namespace DCI;

Nope *nope;
#include "nope_interface.c"

int MAINENTRY () {
  DCI::Interface dci;
  char fname[10] = "OUTSDIF.d";
  Int nvar = 0, ncon = 0, amax = 0;
  Int funit = 42, ierr = 0, status;

  nope = initializeNope();

  FORTRAN_open (&funit, fname, &ierr);
  setFuncs(nope, CUTEST_cdimen, CUTEST_usetup, CUTEST_ufn,
      CUTEST_uofg, CUTEST_uhprod, CUTEST_csetup, CUTEST_cfn,
      CUTEST_cofg, CUTEST_chprod, CUTEST_ccfsg, CUTEST_cdimsj);
  runNope(nope);
  ppDIMEN(nope, &nvar, &ncon);

  dci.cuter ();

  Real x[nvar], bl[nvar], bu[nvar];

  if (ncon == 0) {
    dci.set_ufn(ufn);
    dci.set_uofg(uofg);
    dci.set_uprod(uhprod);
    dci.set_unames(CUTEST_unames);

    runUncSetup(nope, &nvar, x, bl, bu);
    dci.unc_setup(nvar, x, bl, bu);
  } else {
    Real y[ncon], cl[ncon], cu[ncon];
    Bool equatn[ncon], linear[ncon];

    dci.set_cfn(cfn);
    dci.set_cofg(cofg);
    dci.set_cprod(chprod);
    dci.set_ccfsg(ccfsg);
    dci.set_cnames(CUTEST_cnames);

//    dci.set_ccifg (CUTEST_ccifg);

    runConSetup(nope, &nvar, x, bl, bu, &ncon, y, cl, cu, equatn, linear, &amax);
    dci.set_linear (ncon, linear);
    dci.set_amax (amax);
    dci.con_setup (nvar, x, bl, bu, ncon, y, cl, cu, equatn);
  }

  try {
    dci.start ();
    dci.solve ();
    dci.show ();
    dci.printLatex ();
  } catch (const char * ex) {
    std::cout << ex << std::endl;
    return 1;
  } catch (...) {
    std::cout << "Unhandled exception caught" << std::endl;
    return 1;
  }

  Real calls[7], time[2];
  if (ncon == 0)
    CUTEST_ureport(&status, calls, time);
  else
    CUTEST_creport(&status, calls, time);
  std::cout << "Setup time: " << time[0] << std::endl
            << "Execution time: " << time[1] << std::endl;
  FORTRAN_close ((&funit), (&ierr));

  destroyNope(nope);

  return 0;

}
