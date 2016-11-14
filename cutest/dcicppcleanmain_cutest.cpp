#include <iostream>
#include "dci.h"
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <cmath>
extern "C" {
#include "cutest.h"
}

using namespace DCI;

int MAINENTRY () {
  char fname[10] = "OUTSDIF.d";
  Int nvar = 0, ncon = 0, amax = 0;
  Int funit = 42, ierr = 0, io_buffer = 11, status;

  FORTRAN_open (&funit, fname, &ierr);
  CUTEST_cdimen(&status, &funit, &nvar, &ncon);

  Real x[nvar], bl[nvar], bu[nvar];

  DCI::Interface dci;
  dci.cuter ();

  if (ncon == 0) {
    dci.set_ufn(CUTEST_ufn);
    dci.set_uofg(CUTEST_uofg);
    dci.set_uprod(CUTEST_uhprod);
    dci.set_unames(CUTEST_unames);

    CUTEST_usetup(&status, &funit, &ierr, &io_buffer, &nvar, x, bl, bu);
    dci.unc_setup(nvar, x, bl, bu);
  } else {
    Real y[ncon], cl[ncon], cu[ncon];
    Bool equatn[ncon], linear[ncon];
    Int efirst = 0, lfirst = 1, nvfirst = 1;

    dci.set_cfn(CUTEST_cfn);
    dci.set_cofg(CUTEST_cofg);
    dci.set_cprod(CUTEST_chprod);
    dci.set_ccfsg(CUTEST_ccfsg);
    dci.set_cnames(CUTEST_cnames);

    //    dci.set_ccifg (CUTEST_ccifg);

    CUTEST_csetup(&status, &funit, &ierr, &io_buffer, &nvar, &ncon, x, bl, bu,
        y, cl, cu, equatn, linear, &efirst, &lfirst, &nvfirst);
    CUTEST_cdimsj(&status, &amax);
    dci.set_linear (ncon, linear);
    dci.set_amax (amax);
    dci.con_setup (nvar, x, bl, bu, ncon, y, cl, cu, equatn);

  }
  try {
    std::cout << "Starting" << std::endl;
    dci.start ();
    std::cout << "Solving" << std::endl;
    dci.solve ();
    std::cout << "Showing" << std::endl;
    dci.show ();
    std::cout << "Done" << std::endl;
    // dci.printLatex ();
  } catch (const char * ex) {
    std::cout << "Exception: " << ex << std::endl;
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

  return 0;

}
