#include <iostream>
#include "dci.h"
extern "C" {
#include "cutest.h"
}

using namespace DCI;

int MAINENTRY () {
  DCI::Interface dci;
  char fname[10] = "OUTSDIF.d";
  Int nvar = 0, ncon = 0, amax = 0;
  Int funit = 42, ierr = 0, fout = 6, io_buffer = 11, status;

  FORTRAN_open (&funit, fname, &ierr);
  CUTEST_cdimen (&status, &funit, &nvar, &ncon);

  dci.cuter ();

  Real x[nvar], bl[nvar], bu[nvar];

  if (ncon == 0) {
    CUTEST_usetup (&status, &funit, &fout, &io_buffer, &nvar, x, bl, bu);

    dci.set_uofg (CUTEST_uofg);
    dci.set_uprod (CUTEST_uhprod);
    dci.set_ufn (CUTEST_ufn);
    dci.set_unames (CUTEST_unames);

    dci.unc_setup (nvar, x, bl, bu);
  } else {
    Real y[ncon], cl[ncon], cu[ncon];
    Bool equatn[ncon], linear[ncon];
    Int efirst = 0, lfirst = 0, nvfirst = 0;

    CUTEST_csetup (&status, &funit, &fout, &io_buffer, &nvar, &ncon, x, bl, bu,
        y, cl, cu, equatn, linear, &efirst, &lfirst, &nvfirst);
    CUTEST_cdimsj (&status, &amax);
    amax *= 2;
  
    dci.set_cofg (CUTEST_cofg);
    dci.set_cprod (CUTEST_chprod);
    dci.set_cfn (CUTEST_cfn);
    dci.set_ccfsg (CUTEST_ccfsg);
    dci.set_ccifg (CUTEST_ccifg);
    dci.set_cnames (CUTEST_cnames);

    dci.set_linear (ncon, linear);
    dci.set_amax (amax);
    dci.con_setup (nvar, x, bl, bu, ncon, y, cl, cu, equatn);
  }

  dci.start ();
  try {
    dci.solve ();
    dci.show ();
    dci.printLatex ();
  } catch (const char * ex) {
    std::cout << ex << std::endl;
  } catch (...) {
    std::cout << "Unhandled exception caught" << std::endl;
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
