#include <iostream>
#include "dci.h"
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
extern "C" {
#include "nope.h"
#include "cutest.h"
}

using namespace DCI;

Nope *nope;
#include "nope_interface.c"

int MAINENTRY () {
  char fname[10] = "OUTSDIF.d";
  Int nvar = 0, ncon = 0, amax = 0;
  Int funit = 42, ierr = 0, status;

  nope = initializeNope();

  FORTRAN_open (&funit, fname, &ierr);
  setFuncs(nope, CUTEST_cdimen, CUTEST_usetup, CUTEST_unames, CUTEST_ufn,
      CUTEST_uofg, CUTEST_uhprod, CUTEST_csetup, CUTEST_cnames, CUTEST_cfn,
      CUTEST_cofg, CUTEST_chprod, CUTEST_ccfsg, CUTEST_cdimsj);
  runNope(nope);
  ppDIMEN(nope, &nvar, &ncon);

  if (nvar > 0) {
    DCI::Interface dci;
    dci.cuter ();

    Real x[nvar], bl[nvar], bu[nvar];

    if (ncon == 0) {
      dci.set_ufn(ufn);
      dci.set_uofg(uofg);
      dci.set_uprod(uhprod);
      dci.set_unames(unames);

      runUncSetup(nope, &nvar, x, bl, bu);
      dci.unc_setup(nvar, x, bl, bu);
    } else {
      Real y[ncon], cl[ncon], cu[ncon];
      Bool equatn[ncon], linear[ncon];

      dci.set_cfn(cfn);
      dci.set_cofg(cofg);
      dci.set_cprod(chprod);
      dci.set_ccfsg(ccfsg);
      dci.set_cnames(cnames);

      //    dci.set_ccifg (CUTEST_ccifg);

      runConSetup(nope, &nvar, x, bl, bu, &ncon, y, cl, cu, equatn, linear, &amax);
      dci.set_linear (ncon, linear);
      dci.set_amax (amax);
      dci.con_setup (nvar, x, bl, bu, ncon, y, cl, cu, equatn);

    }
#ifndef FAIL_ON_EXCEPTION
    try {
#endif
      dci.start ();
      dci.solve ();
      dci.show ();
      dci.printLatex ();
#ifndef FAIL_ON_EXCEPTION
    } catch (const char * ex) {
      std::cout << ex << std::endl;
      return 1;
    } catch (...) {
      std::cout << "Unhandled exception caught" << std::endl;
      return 1;
    }
#endif
    Real calls[7], time[2];
    if (ncon == 0)
      CUTEST_ureport(&status, calls, time);
    else
      CUTEST_creport(&status, calls, time);
    std::cout << "Setup time: " << time[0] << std::endl
      << "Execution time: " << time[1] << std::endl;
  } else {
    Int status = 0;
    char pname[10], vnames[10*nope->nvar];
    Real f = 1e20;
    if (nope->ncon > 0) {
      char cnames[10*nope->ncon];
      CUTEST_cnames(&status, &nope->nvar, &nope->ncon, pname, vnames, cnames);
      CUTEST_cfn(&status, &nope->nvar, &nope->ncon, nope->x, &f, nope->c);
    } else {
      CUTEST_unames(&status, &nope->nvar, pname, vnames);
      CUTEST_ufn(&status, &nope->nvar, nope->x, &f);
    }
    Real calls[7], time[2];
    if (ncon == 0)
      CUTEST_ureport(&status, calls, time);
    else
      CUTEST_creport(&status, calls, time);
    std::cout << "Setup time: " << time[0] << std::endl
      << "Execution time: " << time[1] << std::endl;

    pname[8] = 0;
    std::string problemName(pname);
    size_t endpos = problemName.find_last_not_of("\t");
    std::ofstream file;
    file.open((problemName.substr(0,endpos+1)+".tab").c_str(),
        std::ios_base::out | std::ios_base::trunc);
    file << problemName << " convergence "
      << std::scientific << std::setprecision(6)
        << (time[0]+time[1] > 0 ? time[0]+time[1] : 1e-6) << " "
      << std::scientific << std::setprecision(6) << f
      << " 0.000000e+00 0.000000e+00 cholok\n";
  }

  FORTRAN_close ((&funit), (&ierr));

  destroyNope(nope);

  return 0;

}
