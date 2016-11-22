#include <iostream>
#include "dci.h"
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <cmath>
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

  if (nvar < nope->nvar)
    std::cout << "NOPE preprocessing removed " << nope->nvar - nvar << " variables" << std::endl;
  if (ncon < nope->ncon)
    std::cout << "NOPE preprocessing removed " << nope->ncon - ncon << " contraints" << std::endl;

  Real x[nvar], bl[nvar], bu[nvar];

  if (nvar > 0) {
    DCI::Interface dci;
    dci.cuter ();

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
    try {
      std::cout << "Starting" << std::endl;
      dci.start ();
      std::cout << "Solving" << std::endl;
      dci.solve ();
      std::cout << "Showing" << std::endl;
      dci.show ();

      if (dci.get_display_level() > 2) {
        for (int i = 0; i < nope->nfix; i++) {
          int j = nope->fixed_index[i];
          printf("fixed x%d = %f\n", j+1, nope->x[j]);
        }
      }

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

  } else {
    // No variables after preprocessing
    Int status = 0;
    char pname[10], vnames[10*nope->nvar];
    Real f = 1e20;
    if (nope->ncon > 0) {
      char cnames[10*nope->ncon];

      CUTEST_cnames(&status, &nope->nvar, &nope->ncon, pname, vnames, cnames);
      CUTEST_cfn(&status, &nope->nvar, &nope->ncon, nope->x, &f, nope->c);

      for (int i = 0; i < ncon; i++)
        if (nope->equatn[i])
          assert(nope->c[i] == 0.0);
        else
          assert((nope->c[i] >= nope->cl[i]) && (nope->c[i] <= nope->cu[i]));
    } else {
      CUTEST_unames(&status, &nope->nvar, pname, vnames);
      CUTEST_ufn(&status, &nope->nvar, nope->x, &f);
    }

    for (int i = 0; i < nope->nvar; i++) {
      assert(nope->x[i] >= nope->original_bl[i]);
      assert(nope->x[i] <= nope->original_bu[i]);
      if (fabs(nope->x[i]) > 0.0)
        printf("x[%i] = %f\n", i+1, nope->x[i]);
    }

    std::cout << "EXIT: After the preprocessing, there no variables" << std::endl
              << "f(x) = " << f << std::endl
              << "Elapsed Time = 0.0" << std::endl
              << "WARNING: Time of the preprocessing not included" << std::endl;

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
