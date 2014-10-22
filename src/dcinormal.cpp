#include "interface.h"
#include <cmath>
#include <fstream>
#include <algorithm>
//#include <cassert>

namespace DCI {
  /* ************************************************
   * ******      NormalStep : Normal Step      ******
   * ************************************************
   *
   * This function finds the normal step. It is supposed
   * to find a point x such that
   *
   * ||h(x)|| <= rho
   *
   * To do this, we will sucessively apply Powell's dogleg
   * method to the problem
   *
   * min 0.5||h(x)+As||^2 - 0.5||h(x)||
   * s.t. ||d|| <= Delta
   *
   * until ||h(x)|| <= rho
   *
   * Todo:
   * - The method itself
   * - Return gp too
   */

  Int Interface::normalStep () {
    Vector gtmp(*env);

    ngp = engp;
    rho = Min(phi1*rhomax*ngp, 0.75*rhomax);
    if ( (rho < csic) && (normc > 100*csic) )
      rho = 0.1*normc;
    else if (ngp <= 5*csig)
      rho = csic;
    rho = Max (rho, csic);
    updateMu ();
    nRest = 0;
    *xc = *x;
    NormalFlag = 0;
    if (iter == 1) {
      Aavail = dciTrue;
      gavail = dciTrue;
    } else {
      Aavail = dciFalse;
      gavail = dciFalse;
    }

#ifdef ITER_MATLAB
    std::string iter_filename("iter_matlab");
    iter_filename += problemName;
    iter_filename += ".m";
    iter_filename.erase(remove_if(iter_filename.begin(),
          iter_filename.end(), isspace), iter_filename.end());
    std::ofstream iter_file(iter_filename.c_str());
    iter_file << "X = [" << xcx[0] << ";" << xcx[1] << "];" << std::endl;
#endif

    Real one[2] = {1,0}, zero[2] = {0,0};
    Vector gradLeastSquare(*env);
    if (ncon) {
      gradLeastSquare.sdmult(*J, 1, one, zero, *c);
//      DeltaV = gradLeastSquare.norm();
    }

    if ( (normc <= rho) && (!Aavail) ) {
      if (ncon > 0) {
        if (!is_linear)
          call_ccfsg_xc (dciTrue); //CuterJacob
        Aavail = dciTrue;

        call_ofg_xc (); //Just g

        gavail = dciTrue;

        if (!is_linear && !use_lsmr) {
          analyzeJacobian ();
          cholesky ();
        }

        LimLbd = dciTrue;
        updateMultipliers ();
        update_yineq ();
      } else {
        call_ofg_xc ();
        *gp = *g;
      }

      normgp = gp->norm ();
      ngp = normgp/(g->norm() + 1);
      rho = Min (phi1*rhomax*ngp, 0.75*rhomax);
      if ( (rho < csic) && (normc > 100*csic) )
        rho = 0.1*normc;
      else if (ngp <= 5*csig)
        rho = csic;
      rho = Max (rho, csic);

      updateMu ();
//      oldAcnt = 0;
    }

    current_time = getTime() - start_time;
    while ( (normc > rho) && (nRest <= maxrest) && (NormalFlag == 0) && (current_time < max_time) ) {

      innerNormalPhase();

      if (!Aavail) {
        if (!is_linear) {
          call_ccfsg_xc (dciTrue); //CuterJacob
          if (!use_lsmr)
            this->cholesky ();
        }
        Aavail = dciTrue;
      }

      if (!gavail) {
        call_ofg_xc (); //Just g

        LimLbd = dciTrue;
        updateMultipliers ();
        update_yineq ();

        normgp = gp->norm ();
        normg = g->norm ();
        ngp = normgp/(normg + 1);
      }

      if (rho > 2*rhomax*ngp)
        rho = Min (phi1*rhomax*ngp, Max (1e-4*rhomax*ngp, 0.75*rhomax) );
      else
        rho = Max (rho, Min (phi1*rhomax*ngp, Max (1e-4*rhomax*ngp, 0.75*rhomax) ) );
      rho = Max (rho, csic);
#ifdef VERBOSE
      if (verbosity_level > 1) {
        std::cout << "After last update"
            << std::endl
            << "|c| = " << normc << std::endl
            << "rho = " << rho << std::endl;
      }
#endif

      updateMu ();
    } //Fim do While

    if (ncon > 0)
      normy = y->norm ();

//    tmp5 = xc - x

    return 0;
  }
}
