#include "interface.h"
#include <cmath>
//#include <cassert>

namespace DCI {
  /* ************************************************
   * ******      VertStep : Vertical Step      ******
   * ************************************************
   *
   * This function finds the vertical step. It is supposed
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

  Int Interface::vertstep () {
//    Vector dn (*env);
    Int fail = 0;
    Int oldAcnt = 1;
    Real oldnormc;
//    Real dnnorm;
//    Int ibfgs = 0;
    Int iout = 0;
//    Bool dnavail = dciFalse;
//    Bool scaleJ = dciTrue;;
    Vector gtmp(*env);

    ngp = engp;
    rho = Min(phi1*rhomax*ngp, 0.75*rhomax);
    if ( (rho < csic) && (normc > 100*csic) )
      rho = 0.1*normc;
    else if (ngp <= 5*csig)
      rho = csic;
    rho = Max (rho, csic);
    update_mu ();
    nRest = 0;
    nbfgs = 0;
    oldnormc = normc;
    *xc = *x;
    if (Ineq)
      *sc = *s;
    VertFlag = 0;
    if (iter == 1) {
      Aavail = dciTrue;
      gavail = dciTrue;
    } else {
      Aavail = dciFalse;
      gavail = dciFalse;
    }

    Real one[2] = {1,0}, zero[2] = {0,0};
    Vector gradLeastSquare(*env);
    if (ncon) {
      gradLeastSquare.sdmult(*J, 1, one, zero, *c);
      DeltaV = gradLeastSquare.norm();
    }


    if ( (normc <= rho) && (!Aavail) ) {
      if (ncon > 0) {
        if (!Linear)
          call_ccfsg_xc (dciTrue); //CuterJacob
        Aavail = dciTrue;

        call_ofg_xc (); //Just g

        gavail = dciTrue;

        if (!Linear && !UseMUMPS) {
          analyze_J ();
          cholesky_J ();
        }

        LimLbd = dciTrue;
        update_lambda ();
        updyineq ();
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

      update_mu ();
      oldAcnt = 0;
    }

    CurrentTime = getTime() - StartTime;
    while ( (normc > rho) && (nRest <= maxrest) && (VertFlag == 0) && (CurrentTime < MaxTime) ) {

      while ( (normc > rho) && (nRest <= maxrest) && (VertFlag == 0) && (CurrentTime < MaxTime) ) {
        
#ifdef VERBOSE
        if (VerboseLevel > 1) {
          std::cout << "Going to dcitrust: nRest " << nRest << std::endl
                    << std::endl
                    << "|c| = " << normc << std::endl
                    << "rho = " << rho << std::endl
                    << std::endl;
          if ( (nvar < 10) && (ncon < 10) ) {
            std::cout <<  "A = " << std::endl;
            full(*J).print_more ();
            std::cout << "xc = " << std::endl;
            xc->print_more ();
            if (Ineq) {
              std::cout << "sc = " << std::endl;
              sc->print_more ();
            }
          }
        }
        GDBSTOP ();
#endif
        nRest++;

        if (UsePorcelli) {
          Porcelli();
        } else if ((!Bounded) || (!UseVertInteriorPoint) )
            dcitrust (oldnormc);
        else
          InteriorPointRestoration ();

#ifndef NDEBUG
        checkInfactibility();
#endif

        gavail = dciFalse;

        if (normc > phi2*oldnormc) {
          fail = fail + 1;
        } else
          fail = 0;

        if (fail >= nfailv) {
#ifdef VERBOSE
          if (VerboseLevel > 1) {
            std::cout << "Going to Safe Guard " << std::endl
                      << std::endl;
            if ( (nvar < 10) && (ncon < 10) ) {
              std::cout <<  "A = " << std::endl;
              full(*J).print_more ();
            }
            GDBSTOP ();
          }
//          std::cout << "Entering fail at dcivert" << std::endl;
#endif
        
          call_ccfsg_xc (dciTrue, ScaleVertical);
          if (calc_feasibilityOpt())
            VertFlag = 2;
//          if (!UseMUMPS)
//            this->cholesky_J();
//          if (UseVertSafeguard)
//            VertFlag = vertSafeguard ();
//          VertFlag = 0;
        } else if ( ( (normc > thetaR*oldnormc) && (oldAcnt > 0) ) || (oldAcnt > 5) || (iout == 5) ) {
          // dcivert failed. Recompute A

          if (!Linear) {
            call_ccfsg_xc (dciTrue, ScaleVertical); //CuterJacob
          }
          Aavail = dciTrue;
          oldAcnt = 0;

          if (!Linear && !UseMUMPS) {
            this->cholesky_J ();
          }

        } else {
          oldAcnt++;
          Aavail = dciFalse;
        }

        oldnormc = normc;
        DeltaV = Max (DeltaV, DeltaMin);

        CurrentTime = getTime() - StartTime;

      } //Fim do While

      if (!Aavail) {
        if (!Linear) {
          call_ccfsg_xc (dciTrue); //CuterJacob
          this->cholesky_J ();
        }
        Aavail = dciTrue;
      }

      if (!gavail) {
        call_ofg_xc (); //Just g

        LimLbd = dciTrue;
        update_lambda ();
        updyineq ();

        normgp = gp->norm ();
        normg = g->norm ();
        ngp = normgp/(normg + 1);
      }
      
      if (rho > 2*rhomax*ngp)
        rho = Min (phi1*rhomax*ngp, Max (1e-4*rhomax*ngp, 0.75*rhomax) );
      else
        rho = Max (rho, Min (phi1*rhomax*ngp, Max (1e-4*rhomax*ngp, 0.75*rhomax) ) );
      rho = Max (rho, csic);

      update_mu ();
    } //Fim do While

    if (ncon > 0)
      normy = y->norm ();

//    tmp5 = xc - x

    return 0;
  }
}
