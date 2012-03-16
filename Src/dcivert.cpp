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
    Real oldnormc, dnnorm;
    Int trflag = 0;
    Int ibfgs = 0, iout = 0;
    Bool dnavail = dciFalse;
    Bool scaleJ = dciTrue;;
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
    xc->scale (*x, 1);
    if (Ineq)
      sc->scale (*s, 1);
    VertFlag = 0;
    if (iter == 1) {
      Aavail = dciTrue;
      gavail = dciTrue;
    } else {
      Aavail = dciFalse;
      gavail = dciFalse;
    }

    if ( (normc <= rho) && (!Aavail) ) {
      if (ncon > 0) {
        if (!lincon)
          call_ccfsg_xc (); //CuterJacob
        Aavail = dciTrue;

        call_ofg_xc (); //Just g

        gavail = dciTrue;

        if (!lincon) {
          analyze_J ();
          cholesky_J ();
        }

        LimLbd = dciTrue;
        update_lambda ();
        updyineq ();
      } else {
        call_ofg_xc ();
        gp->scale (*g, 1);
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
        GDBSTOP ();
#endif
        nRest++;

        trflag = dcitrust (oldnormc);

        checkInfactibility ();

        gavail = dciFalse;

        if (normc > phi2*oldnormc) {
          fail = fail + 1;
        } else
          fail = 0;

        if (fail >= nfailv) {
#ifdef VERBOSE
          std::cout << "Going to Safe Guard " << std::endl
                    << std::endl;
          if ( (nvar < 10) && (ncon < 10) ) {
            std::cout <<  "A = " << std::endl;
            full(*J).print_more ();
          }
          GDBSTOP ();
#endif
          VertFlag = vertSafeguard ();
        } else if ( ( (normc > thetaR*oldnormc) && (oldAcnt > 0) ) || (oldAcnt > 5) || (iout == 5) ) {
          // dcivert failed. Recompute A

          if (!lincon) {
            call_ccfsg_xc (dciTrue, scaleJ); //CuterJacob
          }
          Aavail = dciTrue;
          oldAcnt = 0;

          if (!lincon) {
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
        if (!lincon) {
          call_ccfsg_xc (); //CuterJacob
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
