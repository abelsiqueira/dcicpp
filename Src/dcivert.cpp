#include "interface.h"
#include <cmath>
#include <fstream>
#include <algorithm>
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

  Int Interface::verticalStep () {
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
    updateMu ();
    nRest = 0;
    nbfgs = 0;
    oldnormc = normc;
    *xc = *x;
    if (has_ineq)
      *sc = *s;
    VertFlag = 0;
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

        if (!is_linear) {
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
      oldAcnt = 0;
    }

    current_time = getTime() - start_time;
    while ( (normc > rho) && (nRest <= maxrest) && (VertFlag == 0) && (current_time < max_time) ) {

      while ( (normc > rho) && (nRest <= maxrest) && (VertFlag == 0) && (current_time < max_time) ) {
        
#ifdef VERBOSE
        if (verbosity_level > 1) {
          std::cout << "Going to innerVerticalStep: nRest " << nRest << std::endl
                    << std::endl
                    << "|c| = " << normc << std::endl
                    << "rho = " << rho << std::endl
                    << std::endl;
          if ( (nvar < 10) && (ncon < 10) ) {
            std::cout <<  "A = " << std::endl;
            full(*J).print_more ();
            std::cout << "xc = " << std::endl;
            xc->print_more ();
            if (has_ineq) {
              std::cout << "sc = " << std::endl;
              sc->print_more ();
            }
          }
        }
        GDBSTOP ();
#endif
        nRest++;

        call_ccfsg_xc(dciTrue, dciFalse);
        cholesky();
        infeasible_gradient = 1.0;

        innerVerticalStep(infeasible_gradient);
#ifdef ITER_MATLAB
        iter_file << "X(:,size(X,2)+1) = [" << xcx[0] << ";" << xcx[1] << "];" << std::endl;
#endif

#ifdef VERBOSE
        if (verbosity_level > 1) {
          std::cout << "After innerVerticalStep" << std::endl;
          std::cout << "|c| = " << normc << std::endl;
          std::cout << "rho = " << rho << std::endl;
          if ( (nvar < 10) && (ncon < 10) ) {
            std::cout << "xc = " << std::endl;
            xc->print_more();
            if (has_ineq) {
              std::cout << "sc = " << std::endl;
              sc->print_more();
            }
          }
        }
#endif

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
          if (verbosity_level > 1) {
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
        
          call_ccfsg_xc (dciTrue, dciFalse);
          if (normc > 0 && infeasible_gradient/normc < 1e-6)
            VertFlag = 2;

          Vector ssoc(*env, nvar + nconI);
          Real asoc;
          call_ccfsg_xc(dciTrue, dciTrue);
          cholesky();
          StepFlag = naStep (*c, ssoc);
          scale_xc (ssoc);

          // Arrumar tamanho do ssoc a partir do x
          Real alphassoc = 1;
          pReal ssocx = ssoc.get_doublex();
          for (Int i = 0; i < nvar+nconI; i++) {
            Real xi = xcx[i], bli = l_bndx[i], bui = u_bndx[i], di = ssocx[i];
            if (di == 0)
              continue;
            if (di < 0) {
              Real val = (bli - xi)*(1 - epsmu)/di;
                alphassoc = Min (alphassoc, val);
            } else {
              Real val = (bui - xi)*(1 - epsmu)/di;
                alphassoc = Min (alphassoc, val);
            }
          }
          if (alphassoc < 1)
            ssoc.scale (alphassoc);

          asoc = ssoc.norm (0);
          if (asoc > DeltaV)
            ssoc.scale (DeltaV/asoc);
          for (Int i = 0; i < nvar+nconI; i++)
            xcx[i] += ssocx[i];
          copy_scx();
          call_fn ();
          normc = c->norm ();
          fail = 0;

          if (vertical_fail_reboot && VertFlag == 0) {
            // Has failed but is not infeasible
            Real constr[ncon], funval;
            (*cfn) (&nvar, &ncon, xcx, &funval, &mmax, constr);
            Int numI = nvar;
            for (Int i = 0; i < ncon; i++) {
              if (equatn[i] == dciFalse) {
                if (constr[i] > clx[i] && constr[i] < cux[i])
                  xcx[numI] = constr[i];
                numI++;
              }
            }
            copy_scx();
            normc = c->norm();
          }
        } else if ( ( (normc > thetaR*oldnormc) && (oldAcnt > 0) ) || (oldAcnt > 5) || (iout == 5) ) {
          // dcivert failed. Recompute A

          if (!is_linear) {
            call_ccfsg_xc (dciTrue, dciFalse); //CuterJacob
          }
          Aavail = dciTrue;
          oldAcnt = 0;

          if (!is_linear) {
            this->cholesky ();
          }

        } else {
          oldAcnt++;
          Aavail = dciFalse;
        }

        oldnormc = normc;
        DeltaV = Max (DeltaV, DeltaMin);

        current_time = getTime() - start_time;

      } //Fim do While

      if (!Aavail) {
        if (!is_linear) {
          call_ccfsg_xc (dciTrue); //CuterJacob
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
