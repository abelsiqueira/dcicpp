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

        call_ccfsg_xc(dciTrue, dciFalse);
        cholesky_J();
        infeasible_gradient = 1.0;

        if (UsePorcelli) {
          Porcelli(infeasible_gradient);
        } else if ((!Bounded) || (!UseVertInteriorPoint) )
            dcitrust (oldnormc);
        else
          InteriorPointRestoration ();
#ifdef ITER_MATLAB
    iter_file << "X(:,size(X,2)+1) = [" << xcx[0] << ";" << xcx[1] << "];" << std::endl;
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
        
          call_ccfsg_xc (dciTrue, dciFalse);
          if (normc > 0 && infeasible_gradient/normc < 1e-6)
            VertFlag = 2;

          Vector ssoc(*env, nvar + nconI);
          Real asoc;
          call_ccfsg_xc(dciTrue, dciTrue);
          cholesky_J();
          StepFlag = NAstep (*c, ssoc);
          scale_xc (ssoc);

          // Arrumar tamanho do ssoc a partir do x
          Real alphassoc = 1;
          pReal ssocx = ssoc.get_doublex();
          for (Int i = 0; i < nvar; i++) {
            Real xi = xcx[i], bli = blx[i], bui = bux[i], di = ssocx[i];
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
          for (Int i = 0; i < nconI; i++) {
            Real si = scx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], di = ssocx[nvar + i];
            if (di == 0)
              continue;
            if (di < 0) {
              Real val = (cli - si)*(1 - epsmu)/di;
                alphassoc = Min (alphassoc, val);
            } else {
              Real val = (cui - si)*(1 - epsmu)/di;
                alphassoc = Min (alphassoc, val);
            }
          }
          if (alphassoc < 1)
            ssoc.scale (alphassoc);

          asoc = ssoc.norm (0);
          if (asoc > DeltaV)
            ssoc.scale (DeltaV/asoc);
          for (Int i = 0; i < nvar; i++)
            xcx[i] += ssocx[i];
          for (Int i = 0; i < nconI; i++) {
            scx[i] += ssocx[nvar + i];
          }
          call_fn ();
          normc = c->norm ();
          fail = 0;


          if (RebootOnVertFail && VertFlag == 0) {
            // Has failed but is not infeasible
            Real constr[ncon], funval;
            (*cfn) (&nvar, &ncon, xcx, &funval, &mmax, constr);
            Int numI = 0;
            for (Int i = 0; i < ncon; i++) {
              if (equatn[i] == dciFalse) {
                if (constr[i] > clx[i] && constr[i] < cux[i])
                  scx[numI] = constr[i];
                numI++;
              }
            }
            normc = c->norm();
          }



//          if (!UseMUMPS)
//            this->cholesky_J();
//          if (UseVertSafeguard)
//            VertFlag = vertSafeguard ();
//          VertFlag = 0;
        } else if ( ( (normc > thetaR*oldnormc) && (oldAcnt > 0) ) || (oldAcnt > 5) || (iout == 5) ) {
          // dcivert failed. Recompute A

          if (!Linear) {
            call_ccfsg_xc (dciTrue, dciFalse); //CuterJacob
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
