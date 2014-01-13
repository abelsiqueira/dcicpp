#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  /* *************************************************
   * ******      TangentStep: Tangent Step      ******
   * *************************************************
   *
   * This function finds the Tangent step. It is
   * supposed to find a miniminum to an approximation
   * of f around the current point, with the conditions
   * of being in a trust region, doesn't leave the 
   * cylinder 2*rho, and is in the null space of the
   * Jacobian matrix.
   */

  void Interface::tangentStep (Real & normd) {
    Vector d (*env, nvar + nconI), ssoc (*env, nvar + nconI);
    Real qd = 0;
    Real newnormc, asoc, alphaT;
    Int SteihFlag;
    Real gtd;

    GotH = dciFalse;
    first = dciTrue;
    nRej = 0;
    DLH = -1;
    qd = 2*DLH/eta1;
    newnormc = 2*rho + 1;
    nSoc = 0;
//    Real smlAlpha = 1e-3;
    Bool fail = dciFalse;

#ifdef VERBOSE
    if (verbosity_level > 1) {
      if (ncon > 0 && nvar + ncon < 10) {
        std::cout << "A = " << std::endl;
        full(*J).print_more ();
      }
    }
#endif

    current_time = getTime() - start_time;
    while ( (!fail) && 
            ( (newnormc > zeta1*rho) || 
              (DLH > eta1*qd) ) && 
            (current_time < max_time) ) {

      SteihFlag = innerTangentStep (d, qd, gtd);


      pReal dx = d.get_doublex();
      scale_xc (d);
      for (Int i = 0; i < nvar+nconI; i++) {
        if (l_bndx[i] - u_bndx[i] > -dciEps) {
          dx[i] = 0;
          continue;
        }
        xx[i] = xcx[i] + dx[i];
        if (xx[i] == u_bndx[i])
          xx[i] = u_bndx[i] - dciEps;
        else if (xx[i] == l_bndx[i])
          xx[i] = l_bndx[i] + dciEps;
      }

#ifdef VERBOSE
      if (verbosity_level > 1) {
        std::cout << "SteihFlag = " << SteihFlag 
                  << ", nSteih = " << nSteih << std::endl;
      }
      if ( (nvar + ncon <= 5) && (verbosity_level > 1) ) {
        std::cout << "DeltaH*Diag^-1 = " << std::endl;
        for (Int i = 0; i < nvar+ncon; i++) {
          std::cout << DeltaH/scaling_matrix[i] << std::endl;
        }
        std::cout << "x+ = " << std::endl;
        x->print_more();
      }
#endif

      Vector tmp_c(*env);
      if (ncon > 0)
        tmp_c = *c;
      call_fn ();
#ifndef NDEBUG
      try {
        checkInfactibility ();
      } catch (const char * ex) {
        DeltaH *= 0.1;
        continue;
      }
#endif

      if (ncon > 0)
        newnormc = c->norm ();
      else
        newnormc = 0;

      if ( (SteihFlag == 1) || (SteihFlag == 2) )
        normd = DeltaH;
      else
        normd = d.norm (0);

      if ( first && use_soc && 
         ( (newnormc > Min (zeta1*rho, zeta1*normc + zeta2*rho) ) || 
           ( (normc <= csic) && 
             (newnormc > Max (csic, zeta3*normc) ) ) ) ) {

        pReal ptmp_c = tmp_c.get_doublex();
        for (Int i = 0; i < ncon; i++)
          ptmp_c[i] = cx[i] - ptmp_c[i];
        StepFlag = naStep (tmp_c, ssoc);
        scale_xc (ssoc);

        // Arrumar tamanho do ssoc a partir do x
        Real alphassoc = 1;
        pReal ssocx = ssoc.get_doublex();
        for (Int i = 0; i < nvar+nconI; i++) {
          Real xi = xx[i], bli = l_bndx[i], bui = u_bndx[i], di = ssocx[i];
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

        ssoc.saxpy (d, 1);
        asoc = ssoc.norm (0);
        *x = *xc;
        if (asoc > DeltaH)
          x->saxpy(ssoc, DeltaH/asoc);
        else
          x->saxpy(ssoc, 1.0);

#ifdef VERBOSE
        if (verbosity_level > 1) {
          std::cout << "SOC = " << std::endl;
        }
        if ( (nvar + ncon <= 5) && (verbosity_level > 1) ) {
          std::cout << "x+ = " << std::endl;
          x->print_more();
        }
#endif
        call_fn ();
        newnormc = c->norm ();
        nSoc++;
      }

//      f is already calculated in fn
      if (ncon > 0)
        Ln = *f + y->dot (*c);
      else
        Ln = *f;
      DLH = Ln - Lc;

      if ( (newnormc > zeta1*rho) || (DLH > eta1*qd) ) {
        if (DLH > 0) {
          alphaT = (1 - eta3) * gtd / ( (1 - eta3) * (Lc + gtd) - Ln + eta2*qd);
//          DeltaH = Min (alphaR*normd, Max (alphaT, alphaS) * DeltaH);
          DeltaH = Min( Min (alphaR*normd, Max (alphaT, alphaS) * DeltaH), 0.9*DeltaH);
        } else {
          if (DeltaH == normd * alphaR)
            fail = dciTrue;
          DeltaH = normd * alphaR;
        }
        nRej++;
#ifdef VERBOSE
        if (verbosity_level > 1) {
          std::cout << "Passo rejeitado" << std::endl;
          if (newnormc > zeta1*rho)
            std::cout << "-> Fora do cilindro" << std::endl;
          if (DLH > eta1*qd)
            std::cout << "-> Decrescimo insuficiente" << std::endl;
        }
#endif
      } else if (DLH <= eta2*qd)
        DeltaH = Max (DeltaH, normd*alphaI);

//      if (normc > 0 && d.norm() < 1e-12 &&  infeasible_gradient/normc < 1e-6)
//        NormalFlag = 2;

      first = dciFalse;
      if (DeltaH < 1e-10)
        break;

      current_time = getTime() - start_time;
    }
    updateScaling_x();

    normc = newnormc;

  }
}
