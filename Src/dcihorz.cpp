#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  /* *************************************************
   * ******      HorzStep: Horizontal Step      ******
   * *************************************************
   *
   * This function finds the horizontal step. It is
   * supposed to find a miniminum to an approximation
   * of f around the current point, with the conditions
   * of being in a trust region, doesn't leave the 
   * cylinder 2*rho, and is in the null space of the
   * Jacobian matrix.
   */

  void Interface::horzstep (Real & normd) {
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
    Real smlAlpha = 1e-3;
    Bool fail = dciFalse;

#ifdef VERBOSE
    if (VerboseLevel > 1) {
      if (nvar + ncon < 10) {
        std::cout << "A = " << std::endl;
        full(*J).print_more ();
      }
    }
#endif

    CurrentTime = getTime() - StartTime;
    while ( (!fail) && 
            ( (newnormc > zeta1*rho) || 
              (DLH > eta1*qd) ) && 
            (CurrentTime < MaxTime) ) {

      SteihFlag = dcisteih (d, qd, gtd);

      *x = *xc;
      if (Ineq)
        *s = *sc;
      pReal dx = d.get_doublex();
      scale_xc (d);
      for (Int i = 0; i < nvar; i++)
        xx[i] += dx[i];
      for (Int i = 0; i < nconI; i++)
        sx[i] += dx[nvar + i];

      checkInfactibility ();

      call_fn ();
      if (ncon > 0)
        newnormc = c->norm ();
      else
        newnormc = 0;

      if ( (SteihFlag == 1) || (SteihFlag == 2) )
        normd = DeltaH;
      else
        normd = d.norm ();

      if ( first && 
         ( (newnormc > Min (zeta1*rho, zeta1*normc + zeta2*rho) ) || 
           ( (normc <= csic) && 
             (newnormc > Max (csic, zeta3*normc) ) ) ) ) {

        StepFlag = NAstep (*c, ssoc);
        scale_xc (ssoc);

        // Arrumar tamanho do ssoc a partir do x
        Real alphassoc = 1;
        pReal ssocx = ssoc.get_doublex();
        for (Int i = 0; i < nvar; i++) {
          Real xi = xx[i], bli = blx[i], bui = bux[i], di = ssocx[i];
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
          Real si = sx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], di = ssocx[nvar + i];
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

        ssoc.saxpy (d, 1);
        asoc = ssoc.norm (0);
        if (asoc > DeltaH)
          ssoc.scale (DeltaH/asoc);
        *x = *xc;
        if (Ineq)
          *s = *sc;
        for (Int i = 0; i < nvar; i++)
          xx[i] += ssocx[i];
        for (Int i = 0; i < nconI; i++) {
          sx[i] += ssocx[nvar + i];
        }
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
      } else if (DLH <= eta2*qd)
        DeltaH = Max (DeltaH, normd*alphaI);

      first = dciFalse;
      if (DeltaH < 1e-10)
        break;

      CurrentTime = getTime() - StartTime;
    }

    normc = newnormc;

  }
}
