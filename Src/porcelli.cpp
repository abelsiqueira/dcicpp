#include "interface.h"
//#include <cassert>
#include <cmath>

/* Porcelli
 *
 * This routine computes a trust region vertical step. This step is
 * the solution of the problem:
 *
 * min 0.5*(||h(x) + A*d||^2 - ||h(x)||^2)
 * s.t. ||d|| <= DeltaV
 *
 * The algorithm uses a variant of the dogleg method described by
 * M. Powell in "A hybrid method for nonlinear equations", in
 * Numerical Methods for Nonlinear Algebraic Equation,
 * P. Rabinowitz, ed. 1970.
 */

namespace DCI {
  /* This function must (approximately) solve
   * min m(d) = 0.5*norm */
  Int Interface::LeastSquareTrustRegion (Vector & d) {
    d.reset(nvar + nconI, 0.0);
    Int nLstSqrs = 0, maxLstSqrs = nvar + nconI;
    Real theta, theta0, thetanew, alpha, beta, gamma;
    Vector lsGrad(*env);
    Real one[2] = {1,0}, zero[2] = {0,0}, mone[2] = {-1,0};
    Vector r(*env), p(*env), q(*env), dnew(*env);
    Real gtd, dtq, gtp, ptp, dtd, dtdnew, delta2, dtp, qd = 0;

    delta2 = DeltaV*DeltaV;

    lsGrad.sdmult(*J, 1, mone, zero, *c);
    r = lsGrad;
    theta0 = r.dot(r);
    p = r;
    theta = theta0;
    gtd = 0;
    
    while ( (theta > eps2) && (theta > eps1*theta0) && 
            (nLstSqrs <= maxLstSqrs) && (CurrentTime < MaxTime) ) {
      q.sdmult(*J, 0, one, zero, p);
      q.sdmult(*J, 1, one, zero, q);
      gamma = p.dot(q);
      dtq = d.dot(q);
      gtp = lsGrad.dot(p);
      ptp = p.dot(p);

      if (gamma <= eps3*ptp) {
        // Almost singular matrix.Stops
        return 1;
      }
      alpha = theta/gamma;
      dnew = d;
      dnew.saxpy (p, alpha);
      dtdnew = dnew.dot(dnew);

      if (dtdnew > delta2) {
        // Out of the region
        // Truncate step and stops
        dtp = d.dot(p);
        Real root1 = (-dtp + sqrt(dtp*dtp + (delta2 - dtd)*ptp))/ptp;

        d.saxpy (p, root1);
        gtd += root1*gtp;
        qd += root1*dtq + 0.5*root1*root1*gamma + root1*gtp;

        return 2;
      }

      if (alpha < dciTiny)
        return -2;

      r.saxpy (q, -alpha);

      thetanew = r.dot(r);
      beta = thetanew/theta;
      qd = qd + alpha + dtq + 0.5*alpha*alpha*gamma + alpha*gtp;
      p.scale(beta);
      p.saxpy(r,1);
      d = dnew;
      gtd += alpha*gtp;
      dtd = dtdnew;
      theta = thetanew;

      CurrentTime = getTime() - StartTime;
    }

    return 0;
  }
  
  Int Interface::Porcelli () {
    Real oldnormc = c->norm();
    Vector d(*env), dcp(*env), dn(*env);
    Real ndcp, ndn, normd;
    Real alpha, Ared, Pred, oldDelta;
    Real one[2] = {1,0};
    Real zero[2] = {0,0};
    Int iout, TrustIter, naflag;
    Bool dnavail;
    Vector gtmp (*env); 
    Vector xtmp (*xc), ctmp (*c), stmp (*env);
    Real normgtmp = 0;

    if (Ineq)
      stmp = *sc;

    //Remove later if needed
    call_ccfsg (dciTrue, ScaleVertical);
    Aavail = dciFalse;

    //This method uses the Porcelli scale matrix
    Real scalingMatrix[nvar + nconI];
    Vector Diag(*env);
    Diag.reset(nvar + nconI, 1.0);
    if (ScaleVertical)
      scale_xc(Diag);
    pReal Diagx = Diag.get_doublex();

    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
//    scale_xc(gtmp);
    pReal gtmpx = gtmp.get_doublex();
    for (Int i = 0; i < nvar; i++) {
      Real gi = gtmpx[i], zi = xcx[i], ui = bux[i], li = blx[i];
      if ( (gi < 0) && (ui < dciInf) ) {
        scalingMatrix[i] = ui - zi;
      } else if ( (gi > 0) && (li > -dciInf) ) {
        scalingMatrix[i] = zi - li;
      } else {
        scalingMatrix[i] = 1;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real gi = gtmpx[j], zi = scx[i], ui = cux[ineqIdx[i]], li = clx[ineqIdx[i]];
      if ( (gi < 0) && (ui < dciInf) ) {
        scalingMatrix[j] = ui - zi;
      } else if ( (gi > 0) && (li > -dciInf) ) {
        scalingMatrix[j] = zi - li;
      } else {
        scalingMatrix[j] = 1;
      }
    }
/*     if (penal_trust) {
 *       for (Int i = 0; i < nvar; i++) {
 *         Real val = 0;
 *         if ( (bux[i] < dciInf) && (blx[i] > -dciInf) ) {
 *           if (PartialPenal) {
 *             if ( (xcx[i] - blx[i]) < (bux[i] - xcx[i]) ) {
 *               val = 1;
 *             } else {
 *               val = -1;
 *             }
 *           } else {
 *             val = bux[i] + blx[i] - 2*xcx[i];
 *           }
 *         } else if (bux[i] < dciInf) {
 *           val = -1;
 *         } else if (blx[i] > -dciInf) {
 *           val = 1;
 *         }
 *         gtmpx[i] -= mu*val;
 *       }
 *       for (Int i = 0; i < nconI; i++) {
 *         Real val = 0;
 *         if ( (cux[ineqIdx[i]] < dciInf) && (clx[ineqIdx[i]] > -dciInf) ) {
 *           if (PartialPenal) {
 *             if ( (scx[i] - clx[ineqIdx[i]]) < (cux[ineqIdx[i]] - scx[i]) ) {
 *               val = 1;
 *             } else {
 *               val = -1;
 *             }
 *           } else {
 *             val = cux[ineqIdx[i]] + clx[ineqIdx[i]] - 2*scx[i];
 *           }
 *         } else if (cux[ineqIdx[i]] < dciInf)
 *           val = -1;
 *         else if (clx[ineqIdx[i]] > -dciInf)
 *           val = 1;
 *         gtmpx[nvar + i] -= mu*val;
 *       }
 *     }
 */
    normgtmp = gtmp.norm ();
//    DeltaV = normgtmp;

    if (normgtmp < dciTiny) {
      normc = oldnormc;
      iout = 6;
//      std::cout << "iout = 6" << std::endl;
      return iout;
    }
    //Now with the infinity norm
    Real lower[nvar + nconI], upper[nvar + nconI];
    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], li = blx[i], ui = bux[i];
      lower[i] = (li - zi) * (1 - epsmu)/Diagx[i];
      upper[i] = (ui - zi) * (1 - epsmu)/Diagx[i];
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      lower[j] = (li - zi) * (1 - epsmu)/Diagx[j];
      upper[j] = (ui - zi) * (1 - epsmu)/Diagx[j];
    }
    Vector aux(*env);
    aux.sdmult(*J, 0, one, zero, gtmp);
    d = gtmp;
    pReal dx = 0;
    dx = d.get_doublex();
    for (Int i = 0; i < nvar + nconI; i++) {
      dx[i] *= scalingMatrix[i];
    }
    alpha = normgtmp/aux.norm();
    alpha *= alpha;
    if (!ScaleVertical) {
      assert(0);
    } else {
      for (int i = 0; i < nvar + nconI; i++) {
        Real di = dx[i], ui = upper[i], li = lower[i];
        if (di > 0) {
          alpha = Min(alpha, ui/di);
        } else if (di < 0) {
          alpha = Min(alpha, li/di);
        }
      }
    }
    dcp.scale (gtmp, -alpha);
    pReal dcpx = dcp.get_doublex();
    for (Int i = 0; i < nvar + nconI; i++)
      dcpx[i] *= scalingMatrix[i];
//    if (ScaleVertical) scale_xc (dcp);
    ndcp = dcp.norm(0);

    dnavail = dciFalse;
    ndn = 0;
    Ared = 0;
    Pred = 1;
    oldDelta = DeltaV;
    DeltaV = DeltaV/kappa2;
    normd = DeltaV;
    iout = 0;
    TrustIter = 0;

    // For the inequalities
    Real dcpAlphamu = 1, dnAlphamu = 1;
    pReal dnx = 0;
    Real smlAlphamu = 1e-3;

    //Encontrar dn
    //Ver qual eh melhor
    while ( (Ared < kappa1*Pred) && (Aavail || (TrustIter < 2) ) && (CurrentTime < MaxTime) ) {
      TrustIter++;
      DeltaV = kappa2*normd;
//      DeltaV = Min(kappa2*normd, 0.9*DeltaV);

      // Cauchy step is inside trust region
      // sc + dcps >= epsmu*sc
      dnavail = dciFalse;
      if (!dnavail) {
        naflag = LeastSquareTrustRegion (dn);
//        naflag = NAstep (ctmp, dn); 
        dnavail = dciTrue;
//        if (ScaleVertical) scale_xc (dn);

        dnx = dn.get_doublex();
        //Project this step
        for (Int i = 0; i < nvar + nconI; i++) {
          if (dnx[i] > upper[i])
            dnx[i] = upper[i];
          else if (dnx[i] < lower[i])
            dnx[i] = lower[i];
        }

        if (naflag > 1)
          ndn = 0;
        else
          ndn = dn.norm (0);

        if (ndn > DeltaV)
          dn.scale(DeltaV/ndn);
      }

      Real newtonReduction, cauchyReduction;

      *xc = xtmp;
      if (Ineq)
        *sc = stmp;
      for (Int i = 0; i < nvar; i++)
        xcx[i] += Diagx[i]*(dcpx[i] - dnx[i]);
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        scx[i] += Diagx[j]*(dcpx[j] - dnx[j]);
      }
      call_ccfsg(dciFalse);
      cauchyReduction = oldnormc - c->norm();
      for (Int i = 0; i < nvar; i++)
        xcx[i] += Diagx[i]*dnx[i];
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        scx[i] += Diagx[j]*dnx[j];
      }
      call_ccfsg(dciFalse);
      newtonReduction = oldnormc - c->norm();
      Real factor = 1.0;
      while (newtonReduction/cauchyReduction < 0.5) {
        // Line search among from newton to cauchy
        pReal xtmpx = xtmp.get_doublex(), stmpx = 0;
        if (Ineq)
          stmpx = stmp.get_doublex();
        factor *= 0.9;
        for (Int i = 0; i < nvar; i++)
          xcx[i] = xtmpx[i] + Diagx[i]*(factor*dcpx[i] + (1-factor)*dnx[i]);
        for (Int i = 0; i < nconI; i++) {
          Int j = nvar + i;
          scx[i] = stmpx[i] + Diagx[j]*(factor*dcpx[j] + (1-factor)*dnx[j]);
        }
        call_ccfsg(dciFalse);
        newtonReduction = oldnormc - c->norm();
      }
      
      normc = c->norm ();

      Ared = oldnormc*oldnormc - normc*normc;
      if (iout == 2)
        Pred = oldnormc*oldnormc;
      else {
        gtmp.sdmult (*J, 0, one, zero, dn); // J*d
        Pred = gtmp.norm();
        Pred = - Pred*Pred - 2*gtmp.dot (ctmp);
      }

      CurrentTime = getTime() - StartTime;

    }

    if ( (Ared < kappa1*Pred) && (!Aavail) ) {
//      The algorithm has failed. A must be recomputed.

      *xc = xtmp;
      if (Ineq)
        *sc = stmp;
      *c = ctmp;
      normc = oldnormc;
      DeltaV = oldDelta;
      iout = 5;
    } else if (Ared >= kappa3*Pred)
      DeltaV = Max (kappa4*normd, DeltaV);

//    std::cout << "iout = " << iout << std::endl;
    return iout;

  }

}
