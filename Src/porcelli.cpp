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
  Int Interface::LeastSquareTrustRegion (Vector & d, pReal scalingMatrix) {
    d.reset(nvar + nconI, 0.0);
    Int nLstSqrs = 0, maxLstSqrs = nvar + nconI;
    Real theta, theta0, thetanew, alpha, beta, gamma;
    Vector lsGrad(*env);
    Real one[2] = {1,0}, zero[2] = {0,0}, mone[2] = {-1,0};
    Vector r(*env), p(*env), q(*env), dnew(*env), t(*env);
    Real gtd = 0, dtq = 0, gtp = 0, ptp = 0, dtd = 0, dtdnew = 0, delta2, dtp = 0, qd = 0;

    delta2 = DeltaV*DeltaV;

    lsGrad.sdmult(*J, 1, one, zero, *c);
    r.scale(lsGrad, -1);
    theta0 = r.dot(r);
    p = r;
    theta = theta0;
    gtd = 0;
    Real normGrad0 = lsGrad.norm();
    Real normGrad = normGrad0;
    
    while ( (theta > eps2) && (theta > eps1*theta0) && 
            (normGrad > 0.2*normGrad0) &&
            (nLstSqrs <= maxLstSqrs) && (CurrentTime < MaxTime) ) {
      q.sdmult(*J, 0, one, zero, p);
      q.sdmult(*J, 1, one, zero, q);
      gamma = p.dot(q);
      dtq = d.dot(q);
      gtp = lsGrad.dot(p);
      pReal dx = d.get_doublex(), px = p.get_doublex();
      ptp = p.dot(p);
      ptp = 0.0;
      for (Int i = 0; i < nvar + nconI; i++)
        ptp += pow(px[i]*scalingMatrix[i], 2);


      if (gamma <= eps3*ptp) {
        // Almost singular matrix.Stops
        return 1;
      }
      alpha = theta/gamma;
      dnew = d;
      dnew.saxpy (p, alpha);
      dtdnew = dnew.dot(dnew);
      dtdnew = 0.0;
      pReal dnewx = dnew.get_doublex();
      for (Int i = 0; i < nvar + nconI; i++)
        dtdnew += pow(dnewx[i]*scalingMatrix[i], 2);

      if (dtdnew > delta2) {
        // Out of the region
        // Truncate step and stops
        dtp = d.dot(p);
        dtp = 0;
        for (Int i = 0; i < nvar + nconI; i++)
          dtp += dx[i]*pow(scalingMatrix[i], 2)*px[i];
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

      lsGrad = *c;
      lsGrad.sdmult(*J, 0, one, one, d);
      lsGrad.sdmult(*J, 1, one, zero, lsGrad);
      normGrad = lsGrad.norm();

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
    Real beta1 = 0.1, beta2 = 0.25;

    if (Ineq)
      stmp = *sc;

    //Remove later if needed
//    call_ccfsg_xc (dciTrue, dciFalse);
    call_ccfsg_xc (dciTrue, ScaleVertical);
    Aavail = dciFalse;

    //This method uses the Porcelli scale matrix
    Real scalingMatrix[nvar + nconI];
    Vector Diag(*env);
    Diag.reset(nvar + nconI, 1.0);
    if (ScaleVertical) scale_xc(Diag);
    pReal Diagx = Diag.get_doublex();

    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
//    scale_xc(gtmp);
    pReal gtmpx = gtmp.get_doublex();
    for (Int i = 0; i < nvar; i++) {
      Real gi = gtmpx[i], zi = xcx[i], ui = bux[i], li = blx[i];
      if ( (gi < 0) && (ui < dciInf) ) {
        scalingMatrix[i] = 1.0/sqrt(ui - zi);
      } else if ( (gi > 0) && (li > -dciInf) ) {
        scalingMatrix[i] = 1.0/sqrt(zi - li);
      } else {
        scalingMatrix[i] = 1;
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real gi = gtmpx[j], zi = scx[i], ui = cux[ineqIdx[i]], li = clx[ineqIdx[i]];
      if ( (gi < 0) && (ui < dciInf) ) {
        scalingMatrix[j] = 1.0/sqrt(ui - zi);
      } else if ( (gi > 0) && (li > -dciInf) ) {
        scalingMatrix[j] = 1.0/sqrt(zi - li);
      } else {
        scalingMatrix[j] = 1;
      }
    }
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
      lower[i] = (li > -dciInf ? (li - zi) * (1 - epsmu)/Diagx[i] : -dciInf);
      upper[i] = (ui < dciInf ? (ui - zi) * (1 - epsmu)/Diagx[i] : dciInf);
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      lower[j] = (li > -dciInf ? (li - zi) * (1 - epsmu)/Diagx[j] : -dciInf);
      upper[j] = (ui < dciInf ? (ui - zi) * (1 - epsmu)/Diagx[j] : dciInf);
    }
    Vector aux(*env);
    d = gtmp;
    pReal dx = 0; 
    gtmpx = gtmp.get_doublex();
    dx = d.get_doublex();
    for (Int i = 0; i < nvar + nconI; i++) {
      gtmpx[i] /= scalingMatrix[i];
      dx[i] = -dx[i]/pow(scalingMatrix[i], 2);
    }

    if (ScaleVertical) scale_xc (d);

    aux.sdmult(*J, 0, one, zero, d);
    alpha = Min(gtmp.dot(gtmp)/aux.dot(aux), DeltaV/gtmp.norm());
    dcp.scale (d, alpha);
    pReal dcpx = dcp.get_doublex();

//    alpha = -d.dot(gtmp)/aux.dot(aux);
//    alpha = Min(alpha, DeltaV/d.norm());
    alpha = 1.0;
    for (int i = 0; i < nvar + nconI; i++) {
      Real di = dcpx[i], ui = upper[i], li = lower[i];
      if (di > 0) {
        alpha = Min(alpha, ui/(Diagx[i]*di));
      } else if (di < 0) {
        alpha = Min(alpha, li/(Diagx[i]*di));
      }
    }
    Real theta = 0.99995;
    if (alpha < 1) {
      alpha = Max(theta, 1 - dcp.norm())*alpha;
      dcp.scale(alpha);
    }
/*     for (Int i = 0; i < nvar + nconI; i++)
 *       dcpx[i] *= scalingMatrix[i];
 */
    if (ScaleVertical) scale_xc (dcp);
    ndcp = dcp.norm();

    dnavail = dciFalse;
    ndn = 0;
    Ared = 0;
    Pred = 1;
    oldDelta = DeltaV;
//    DeltaV = DeltaV/kappa2;
    normd = DeltaV;
    iout = 0;
    TrustIter = 0;

    // For the inequalities
    pReal dnx = 0;

    //Encontrar dn
    //Ver qual eh melhor
//    while ( (Ared < kappa1*Pred) && (Aavail || (TrustIter < 20) ) && (CurrentTime < MaxTime) ) {
//    while ( (Ared < kappa1*Pred) && (TrustIter < 100000) ) {
//      TrustIter++;
//      DeltaV = kappa2*normd;
//      DeltaV = Min(kappa2*normd, 0.9*DeltaV);

      // Cauchy step is inside trust region
      // sc + dcps >= epsmu*sc
      dnavail = dciFalse;
      if (!dnavail) {
        naflag = LeastSquareTrustRegion (dn, scalingMatrix);
/*         naflag = NAstep (ctmp, dn); 
 *         if (dn.norm() > DeltaV) {
 *           dn.scale(DeltaV/dn.norm());
 *         }
 */

        dnavail = dciTrue;
        if (ScaleVertical) scale_xc (dn);

        dnx = dn.get_doublex();
        //Project this step
        alpha = 1.0;
        for (Int i = 0; i < nvar + nconI; i++) {
          Real di = dnx[i], ui = upper[i], li = lower[i];
          if (di > 0) {
            alpha = Min(alpha, ui/(Diagx[i]*di));
          } else if (di < 0) {
            alpha = Min(alpha, li/(Diagx[i]*di));
          }
        }
        if (alpha < 1) {
          alpha = Max(theta, 1 - dn.norm())*alpha;
          dn.scale(alpha);
        }

        if (naflag > 1)
          ndn = 0;
        else
          ndn = dn.norm ();
        assert(ndn <= DeltaV || "ndn > DeltaV");
      }

      /* ||a + b||^2 = <a+b,a+b> = <a,a> + 2*<a,b> + <b,b> */
      /* m(d) = 0.5*||J*d + h||^2 
       * dtr = t*dn + (1 - t)*dcp 
       * m(dtr) = 0.5*||J*(t*dn + (1-t)*dcp) + h||^2 
       *   = 0.5*||J*dcp + h + t*J*(dn - dcp)||^2 
       *   = 0.5*||J*dcp + h||^2 + t*(J*dcp + h)'*J*(dn - dcp) + 0.5*t^2*||J*(dn - dcp)||^2 */
      Vector Adcph(*c), difdcdn(dn), Adif(*env);
      Adcph.sdmult(*J, 0, one, one, dcp);
      difdcdn.saxpy(dcp, -1);
      Adif.sdmult(*J, 0, one, zero, difdcdn);

      Real objValAdcph = 0.5*Adcph.dot(Adcph);
      Real dotAdcphAdif = Adcph.dot(Adif);
      Real halfSqrNormAdif = 0.5*Adif.dot(Adif);

      Real cauchyReduction = 0.5*oldnormc*oldnormc - objValAdcph;
      Real newtonReduction = cauchyReduction - dotAdcphAdif - halfSqrNormAdif;

      Real factor = 1.0;
      while (newtonReduction/cauchyReduction < beta1) {
        // Line search among from newton to cauchy
        factor *= 0.9;
        newtonReduction = cauchyReduction - factor*dotAdcphAdif - pow(factor,2)*halfSqrNormAdif;
      }

      Vector xtemp(*xc), stemp((nconI ? *sc : *env));

      for (Int i = 0; i < nvar; i++)
        xcx[i] += Diagx[i]*(factor*dnx[i] + (1 - factor)*dcpx[i]);
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        scx[i] += Diagx[j]*(factor*dnx[j] + (1 - factor)*dcpx[j]);
      }
      
      checkInfactibility();

      call_ccfsg_xc(dciFalse);
      normc = c->norm ();

      Ared = 0.5*(oldnormc*oldnormc - normc*normc);
      Pred = newtonReduction;

      if (Ared/Pred < beta2) {
        DeltaV /= 4;
        *xc = xtmp;
        if (Ineq) *sc = stmp;
        call_ccfsg_xc(dciFalse);
        normc = c->norm();
//        std::cout << "Porcelli: Bad point" << std::endl;
      } else if (Ared/Pred > 0.75) {
        DeltaV *= 2;
      }

      if (normc < rho)
        return 0;

      CurrentTime = getTime() - StartTime;

//    }

  }

}
