#include "interface.h"
//#include <cassert>
#include <cmath>

/* innerVerticalStep
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
  Int Interface::innerVerticalStep (Real & infeasible_gradient) {
    Real oldnormc = c->norm();
    Vector d(*env), dcp(*env), dn(*env);
    Real ndn;
    Real alpha, Ared, Pred;
    Real one[2] = {1,0};
    Real zero[2] = {0,0};
    Int iout, naflag;
    Bool dnavail;
    Vector gtmp (*env); 
    Vector xtmp (*xc), ctmp (*c), stmp (*env);
    Real normgtmp = 0;
    Real beta1 = 0.1, beta2 = 0.25;

    Aavail = dciFalse;

    Real scalingMatrix[nvar + nconI];
    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
    pReal gtmpx = gtmp.get_doublex();
    for (Int i = 0; i < nvar+nconI; i++) {
      Real gi = gtmpx[i], zi = xcx[i], ui = u_bndx[i], li = l_bndx[i];
      if ( (gi < 0) && (ui < dciInf) ) {
        scalingMatrix[i] = 1.0/sqrt(ui - zi);
      } else if ( (gi > 0) && (li > -dciInf) ) {
        scalingMatrix[i] = 1.0/sqrt(zi - li);
      } else {
        scalingMatrix[i] = 1;
      }
    }
    normgtmp = gtmp.norm ();
    Vector gtmp_proj(*env);
    gtmp_proj.scale(gtmp, -1.0);
    projectBounds_xc(gtmp_proj);
    infeasible_gradient = gtmp_proj.norm();
//    DeltaV = normgtmp;

    if (normgtmp < dciTiny) {
      normc = oldnormc;
      iout = 6;
//      std::cout << "iout = 6" << std::endl;
      return iout;
    }
    //Now with the infinity norm
    Real lower[nvar + nconI], upper[nvar + nconI];
    for (Int i = 0; i < nvar+nconI; i++) {
      Real zi = xcx[i], li = l_bndx[i], ui = u_bndx[i];
      lower[i] = Max( -DeltaV, (li > -dciInf ? (li - zi) * (1 - epsmu) :
            -dciInf) );
      upper[i] = Min(  DeltaV, (ui <  dciInf ? (ui - zi) * (1 - epsmu) : 
             dciInf) );
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

    aux.sdmult(*J, 0, one, zero, d);
    alpha = Min(gtmp.dot(gtmp)/aux.dot(aux), DeltaV/gtmp.norm());
    dcp.scale (d, alpha);
    pReal dcpx = dcp.get_doublex();

//    alpha = -d.dot(gtmp)/aux.dot(aux);
//    alpha = Min(alpha, DeltaV/d.norm());
    alpha = 1.0;
    for (int i = 0; i < nvar + nconI; i++) {
      Real di = dcpx[i], ui = upper[i], li = lower[i];
      if (li - ui > -dciEps) {
        dcpx[i] = 0;
        continue;
      }
      if (di > dciEps) {
        alpha = Min(alpha, ui/(di));
      } else if (di < -dciEps) {
        alpha = Min(alpha, li/(di));
      } else
        dcpx[i] = 0;
    }
    Real theta = 0.99995;
    if (alpha < 1) {
      alpha = Max(theta, 1 - dcp.norm())*alpha;
      dcp.scale(alpha);
    }
/*     for (Int i = 0; i < nvar + nconI; i++)
 *       dcpx[i] *= scalingMatrix[i];
 */
    dnavail = dciFalse;
    ndn = 0;
    Ared = 0;
    Pred = 1;
//    DeltaV = DeltaV/kappa2;
    iout = 0;

    // For the inequalities
    pReal dnx = 0;

    //Encontrar dn
    //Ver qual eh melhor
//    while ( (Ared < kappa1*Pred) && (Aavail || (TrustIter < 20) ) && (current_time < max_time) ) {
//    while ( (Ared < kappa1*Pred) && (TrustIter < 100000) ) {
//      TrustIter++;

//    Int iter = 0;

//    while (normc > rho && iter < 2) {
      // Cauchy step is inside trust region
      // sc + dcps >= epsmu*sc
      dnavail = dciFalse;
      if (!dnavail) {
//        naflag = leastSquaresTrustRegion (dn, scalingMatrix, lower, upper);
        naflag = naStep (*c, dn);
/*         naflag = naStep (ctmp, dn); 
 *         if (dn.norm() > DeltaV) {
 *           dn.scale(DeltaV/dn.norm());
 *         }
 */

        dnavail = dciTrue;

        dnx = dn.get_doublex();
        //Project this step
        alpha = 1.0;
        for (Int i = 0; i < nvar + nconI; i++) {
          Real di = dnx[i], ui = upper[i], li = lower[i];
          if (li - ui > -dciEps) {
            dnx[i] = 0;
            continue;
          }
          if (di > dciEps) {
            alpha = Min(alpha, ui/di);
          } else if (di < -dciEps) {
            alpha = Min(alpha, li/di);
          } else
            di = 0;
        }
        if (alpha < 1) {
          alpha = Max(theta, 1 - dn.norm())*alpha;
          dn.scale(alpha);
        }

        if (naflag > 1)
          ndn = 0;
        else
          ndn = dn.norm (0);
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
        if (factor < 1e-8) {
          factor = 0;
          break;
        }
      }

      Vector xtemp(*xc);

      for (Int i = 0; i < nvar+nconI; i++) {
        if (l_bndx[i] - u_bndx[i] > -dciEps)
          continue;
        xcx[i] += (factor*dnx[i] + (1 - factor)*dcpx[i]);
        if (xcx[i] >= u_bndx[i])
          xcx[i] = u_bndx[i] - dciEps;
        else if (xcx[i] <= l_bndx[i])
          xcx[i] = l_bndx[i] + dciEps;
      }
      
#ifndef NDEBUG
      checkInfactibility();
#endif

      call_ccfsg_xc(dciFalse);
      normc = c->norm ();

      Ared = 0.5*(oldnormc*oldnormc - normc*normc);
      Pred = newtonReduction;

      if (Ared/Pred < beta2) {
        DeltaV /= 4;
        *xc = xtmp;
        call_ccfsg_xc(dciFalse);
        normc = c->norm();
      } else if (Ared/Pred > 0.75) {
        DeltaV *= 2;
      }

      if (normc < rho)
        return 0;

//      for (Int i = 0; i < nvar + nconI; i++) {
//        lower[i] -= factor*dnx[i] + (1 - factor) * dcpx[i];
//        upper[i] -= factor*dnx[i] + (1 - factor) * dcpx[i];
//      }

//      iter++;

      current_time = getTime() - start_time;

//    }

      return 0;
  }

}
