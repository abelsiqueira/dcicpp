#include "interface.h"
//#include <cassert>
#include <cmath>

/* This routine computes the horizontal step by using Steihaug's algorithm.
 *
 * The horizontal step is an approximate solution for:
 *
 * min 0.5*d'*B*d + g_p'*d
 * s.t. A*d = 0
 *      ||d|| < delta
 */

namespace DCI {
  Int Interface::dcisteih (Vector & d, Real & qd, Real & gtd) {
    Real delta2, dtd, theta0, theta, thetanew, alpha, beta, gamma;
    Real dtq, gtp, dtp, ptp, root1, root2, dtdnew, qdnew;
    Int SteihFlag;
    Vector r(*env), p(d), q(d), dnew(*env), v(*env), tmp(*env);
    Real lower[nvar + nconI], upper[nvar + nconI];

    pReal qx = q.get_doublex (), px = p.get_doublex ();
    pReal dx = d.get_doublex();

    for (Int i = 0; i < nvar + nconI; i++)
      dx[i] = 0;

    if (DeltaH < 1e-12)
      return -1;
    delta2 = DeltaH*DeltaH;
    qd = 0;
    dtd = 0;
    r.scale (*gp, -1);
    theta0 = r.dot (r);
    p.scale (r, 1);
    theta = theta0;
    gtd = 0;
    nSteih = 0;
//    std::vector<Real> tmpDiag (nvar + nconI, 1);
//    Vector Diag(*env, tmpDiag);
//    pReal Diagx = Diag.get_doublex();
//    scale_xc (Diag);
    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], bli = blx[i], bui = bux[i];
//      lower[i] = (bli - xi) * (1 - epsmu) / Diagx[i];
//      upper[i] = (bui - xi) * (1 - epsmu) / Diagx[i];
      lower[i] = (bli - xi) * (1 - epsmu);
      upper[i] = (bui - xi) * (1 - epsmu);
      if ( (PartialPenal) && (bli > -dciInf) && (bui < dciInf) ) {
        if ( (xi - bli) < (bui - xi) ) {
          lower[i] = lower[i]/(xi - bli);
          upper[i] = upper[i]/(xi - bli);
        } else {
          lower[i] = lower[i]/(bui - xi);
          upper[i] = upper[i]/(bui - xi);
        }
        continue;
      }
      if (bli > -dciInf) {
        lower[i] = lower[i]/(xi - bli);
        upper[i] = upper[i]/(xi - bli);
      }
      if (bui < dciInf) {
        lower[i] = lower[i]/(bui - xi);
        upper[i] = upper[i]/(bui - xi);
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Real si = scx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
      Int j = nvar + i;
//      lower[j] = (cli - si) * (1 - epsmu) / Diagx[j];
//      upper[j] = (cui - si) * (1 - epsmu) / Diagx[j];
      lower[j] = (cli - si) * (1 - epsmu);
      upper[j] = (cui - si) * (1 - epsmu);
      if ( (PartialPenal) && (cli > -dciInf) && (cui < dciInf) ) {
        if ( (si - cli) < (cui - si) ) {
          lower[j] = lower[j]/(si - cli);
          upper[j] = upper[j]/(si - cli);
        } else {
          lower[j] = lower[j]/(cui - si);
          upper[j] = upper[j]/(cui - si);
        }
        continue;
      }
      if (cli > -dciInf) {
        lower[j] = lower[j]/(si - cli);
        upper[j] = upper[j]/(si - cli);
      }
      if (cui < dciInf) {
        lower[j] = lower[j]/(cui - si);
        upper[j] = upper[j]/(cui - si);
      }
    }


    while ( (theta > eps2) && (theta > eps1*theta0) && (nSteih <= maxitSteih) && (CurrentTime < MaxTime) ) {
      nSteih++;
      call_prod_xc (GotH, px, qx);
      nHprod++;
      gamma = p.dot (q);
      dtq = d.dot (q);
      gtp = g->dot (p);
      ptp = p.dot (p);

      if (gamma <= eps3*ptp) {
        // Negative Curvature
        dtp = d.dot (p);
        root1 = sqrt (dtp*dtp + (delta2 - dtd)*ptp);
        root2 = (-dtp - root1)/ptp;
        root1 = (-dtp + root1)/ptp;

        Real root1max = 1, root2max = 1;
        for (Int i = 0; i < nvar + nconI; i++) {
          Real pxi = px[i], lowi = lower[i], uppi = upper[i];
          if (pxi == 0)
            continue;
          pxi = root1*pxi;
          if (pxi < 0)
            root1max = Min (root1max, lowi/pxi);
          else
            root1max = Min (root1max, uppi/pxi);
          pxi = root2*px[i];
          if (pxi < 0)
            root2max = Min (root2max, lowi/pxi);
          else
            root2max = Min (root2max, uppi/pxi);
        }

        root1 *= root1max;
        root2 *= root2max;

        dnew.scale (d, 1);
        dnew.saxpy (p, root2);
        qdnew = qd + root2*dtq + 0.5 * root2*root2*gamma + root2*gtp;
        d.saxpy (p, root1);
        qd = qd + root1*dtq + 0.5 * root1*root1*gamma + root1*gtp;

        if (qdnew < qd) {
          d.scale (dnew, 1);
          qd = qdnew;
          gtd = gtd + root2*gtp;
        } else
          gtd = gtd + root1*gtp;

        SteihFlag = 1;
        return SteihFlag;
      }

      alpha = theta/gamma;
      dnew.scale (d, 1);
      dnew.saxpy (p, alpha);
      dtdnew = dnew.dot (dnew);

      if (dtdnew > delta2) {
        dtp = d.dot (p);
        root1 = (-dtp + sqrt(dtp*dtp + (delta2 - dtd)*ptp))/ptp;

        Real root1max = 1;
        for (Int i = 0; i < nvar + nconI; i++) {
          Real pxi = px[i], lowi = lower[i], uppi = upper[i];
          pxi = pxi*root1;
          if (pxi == 0)
            continue;
          if (pxi < 0)
            root1max = Min (root1max, lowi/pxi);
          else
            root1max = Min (root1max, uppi/pxi);
        }

        root1 *= root1max;

        d.saxpy (p, root1);
        gtd = gtd + root1*gtp;
        qd = qd + root1*dtq + 0.5*root1*root1*gamma + root1*gtp;
        SteihFlag = 2;
        return SteihFlag;
      }

      if (alpha < dciTiny)
        return -2;

      Real alphamax = 1;
      for (Int i = 0; i < nvar + nconI; i++) {
        Real pxi = px[i], lowi = lower[i], uppi = upper[i];
        if (pxi == 0)
          continue;
        pxi = alpha*pxi;
        if (pxi < 0)
          alphamax = Min (alphamax, lowi/pxi);
        else
          alphamax = Min (alphamax, uppi/pxi);
      }

      if (alphamax < 1) {
        alpha *= alphamax;
        d.saxpy (p, alpha);
        qd = qd + alpha*dtq + 0.5*alpha*alpha*gamma + alpha*gtp;
        SteihFlag = 10;
        return SteihFlag;
      }

      for (Int i = 0; i < nvar + nconI; i++) {
        lower[i] -= alpha * px[i];
        upper[i] -= alpha * px[i];
      }

      r.saxpy (q, -alpha);
      v.scale (r, 1);

      LimLbd = dciFalse;
      LbdMax = 1e20;
      if (ncon > 0)
        NAproj (v, r, tmp);
      else
        r.scale (v, 1);

      thetanew = r.dot (r);
      beta = thetanew/theta;
      qd = qd + alpha * dtq + 0.5*alpha*alpha*gamma + alpha*gtp;
      p.scale (beta);
      p.saxpy (r, 1);
      d.scale (dnew, 1);
      gtd = gtd + alpha*gtp;
      dtd = dtdnew;
      theta = thetanew;

      CurrentTime = getTime() - StartTime;
    }

    if (nSteih > maxitSteih)
      SteihFlag = 3;
    else
      SteihFlag = 0;

    return SteihFlag;
  }
}
