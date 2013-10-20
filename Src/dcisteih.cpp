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
  Int Interface::innerHorizontalStep (Vector & d, Real & qd, Real & gtd) {
    Real theta0, theta, thetanew, alpha, beta, gamma;
    Real dtHp, gtp, ptp, root1, root2, qdnew;
    Int SteihFlag;
    Vector r(*env), p(d), Hp(d), dnew(*env), v(*env), tmp(*env);
    Real lower[nvar + nconI], upper[nvar + nconI];

    pReal Hpx = Hp.get_doublex (), px = p.get_doublex ();
    pReal dx = d.get_doublex();

#ifdef VERBOSE
    if (VerboseLevel > 2) {
      std::cout << "DeltaH = " << DeltaH << std::endl;
    }
#endif
    for (Int i = 0; i < nvar + nconI; i++)
      dx[i] = 0;

    if (DeltaH < 1e-12)
      return -1;
    qd = 0;

    r.scale (*gp, -1);
    theta0 = r.dot (r);
    p = r;
    theta = theta0;
    gtd = 0;
    nSteih = 0;

    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], bli = blx[i], bui = bux[i];
      lower[i] = Max( (bli - xi) * (1 - epsmu)/Lambda[i], -DeltaH/Lambda[i] );
      upper[i] = Min( (bui - xi) * (1 - epsmu)/Lambda[i], DeltaH/Lambda[i] );
    }
    for (Int i = 0; i < nconI; i++) {
      Real si = scx[i], cli = clx[ineq_index[i]], cui = cux[ineq_index[i]];
      Int j = nvar + i;
      lower[j] = Max( (cli - si) * (1 - epsmu)/Lambda[j], -DeltaH/Lambda[j] );
      upper[j] = Min( (cui - si) * (1 - epsmu)/Lambda[j], DeltaH/Lambda[j] );
    }


    while ( (theta > eps2) && (theta > eps1*theta0) && (nSteih <= maxitSteih) && (CurrentTime < MaxTime) ) {
      nSteih++;
      call_prod_xc (GotH, px, Hpx);
      nHprod++;
      gamma = p.dot (Hp);
      dtHp = d.dot (Hp);
      gtp = g->dot (p);
      ptp = p.dot (p);
      ptp = 0.0;
      for (Int i = 0; i < nvar + nconI; i++) {
        ptp += pow(px[i]*Lambda[i], 2);
      }

      if (gamma <= eps3*ptp) {
        // Negative Curvature
        root1 = dciInf;
        root2 = dciInf;
        for (Int i = 0; i < nvar + nconI; i++) {
          if (px[i] > 0) {
            root1 = Min( root1, (upper[i] - dx[i])/px[i] );
            root2 = Min( root2, -(lower[i] - dx[i])/px[i] );
          } else if (px[i] < 0) {
            root1 = Min( root1, (lower[i] - dx[i])/px[i] );
            root2 = Min( root2, -(upper[i] - dx[i])/px[i] );
          }
        }

        root2 = -root2;
        dnew = d;
        dnew.saxpy (p, root2);
        qdnew = qd + root2*dtHp + 0.5 * root2*root2*gamma + root2*gtp;
        d.saxpy (p, root1);
        qd = qd + root1*dtHp + 0.5 * root1*root1*gamma + root1*gtp;

        if (qdnew < qd) {
          d = dnew;
          qd = qdnew;
          gtd = gtd + root2*gtp;
        } else
          gtd = gtd + root1*gtp;

        SteihFlag = 1;
        return SteihFlag;
      }

      alpha = theta/gamma;
      dnew = d;
      dnew.saxpy (p, alpha);
      pReal dnewx = dnew.get_doublex();
      bool outsideRegion = false;
      for (Int i = 0; i < nvar + nconI; i++)  {
        if ( (dnewx[i] > upper[i]) || (dnewx[i] < lower[i]) ) {
          outsideRegion = true;
          break;
        }
      }

      if (outsideRegion) {
        root1 = dciInf;
        for (Int i = 0; i < nvar + nconI; i++) {
          if (px[i] > 0) {
            root1 = Min( root1, (upper[i] - dx[i])/px[i] );
          } else if (px[i] < 0) {
            root1 = Min( root1, (lower[i] - dx[i])/px[i] );
          }
        }

        d.saxpy (p, root1);
        gtd = gtd + root1*gtp;
        qd = qd + root1*dtHp + 0.5*root1*root1*gamma + root1*gtp;
        SteihFlag = 2;
        return SteihFlag;
      }

      if (alpha < dciTiny)
        return -2;

      r.saxpy (Hp, -alpha);
      v = r;

      LimLbd = dciFalse;
      LbdMax = 1e20;
      if (ncon > 0)
        naProj (v, r, tmp);
      else
        r = v;

      thetanew = r.dot (r);
      beta = thetanew/theta;
      qd = qd + alpha * dtHp + 0.5*alpha*alpha*gamma + alpha*gtp;
      p.scale (beta);
      p.saxpy (r, 1);
      d = dnew;
      gtd = gtd + alpha*gtp;
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
