#include "interface.h"
//#include <cassert>
#include <cmath>

/* DCI Trust
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
  Int Interface::dcitrust (Real oldnormc) {
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
      stmp.scale(*sc, 1);


    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
    normgtmp = gtmp.norm ();

    if (normgtmp < dciTiny) {
      normc = oldnormc;
      iout = 6;
      return iout;
    }
    d.sdmult (*J, 0, one, zero, gtmp);
    alpha = normgtmp/d.norm();
    alpha *= alpha;
    alpha = Min (alpha, DeltaV/normgtmp);
    dcp.scale (gtmp, -alpha);
    scale_xc (dcp);
    ndcp = dcp.norm();


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
    pReal dcpx = dcp.get_doublex();
    pReal dx = 0;
    pReal dnx = 0;
    Real smlAlphamu = 1e-3;

    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], dcpi = dcpx[i], bli = blx[i], bui = bux[i];
      if (dcpi == 0)
        continue;
      if (dcpi < 0) {
        Real val = (bli - xi)*(1 - epsmu)/dcpi;
        if (val < 1)
          dcpx[i] *= val;
//        dcpAlphamu = Min (dcpAlphamu, val);
      } else {
        Real val = (bui - xi)*(1 - epsmu)/dcpi;
        if (val < 1)
          dcpx[i] *= val;
//        dcpAlphamu = Min (dcpAlphamu, val);
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real si = scx[i], dcpi = dcpx[j], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
      if (dcpi == 0)
        continue;
      if (dcpi < 0) {
        Real val = (cli - si)*(1 - epsmu)/dcpi;
        if (val < 1)
          dcpx[j] *= val;
//          dcpAlphamu = Min (dcpAlphamu, val);
      } else {
        Real val = (cui - si)*(1 - epsmu)/dcpi;
        if (val < 1)
          dcpx[j] *= val;
//          dcpAlphamu = Min (dcpAlphamu, val);
      }
    }

    while ( (Ared < kappa1*Pred) && (Aavail || (TrustIter < 2) ) && (CurrentTime < MaxTime) ) {
      TrustIter++;
      DeltaV = kappa2*normd;

      if ( (ndcp < DeltaV) && (dcpAlphamu == 1) ) {
        // Cauchy step is inside trust region
        // sc + dcps >= epsmu*sc
        if (!dnavail) {
          naflag = NAstep (ctmp, dn); 
          dnavail = dciTrue;
          scale_xc (dn);

          if (naflag > 1)
            ndn = 0;
          else
            ndn = dn.norm ();
          
          dnx = dn.get_doublex();
          for (Int i = 0; i < nvar; i++) {
            Real xi = xcx[i], dni = dnx[i], bli = blx[i], bui = bux[i];
            if (dni == 0)
              continue;
            if (dni < 0) {
              Real val = (bli - xi)*(1 - epsmu)/dni;
                dnAlphamu = Min (dnAlphamu, val);
            } else {
              Real val = (bui - xi)*(1 - epsmu)/dni;
                dnAlphamu = Min (dnAlphamu, val);
            }
          }
          for (Int i = 0; i < nconI; i++) {
            Real si = scx[i], dni = dnx[nvar + i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
            if (dni == 0)
              continue;
            if (dni < 0) {
              Real val = (cli - si)*(1 - epsmu)/dni;
                dnAlphamu = Min (dnAlphamu, val);
            } else {
              Real val = (cui - si)*(1 - epsmu)/dni;
                dnAlphamu = Min (dnAlphamu, val);
            }
          }
        }

        if (ndn <= ndcp) {
          // dn is too small. Use dcp.
          d.scale (dcp, 1);
          normd = ndcp;
          iout = 1;
        } else if ( (ndn <= DeltaV) && (dnAlphamu == 1) ) {
          d.scale (dn, 1);
          normd = ndn;
          iout = 2;
        } else if ( (ndn <= DeltaV) && (dnAlphamu < 1) ) {
          d.scale (dn, dnAlphamu);
          normd = ndn*dnAlphamu;
          iout = 10;
        } else if ( (ndn > DeltaV) && (dnAlphamu == 1) ) {
          Real dntdcp = dcp.dot(dn);
          Real nwsqr = ndcp*ndcp,
               wtv = dntdcp - nwsqr,
               nvsqr = ndn*ndn - 2*dntdcp + nwsqr;
          alpha = (-wtv + sqrt(wtv*wtv - nvsqr*(nwsqr - DeltaV*DeltaV)))/nvsqr;
          d.scale (dn, alpha);
          dx = d.get_doublex();
          alpha = 1 - alpha;
          for (Int i = 0; i < nvar + nconI; i++)
            dx[i] += alpha * dcpx[i];
          normd = DeltaV;
          iout = 3;
        } else {
          d.scale (dcp, 1);
          normd = ndcp;
          /*dn.scale (dnAlphamu);
          ndn *= dnAlphamu;
          Real dntdcp = dcp.dot(dn);
          Real nwsqr = ndcp*ndcp,
               wtv = dntdcp - nwsqr,
               nvsqr = ndn*ndn - 2*dntdcp + nwsqr;
          alpha = (-wtv + sqrt(wtv*wtv - nvsqr*(nwsqr - DeltaV*DeltaV)))/nvsqr;
          d.scale (dn, alpha);
          dx = d.get_doublex();
          alpha = 1 - alpha;
          for (Int i = 0; i < nvar + nconI; i++)
            dx[i] += alpha * dcpx[i];

          Real dAlphamu = 1;
          for (Int i = 0; i < nvar; i++) {
            Real xi = xcx[i], di = dx[i], bli = blx[i], bui = bux[i];
            if (di == 0)
              continue;
            if (di < 0) {
              Real val = (bli - xi)*(1 - epsmu)/di;
                dAlphamu = Min (dAlphamu, val);
            } else {
              Real val = (bui - xi)*(1 - epsmu)/di;
                dAlphamu = Min (dAlphamu, val);
            }
          }
          for (Int i = 0; i < nconI; i++) {
            Real si = scx[i], di = dx[nvar + i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
            if (di == 0)
              continue;
            if (di < 0) {
              Real val = (cli - si)*(1 - epsmu)/di;
                dAlphamu = Min (dAlphamu, val);
            } else {
              Real val = (cui - si)*(1 - epsmu)/di;
                dAlphamu = Min (dAlphamu, val);
            }
          }

          d.scale (dAlphamu);
          normd = dAlphamu*DeltaV;
          */
          iout = 11;
        }
      } else if ( (ndcp >= DeltaV) && (dcpAlphamu == 1) ) {
        d.scale (dcp, DeltaV/ndcp);
        normd = DeltaV;
        iout = 4;
      } else if ( (ndcp < DeltaV) && (dcpAlphamu < 1) ) {
        d.scale (dcp, dcpAlphamu);
        normd = ndcp*dcpAlphamu;
      } else {
        d.scale (dcp, Min(dcpAlphamu, DeltaV/ndcp) );
        normd = Min (DeltaV, ndcp*dcpAlphamu);
        iout = 13;
      }

      xc->scale (xtmp, 1);
      if (Ineq)
        sc->scale (stmp, 1);
      dx = d.get_doublex();
      for (Int i = 0; i < nvar; i++)
        xcx[i] += dx[i];
      for (Int i = 0; i < nconI; i++) {
        scx[i] += dx[nvar + i];
      }
      call_ccfsg_xc (dciFalse);
      normc = c->norm ();

      Ared = oldnormc*oldnormc - normc*normc;
      if (iout == 2)
        Pred = oldnormc*oldnormc;
      else {
        gtmp.sdmult (*J, 0, one, zero, d); // J*d
        Pred = gtmp.norm();
        Pred = - Pred*Pred - 2*gtmp.dot (ctmp);
      }

      CurrentTime = getTime() - StartTime;

    }

    if ( (Ared < kappa1*Pred) && (!Aavail) ) {
//      The algorithm has failed. A must be recomputed.

      xc->scale (xtmp,1);
      if (Ineq)
        sc->scale (stmp,1);
      c->scale (ctmp,1);
      normc = oldnormc;
      DeltaV = oldDelta;
      iout = 5;
    } else if (Ared >= kappa3*Pred)
      DeltaV = Max (kappa4*normd, DeltaV);

    return iout;

  }

}
