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
      stmp = *sc;


    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
    pReal gtmpx = gtmp.get_doublex();
    if (penal_trust) {
      for (Int i = 0; i < nvar; i++) {
        Real val = 0;
        if ( (bux[i] < dciInf) && (blx[i] > -dciInf) ) {
          if (PartialPenal) {
            if ( (xcx[i] - blx[i]) < (bux[i] - xcx[i]) ) {
              val = 1;
            } else {
              val = -1;
            }
          } else {
            val = bux[i] + blx[i] - 2*xcx[i];
          }
        } else if (bux[i] < dciInf) {
          val = -1;
        } else if (blx[i] > -dciInf) {
          val = 1;
        }
        gtmpx[i] -= mu*val;
      }
      for (Int i = 0; i < nconI; i++) {
        Real val = 0;
        if ( (cux[ineqIdx[i]] < dciInf) && (clx[ineqIdx[i]] > -dciInf) ) {
          if (PartialPenal) {
            if ( (scx[i] - clx[ineqIdx[i]]) < (cux[ineqIdx[i]] - scx[i]) ) {
              val = 1;
            } else {
              val = -1;
            }
          } else {
            val = cux[ineqIdx[i]] + clx[ineqIdx[i]] - 2*scx[i];
          }
        } else if (cux[ineqIdx[i]] < dciInf)
          val = -1;
        else if (clx[ineqIdx[i]] > -dciInf)
          val = 1;
        gtmpx[nvar + i] -= mu*val;
      }
    }
    normgtmp = gtmp.norm ();

    if (normgtmp < dciTiny) {
      normc = oldnormc;
      iout = 6;
//      std::cout << "iout = 6" << std::endl;
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
        if (project_dcp) {
          if (val < 1)
            dcpx[i] *= val;
        } else {
          dcpAlphamu = Min (dcpAlphamu, val);
        }
      } else {
        Real val = (bui - xi)*(1 - epsmu)/dcpi;
        if (project_dcp) {
          if (val < 1)
            dcpx[i] *= val;
        } else {
          dcpAlphamu = Min (dcpAlphamu, val);
        }
      }
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real si = scx[i], dcpi = dcpx[j], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
      if (dcpi == 0)
        continue;
      if (dcpi < 0) {
        Real val = (cli - si)*(1 - epsmu)/dcpi;
        if (project_dcp) {
          if (val < 1)
            dcpx[j] *= val;
        } else {
          dcpAlphamu = Min (dcpAlphamu, val);
        }
      } else {
        Real val = (cui - si)*(1 - epsmu)/dcpi;
        if (project_dcp) {
          if (val < 1)
            dcpx[j] *= val;
        } else {
          dcpAlphamu = Min (dcpAlphamu, val);
        }
      }
    }
    ndcp = dcp.norm();

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
              if (project_dn) {
                if (val < 1)
                  dnx[i] *= val;
              } else {
                dnAlphamu = Min (dnAlphamu, val);
              }
            } else {
              Real val = (bui - xi)*(1 - epsmu)/dni;
              if (project_dn) {
                if (val < 1)
                  dnx[i] *= val;
              } else {
                dnAlphamu = Min (dnAlphamu, val);
              }
            }
          }
          for (Int i = 0; i < nconI; i++) {
            Real si = scx[i], dni = dnx[nvar + i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
            if (dni == 0)
              continue;
            if (dni < 0) {
              Real val = (cli - si)*(1 - epsmu)/dni;
              if (project_dn) {
                if (val < 1)
                  dnx[nvar + i] *= val;
              } else {
                dnAlphamu = Min (dnAlphamu, val);
              }
            } else {
              Real val = (cui - si)*(1 - epsmu)/dni;
              if (project_dn) {
                if (val < 1)
                  dnx[nvar + i] *= val;
              } else {
                dnAlphamu = Min (dnAlphamu, val);
              }
            }
          }
        }

        ndn = dn.norm();
        if (ndn <= ndcp) {
          // dn is too small. Use dcp.
          d = dcp;
          normd = ndcp;
          iout = 1;
        } else if ( (ndn <= DeltaV) && (dnAlphamu == 1) ) {
          d = dn;
          normd = ndn;
          iout = 2;
        } else if ( (ndn <= DeltaV) && (dnAlphamu < 1) ) {
          if (trustConvexBox) {
            Real convAux = 1.0;
            for (Int i = 0; i < nvar; i++) { 
              Real difx = dnx[i] - dcpx[i], zi = xcx[i],
                   lx = (blx[i] - zi)*(1 - epsmu) - dcpx[i], 
                   ux = (bux[i] - zi)*(1 - epsmu) - dcpx[i];
              if (difx > 0)
                convAux = Min (convAux, ux/difx);
              else if (difx < 0)
                convAux = Min (convAux, lx/difx);
            }
            for (Int i = 0; i < nconI; i++) { 
              Real difx = dnx[nvar + i] - dcpx[nvar + i], zi = scx[i],
                   lx = (clx[ineqIdx[i]] - zi)*(1 - epsmu) - dcpx[nvar + i], 
                   ux = (cux[ineqIdx[i]] - zi)*(1 - epsmu) - dcpx[nvar + i];
              if (difx > 0)
                convAux = Min (convAux, ux/difx);
              else if (difx < 0)
                convAux = Min (convAux, lx/difx);
            }
            d.scale (dn, convAux);
            dx = d.get_doublex();
            convAux = 1 - convAux;
            for (Int i = 0; i < nvar + nconI; i++)
              dx[i] += convAux * dcpx[i];
            normd = d.norm();
          } else {
            d.scale (dn, dnAlphamu);
            normd = ndn*dnAlphamu;
          }
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
          if (trustWorstdn) {
            dn.scale (dnAlphamu);
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
          } else {
            d = dcp;
            normd = ndcp;
          }
          iout = 11;
        }
      } else if ( (ndcp >= DeltaV) && (dcpAlphamu == 1) ) {
        d.scale (dcp, DeltaV/ndcp);
        normd = DeltaV;
        iout = 4;
      } else if ( (ndcp < DeltaV) && (dcpAlphamu < 1) ) {
        d.scale (dcp, dcpAlphamu);
        normd = ndcp*dcpAlphamu;
        iout = 14;
      } else {
        d.scale (dcp, Min(dcpAlphamu, DeltaV/ndcp) );
        normd = Min (DeltaV, ndcp*dcpAlphamu);
        iout = 13;
      }

      *xc = xtmp;
      if (Ineq)
        *sc = stmp;
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
