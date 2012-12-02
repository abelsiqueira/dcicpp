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

    //Remove later if needed
    call_ccfsg (dciTrue, ScaleVertical);
    Aavail = dciFalse;

    gtmp.sdmult (*J, 1, one, zero, ctmp); // g = J'*c
//    scale_xc(gtmp);
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
      lower[i] = Max((li - zi) * (1 - epsmu), -DeltaV);
      upper[i] = Min((ui - zi) * (1 - epsmu),  DeltaV);
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real zi = scx[i], li = clx[ineqIdx[i]], ui = cux[ineqIdx[i]];
      lower[j] = Max((li - zi) * (1 - epsmu), -DeltaV);
      upper[j] = Min((ui - zi) * (1 - epsmu),  DeltaV);
    }
    d.sdmult(*J, 0, one, zero, gtmp);
    alpha = normgtmp/d.norm();
    alpha *= alpha;
    if (!ScaleVertical) {
      for (Int i = 0; i < nvar + nconI; i++) {
        Real gtmpi = gtmpx[i];
        if (gtmpi > 0)
          alpha = Min (alpha, upper[i]/gtmpi);
        else if (gtmpi < 0)
          alpha = Min (alpha, lower[i]/gtmpi);
        assert (alpha > 0);
      } 
    } else {
      for (Int i = 0; i < nvar + nconI; i++) {
        Real gtmpi = gtmpx[i];
        if (gtmpi > 0)
          alpha = Min (alpha, upper[i]/(Lambda[i]*gtmpi));
        else if (gtmpi < 0)
          alpha = Min (alpha, lower[i]/(Lambda[i]*gtmpi));
        assert (alpha > 0);
      } 
    }
    dcp.scale (gtmp, -alpha);
    if (ScaleVertical)
      scale_xc (dcp);
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
    pReal dcpx = dcp.get_doublex();
    pReal dx = 0;
    pReal dnx = 0;
//    Real smlAlphamu = 1e-3;
    bool dcpOutsideRegion = false;
    for (Int i = 0; i < nvar + nconI; i++) {
      if ( (dcpx[i] >= upper[i]) || (dcpx[i] <= lower[i]) ) {
        dcpOutsideRegion = true;
        break;
      }
    }

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
    ndcp = dcp.norm(0);

    while ( (Ared < kappa1*Pred) && (Aavail || (TrustIter < 2) ) && (CurrentTime < MaxTime) ) {
      TrustIter++;
      DeltaV = kappa2*normd;
//      DeltaV = Min(kappa2*normd, 0.9*DeltaV);

      if (!dcpOutsideRegion) {
        // Cauchy step is inside trust region
        // sc + dcps >= epsmu*sc
        if (!dnavail) {
          naflag = NAstep (ctmp, dn); 
          dnavail = dciTrue;
          if (ScaleVertical)
            scale_xc (dn);

          if (naflag > 1)
            ndn = 0;
          else
            ndn = dn.norm (0);
          
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

        ndn = dn.norm(0);
        bool dnOutsideRegion = false;
        for (Int i = 0; i < nvar + nconI; i++) {
          if ( (dnx[i] > upper[i]) || (dnx[i] < lower[i]) ) {
            dnOutsideRegion = true;
            break;
          }
        }
        if (ndn <= ndcp) {
          // dn is too small. Use dcp.
          d = dcp;
          normd = ndcp;
          iout = 1;
        } else if (!dnOutsideRegion) {
          d = dn;
          normd = ndn;
          iout = 2;
        } else {
          // dn outside region and dcp inside
          Real convAux = 1.0;
          for (Int i = 0; i < nvar; i++) {
            Real difx = dnx[i] - dcpx[i];
            if (difx > 0)
              convAux = Min(convAux, (upper[i] - dcpx[i])/difx);
            else if (difx < 0)
              convAux = Min(convAux, (lower[i] - dcpx[i])/difx);
          }
          for (Int j = 0; j < nconI; j++) {
            Int i = nvar + j;
            Real difx = dnx[i] - dcpx[i];
            if (difx > 0)
              convAux = Min(convAux, (upper[i] - dcpx[i])/difx);
            else if (difx < 0)
              convAux = Min(convAux, (lower[i] - dcpx[i])/difx);
          }
          d.scale (dn, convAux);
          dx = d.get_doublex();
          convAux = 1 - convAux;
          for (Int i = 0; i < nvar + nconI; i++)
            dx[i] += convAux * dcpx[i];
          normd = d.norm(0);
          iout = 10;
        }
      } else {
        Real alphadcp = 1.0;
        for (Int i = 0; i < nvar + nconI; i++) {
          Real di = dcpx[i];
          if (di > 0)
            alphadcp = Min (alphadcp, upper[i]/di);
          else if (di < 0)
            alphadcp = Min (alphadcp, lower[i]/di);
          assert (alphadcp > 0);
        } 
        d.scale(dcp, alphadcp);
        normd = ndcp;
        iout = 13;
      }

//      std::cout << "iout = " << iout << std::endl;
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
