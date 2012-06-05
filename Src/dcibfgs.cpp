#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  Int Interface::dcibfgs (const Vector & dn, Int & ibfgs) {
    Vector gtmp (*env);
    Vector p (*env, nvar + nconI), xold (*x), sold (*env), gold (*env, nvar + nconI), pold (*env, nvar + nconI), d (*env, nvar + nconI);
    Vector Y (*env, (nvar + nconI) * (bfgsupd + 1) );
    Vector D (*env, (nvar + nconI) * (bfgsupd + 1) );
    Vector U (*env, (nvar + nconI) * bfgsupd );
    pReal Yx = Y.get_doublex (), Dx = D.get_doublex (), Ux = U.get_doublex (),
          gtmpx = 0, goldx = gold.get_doublex (), 
          px = p.get_doublex (), poldx = pold.get_doublex (),
          xoldx = xold.get_doublex (), soldx = 0;
    Real one[2] = {1,0}, zero[2] = {0,0};
    Real angle, objfun, gnorm, gtd, pnorm, znorm;
    Real dnorm = dn.norm ();
    Real alpha[bfgsupd], csi[bfgsupd];
    Int iout;

    if (Ineq) {
      sold = *s;
      soldx = sold.get_doublex();
    }

    checkInfactibility ();

    objfun = 0.5 * normc * normc;
    gtmp.sdmult (*J, 1, one, zero, *c);
    gtmpx = gtmp.get_doublex ();
    if (penal_bfgs) {
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
    gnorm = gtmp.norm ();
    ibfgs = 0;
    iout = 0;

    checkInfactibility ();

    if (normc < rho)
      return -1;

    if (gnorm < 1e-10) {
      iout = 3;
      return iout;
    }

    p.scale (gtmp, -1);
    pnorm = gnorm;
    // davail = dciTrue always
    if (dnorm < 1e-6)
      angle = 0;
    else
      angle = p.dot (dn) / (dnorm * pnorm);
    if ( (angle >= 0) && (angle < 1e-4) )
      d = p;
    else
      d = dn;

    gtd = d.dot (gtmp);
    xold = *xc;
    if (Ineq)
      sold = *sc;

    gold = gtmp;

    iout = linesearch (xold, sold, d, objfun, gtd, gtmp);

    gtmp.sdmult (*J, 1, one, zero, *c);
    
    dnorm = 0;
    for (Int i = 0; i < nvar; i++) {
      Real Dxi;
      Yx[i] = gtmpx[i] - goldx[i];
      Dxi = xcx[i] - xoldx[i];
      Dx[i] = Dxi;
      dnorm += Dxi*Dxi;
    }
    if (Ineq) {
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        Real Dxi;
        Yx[j] = gtmpx[j] - goldx[j];
        Dxi = scx[i] - soldx[i];
        Dx[j] = Dxi;
        dnorm += Dxi*Dxi;
      }
    }
    dnorm = sqrt(dnorm);
    normc = c->norm ();
    if (Ineq) {
      Real xnorm = xc->norm (), snorm = sc->norm ();
      znorm = sqrt(xnorm*xnorm + snorm*snorm);
    } else
      znorm = xc->norm ();

    CurrentTime = getTime() - StartTime;
    while ( (normc > rho) && (ibfgs < bfgsupd) && (iout == 0) && (dnorm >= znorm*1e-17) && (CurrentTime < MaxTime) ) {
      alpha [ibfgs] = 0;
      for (Int i = 0; i < nvar + nconI; i++)
        alpha [ibfgs] += Yx[ibfgs*(nvar + nconI) + i] * Dx[ibfgs*(nvar + nconI) + i];
      if (alpha[ibfgs] < 1e-12)
        alpha[ibfgs] = 1e-12;

      pold = p;
      p = gtmp;

      for (Int i = 0; i < ibfgs; i++)
        Hiprod (i, alpha[i], csi[i], Dx, Ux, Yx, gtmpx, px);

      for (Int i = 0; i < nvar + nconI; i++)
        Ux[ibfgs*(nvar + nconI) + i] = px[i] + poldx[i];

      csi[ibfgs] = 0;
      for (Int i = 0; i < nvar + nconI; i++)
        csi[ibfgs] += Yx[ibfgs*(nvar + nconI) + i] * Ux[ibfgs*(nvar + nconI) + i];
      csi[ibfgs] /= alpha[ibfgs];
      csi[ibfgs] += 1;

      Hiprod (ibfgs, alpha[ibfgs], csi[ibfgs], Dx, Ux, Yx, gtmpx, px);
      p.scale (-1);

      gtd = gtmp.dot(p);

      xold = *xc;
      if (Ineq)
        sold = *sc;
      gold = gtmp;

      iout = linesearch (xold, sold, p, objfun, gtd, gtmp);

      checkInfactibility();
      
      gtmp.sdmult (*J, 1, one, zero, *c);

      ibfgs++;

      dnorm = 0;
      for (Int i = 0; i < nvar; i++) {
        Yx[ibfgs*(nvar + nconI) + i] = gtmpx[i] - goldx[i];
        Real Dxi = xcx[i] - xoldx[i];
        Dx[ibfgs*(nvar + nconI) + i] = Dxi;
        dnorm += Dxi*Dxi;
      }
      if (Ineq) {
        for (Int i = 0; i < nconI; i++) {
          Int j = nvar + i;
          Yx[ibfgs*(nvar + nconI) + j] = gtmpx[j] - goldx[j];
          Real Dxi = scx[i] - soldx[i];
          Dx[ibfgs*(nvar + nconI) + j] = Dxi;
          dnorm += Dxi*Dxi;
        }
      }
      dnorm = sqrt (dnorm);
      normc = c->norm ();
      if (Ineq) {
        Real xnorm = xc->norm (), snorm = sc->norm ();
        znorm = sqrt(xnorm*xnorm + snorm*snorm);
      } else
        znorm = xc->norm ();

      CurrentTime = getTime() - StartTime;

    }

    ibfgs++;
    if (normc <= rho)
      iout = 0;
    else if (ibfgs > bfgsupd)
      iout = 2;
    else if (dnorm < znorm*1e-10)
      iout = 4;

    return iout;

  }

  void Interface::Hiprod (Int i, Real alpha, Real csi, pReal d, pReal u, pReal y,
      pReal v, pReal Hv) {

    Real cte1 = 0, cte2 = 0;

    for (Int j = 0; j < nvar; j++) {
      cte1 += d[i*nvar + j] * v[j];
      cte2  += y[i*nvar + j] * Hv[j];
    }
    cte1 /= alpha;
    cte2 /= alpha;

    for (Int j = 0; j < nvar; j++) {
      Hv[j] += d[i*nvar + j] * (csi*cte1 - cte2) - cte1 * u[i*nvar + j];
    }

  }

  Int Interface::linesearch (const Vector & x0, const Vector & s0, const Vector & d, Real & objfun, Real & gtd, Vector & gtmp) {
    Real f0 = objfun;
    Real fold = f0;
    Real one[2] = {1,0}, zero[2] = {0,0};
    Real gtd0 = gtd;
    Bool ishi;
    Real lbdold = 0;
    Real lambda = 1;
    pReal dx = d.get_doublex();
    Real maxStepS = dciInf;
    Real lbdmax = dciInf;
    Real smlStepS = 1e-3;
    Bool bfgsfirst = dciTrue;

    *xc = x0;
    if (Ineq)
      *sc = s0;

    Vector Diag(*env);
    Diag.reset (nvar + nconI, 1);
    pReal Diagx = Diag.get_doublex();
    if (ScaleVertical)
      scale_xc (Diag);

    for (Int i = 0; i < nvar; i++) {
      Real xi = xcx[i], di = Diagx[i]*dx[i], bli = blx[i], bui = bux[i];
      if (fabs(di) < dciEps) {
        dx[i] = 0.0;
        continue;
      }
      if (di == 0)
        continue;
      if (di < 0) {
        Real val = (bli - xi)*(1 - epsmu)/di;
        if (project_bfgs) {
          if (val < 1)
            dx[i] *= val/Diagx[i];
        } else {
          maxStepS = Min (maxStepS, val);
        }
      } else {
        Real val = (bui - xi)*(1 - epsmu)/di;
        if (project_bfgs) {
          if (val < 1)
            dx[i] *= val/Diagx[i];
        } else {
          maxStepS = Min (maxStepS, val);
        }
      }
    }

    for (Int i = 0; i < nconI; i++) {
      Real si = scx[i], di = Diagx[nvar + i] * dx[nvar + i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
      if (fabs(di) < dciEps) {
        dx[nvar + i] = 0.0;
        continue;
      }
      if (di == 0)
        continue;
      if (di < 0) {
        Real val = (cli - si)*(1 - epsmu)/di;
        if (project_bfgs) {
          if (val < 1)
            dx[nvar + i] *= val/Diagx[nvar + i];
        } else {
          maxStepS = Min (maxStepS, val);
        }
      } else {
        Real val = (cui - si)*(1 - epsmu)/di;
        if (project_bfgs) {
          if (val < 1)
            dx[nvar + i] *= val/Diagx[nvar + i];
        } else {
          maxStepS = Min (maxStepS, val);
        }
      }
    }

    lbdmax = Min (lbdmax, maxStepS);
    if (maxStepS < 1)
      lambda = lbdmax/2;

    while (dciTrue) {
      *xc = x0;
      if (Ineq)
        *sc = s0;

      for (Int i = 0; i < nvar; i++) {
        xcx[i] += lambda*Diagx[i]*dx[i];
        if (xcx[i] == bux[i])
          xcx[i] = bux[i] - dciEps;
        else if (xcx[i] == blx[i])
          xcx[i] = blx[i] + dciEps;
      }
      for (Int i = 0; i < nconI; i++) {
        scx[i] += lambda*Diagx[nvar + i]*dx[nvar + i];
        if (scx[i] == cux[i])
          scx[i] = cux[i] - dciEps;
        else if (scx[i] == clx[i])
          scx[i] = clx[i] + dciEps;
      }

      
      checkInfactibility ();

      call_ccfsg_xc (dciFalse);
      normc = c->norm ();
      objfun = 0.5 * normc * normc;

      if ( (objfun > f0 + c1*lambda*gtd0) || ( (objfun >= fold) && (!bfgsfirst) ) ) {
        ishi = dciTrue;
        zoom (x0, s0, d, gtmp, f0, lbdold, lambda, objfun, gtd, objfun, fold, ishi, lambda);
        return 2;
      }

      gtmp.sdmult (*J, 1, one, zero, *c);
      gtd = gtmp.dot(d);

      if (fabs(gtd) <= -c2*gtd0)
        return 3;

      if (gtd >= 0.0) {
        ishi = dciFalse;
        zoom (x0, s0, d, gtmp, f0, lbdold, lambda, objfun, gtd, objfun, fold, ishi, lambda);
        return 4;
      }

      if (lambda == lbdmax)
        return 5;

      lbdold = lambda;
      fold = objfun;
      lambda = Min (2*lambda, lbdmax);
      bfgsfirst = dciFalse;
    }

    return 0;

  }

  Int Interface::zoom (const Vector & x0, const Vector & s0, const Vector & d, Vector & gtmp, Real f0, Real lbdlo0, Real lbdhi0, Real & objfun, Real & gtd, Real flo0, Real fhi0, Real ishi0, Real lambda) {
    Bool scaleJ = dciTrue;

    if (lbdlo0 == lbdhi0) {
      lambda = lbdlo0;
      return 0;
    }

    Int itermax = 10;
    Int zoomiter = 1;
    Bool enough = dciFalse;
    Real flo = flo0, fhi = fhi0, lbdlo = lbdlo0, lbdhi = lbdhi0, ishi = ishi0;
    pReal x0x = x0.get_doublex(), s0x = 0, dx = d.get_doublex ();
    Real one[2] = {1,0}, zero[2] = {0,0};
    Real gtd0 = gtd;
    Bool zoomfirst;
    if (Ineq)
      s0x = s0.get_doublex();


    Vector Diag(*env);
    Diag.reset (nvar + nconI, 1);
    pReal Diagx = Diag.get_doublex();
    if (ScaleVertical)
      scale_xc (Diag);

    while ( (zoomiter <= itermax) && (!enough) ) {

      zoomfirst = ( (zoomiter == 1) || (lbdlo == 0.0) || (lbdhi == 0.0) );
      if (ishi)
        lambda = interpolate (f0, gtd0, fhi, flo, lbdhi, lbdlo, zoomfirst);
      else
        lambda = interpolate (f0, gtd0, flo, fhi, lbdlo, lbdhi, zoomfirst);

      Real diff = fabs(lbdlo - lbdhi);
      lambda = Max (Min (lbdlo, lbdhi) + 0.01*diff, 
               Min (lambda, Max(lbdlo, lbdhi) - 0.01*diff) );

      for (Int i = 0; i < nvar; i++)
        xcx[i] = x0x[i] + lambda * Diagx[i]*dx[i];
      for (Int i = 0; i < nconI; i++)
        scx[i] = s0x[i] + lambda * Diagx[nvar + i]*dx[nvar + i];
      call_ccfsg_xc (dciFalse);
      normc = c->norm ();
      objfun = 0.5 * normc * normc;

      checkInfactibility ();

      if ( (objfun > f0 + c1*lambda*gtd0) || (objfun >= flo) ) {
        lbdhi = lambda;
        fhi = objfun;
        ishi = dciTrue;
        zoomiter = zoomiter + 1;
      } else {
        call_ccfsg_xc (dciTrue, ScaleVertical);
        gtmp.sdmult (*J, 1, one, zero, *c);
        gtd = gtmp.dot(d);

        if (fabs(gtd) <= -c2*gtd0)
          enough = dciTrue;
        else {
          if (gtd*(lbdhi - lbdlo) >= 0.0) {
            lbdhi = lbdlo;
            fhi = flo;
          }

          lbdlo = lambda;
          flo = objfun;
          ishi = dciFalse;
          zoomiter++;
        }
      }
    }

    if (zoomiter > itermax) {

      if (ishi) {
        call_ccfsg_xc (dciTrue, ScaleVertical);
        gtmp.sdmult (*J, 1, one, zero, *c);
        gtd = gtmp.dot(d);
      }
      return 1;

    }

    return 0;

  }

  Real Interface::interpolate (Real f0, Real fl0, Real fm1, Real fm2, Real lm1, Real lm2, Bool zoomfirst) {

    Real lambda;

    if (zoomfirst) {
      // Quadratic function
      lambda = -0.5*fl0*lm1*lm1/(fm1 - f0 - fl0 * lm1);
    } else {
      // Cubic function
      Real cte1 = ( -(lm2/(lm1*lm1))*(fm1-f0-fl0*lm1)+
          (lm1/(lm2*lm2))*(fm2-f0-fl0*lm2))/(lm1-lm2);
      Real cte2 = ( (fm1-f0-fl0*lm1)/(lm1*lm1) - (fm2-f0-fl0*lm2)/(lm2*lm2) )/(lm1 - lm2);
      Real delta = cte1*cte1 - 3*fl0*cte2;

      if (delta >= 0.0)
        lambda = (-cte1 + sqrt(delta))/(3*cte2);
      else
        lambda = -0.5*fl0*lm1*lm1/(fm1-f0-fl0*lm1);
    }

    if ( (fabs(lambda - lm1) < lm1*1e-3) || (lambda < lm1*1e-3) )
      lambda = 0.5*lm1;

    return lambda;

  }

}

