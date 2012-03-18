#include "interface.h"
//#include <cassert>

namespace DCI {

  void Interface::InitialParameters () {
    DeltaMax = 1e6;
    maxrest = 20000;
    maxit = 200000;
    maxssmll = 5;
    maxitSteih = 100000;
    minitSteih = 100;
    relitSteih = 10;
    nfailv = 3;
    csic = 1e-6;
    csig = 1e-6;
    rhomin = 1e-8;
    phi1 = 1;
    phi2 = 0.95;
    kappa1 = 1e-4;
    kappa2 = 0.25;
    kappa3 = 0.7;
    kappa4 = 2.5;
    zeta1 = 2;
    zeta2 = 1;
    zeta3 = 5;
    alphaR = 0.25;
    alphaI = 2.5;
    alphaS = 6.25e-2;
    eta1 = 1e-4;
    eta2 = 0.7;
    eta3 = 0.1;
    DeltaMin = 1e-4;
    DeltaInf = 1e-10;
    minstep = 1e-3;
    Delta0 = 1e5;
    thetaR = 0.9;
    LbdMax = 1e6;
    eps1 = 1e-6;
    eps2 = 1e-14;
    eps3 = 1e-8;
    epsmu = 1e-6;
    epsgap = 1e-4;
    bfgsupd = 5;
    c1 = 0.5;
    c2 = 5e-1;
    MaxTime = 5 * 60; // 5 minutes
    minBk = 1e-12;
  }

  void Interface::InitialValues () {
    Real one[2] = {1,0};

    Solved = 0;
    VertFlag = 0;
    iter = 0;
    itssmll = 0;
    tRest = 0;
    tSteih = 0;
    tRej = 0;
    tSoc = 0;
    tbfgs = 0;
    nSoc = 0;
    nbfgs = 0;
    nRej = 0;
    nSteih = 0;
    nHprod = 0;
    nRest = 0;
    Lref = dciInf;
    DLH = dciInf;
    DLV = 0;
    normck = dciInf;
    mu = 1;
    murho = 1;
    mugap = 1;

    maxitSteih = Max (minitSteih, Min(maxitSteih, Int(nvar*relitSteih+5) ) );
    minstep *= Min (csig, csic);
    
    call_fn ();
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i];
      if ( (xx[i] > bli) && (xx[i] < bui) ) {
        xcx[i] = xx[i];
        continue;
      }
      // Fixed variables WORKAROUND:
      if (bli == bui) {
        blx[i] -= 1e-12;
        bux[i] += 1e-12;
        bli = blx[i];
        bui = bux[i];
      }
      Real smldelta = Min ( 1.0, (bui - bli)/100);
      xx[i] = Max ( Min ( xx[i], bui - smldelta ), bli + smldelta );
      xcx[i] = Max ( Min ( xx[i], bui - smldelta ), bli + smldelta );
      assert (xx[i] > bli);
      assert (xx[i] < bui);
    }

    if (ncon > 0) {
      for (Int i = 0; i < nconI; i++) {
        Real cxi = cx[ineqIdx[i]], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        if ( (sx[i] > cli) && (sx[i] < cui) ) {
          scx[i] = sx[i];
          continue;
        }
        Real smldelta = Min ( 1.0, (cui - cli)/100);
        assert (smldelta > 0);
        sx[i] = Max ( Min ( cxi, cui - smldelta ), cli + smldelta );
        scx[i] = Max ( Min ( cxi, cui - smldelta ), cli + smldelta );
      }

      Running = dciTrue;
      call_ccfsg ();
      normc = c->norm ();
      call_ofg ();

      this->analyze_J ();
      this->cholesky_J ();

      if (y == 0) {
        y = new Vector (*env, ncon);
        yx = y->get_doublex();
        for (Int i = 0; i < ncon; i++)
          yx[i] = 0;
      }
      update_lambda ();
//      gp->scale (*g, 1);
//      gp->sdmult (*J, 1, one, one, *y);
      Ln = *f + y->dot (*c);
    } else {
      call_ofg ();
      normc = 0;
      Ln = *f;
      gp->scale (*g, 1);
    }

    normgp = gp->norm ();
    ngp = normgp / (g->norm () + 1);
    engp = ngp;
    rhomax = Max (csic, Max (5.1*normc, 50*ngp) );
    rho = Min (phi1*rhomax*ngp, 0.75*rhomax);
    if ( (rho < csic) && (normc > 100*csic) )
      rho = 0.1*normc;
    else if (ngp <= 5*csig)
      rho = csic;

    DeltaV = Max (10*x->norm (), Delta0);
    DeltaH = DeltaV;

    Running = dciFalse;
      
  }

}
