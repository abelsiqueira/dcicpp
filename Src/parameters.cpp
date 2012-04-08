#include "interface.h"
//#include <cassert>

namespace DCI {

  void Interface::DefineParameters () {
    DeltaMax = 1e6;
    maxrest = 200;
    maxit = 200;
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

  void Interface::Initialization () {
    ufn = 0;
    cfn = 0;
    uofg = 0;
    cofg = 0;
    uprod = 0;
    cprod = 0;
    ccfsg = 0;
    ccifg = 0;
    unames = 0;
    cnames = 0;
    f = 0;
    fxc = 0;
    g = 0;
    H = 0;
    Htrip = 0;
    c = 0;
    J = 0;
    Jtrip = 0;
    LJ = 0;
    gp = 0;
    normgp = 0;
    normg = 0;
    normc = 0;
    x = 0;
    solx = 0;
    bl = 0;
    bu = 0;
    y = 0;
    yineq = 0;
    cl = 0;
    cu = 0;
    s = 0;
    sols = 0;
    xc = 0;
    sc = 0;
    feasOpt = 0;
    equatn = 0;
    linear = 0;
    nmax = 0;
    mmax = 0;
    amax = 0;
    nvar = 0;
    ncon = 0;
    nconE = 0;
    nconI = 0;
    ineqIdx = 0;
    CurrentTime = 0;
    MaxTime = dciInf;
    DLH = dciInf;
    DLV = 0;
    Lref = dciInf;
    StartAtOne = dciFalse;
    Initialized = dciFalse;
    Running = dciFalse;
    Solved = dciFalse;
    Unlimited = dciFalse;
    Solved = 0;
    VertFlag = 0;
    iter = 0;
    itssmll = 0;
    tRest = 0;
    tSteih = 0;
    tRej = 0;
    tSoc = 0;
    tbfgs = 0;
    cholFacs = 0;
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

    //Program properties
    Ineq = dciFalse; //Has some inequalities
    Linear = dciFalse; //Only linear constraints
    Bounded = dciFalse; //Has bounds
  }

  void Interface::InitialValues () {
    Real one[2] = {1,0};

    maxitSteih = Max (minitSteih, Min(maxitSteih, Int(nvar*relitSteih+5) ) );
    minstep *= Min (csig, csic);
    
    // Calculating the function value and c without the barrier.
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
      Real smldelta = Min ( 1.0, (bui - bli)/100.0);
      xx[i] = Max ( Min ( xx[i], bui - smldelta ), bli + smldelta );
      xcx[i] = xx[i];
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
        scx[i] = sx[i];
      }

      Running = dciTrue;
      // Now, adding s.
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
    } else { //No constraints, may have bounds on the variables
      Running = dciTrue; //If there are bounds, this will get them.
      call_ofg ();
      normc = 0;
      Ln = *f;
      *gp = *g;
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
