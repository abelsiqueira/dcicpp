#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  void Interface::update_lambda () {
    if (ncon == 0) {
      gp->scale (*g, 1);
      return;
    }
    double one[2] = {1,0}, zero[2] = {0, 0};
    Vector ytmp (*env);
    NAproj (*g, *gp, ytmp);
    y->scale (ytmp, 1);

    bool ForceSign = dciTrue;
    //Force sign
    if (!ForceSign)
      return;
    double coef = 100, expo = 0.5;
    double limiter = coef * pow (mu, expo);
    for (int i = 0; i < nconI; i++) {
      int j = ineqIdx[i];
      if (clx[j] <= -dciInf)
        yx[j] = Max (yx[j], -limiter);
      else if (cux[j] >= dciInf)
        yx[j] = Min (yx[j], limiter);
      else {
        double middle = (clx[j] + cux[j])/2;
        if (scx[j] > middle)
          yx[j] = Max (yx[j], -limiter);
        else
          yx[j] = Min (yx[j], limiter);
      }
    }
  }

  void Interface::LstSqrCG (Bool transp, const Vector & rhs, Vector & sol, Vector & res) {
    /* Solves 
     * min 0.5 * || J*sol - rhs ||^2
     * or
     * min 0.5 * || J'*sol - rhs ||^2
     */
    double * solx = sol.get_doublex();
    int dimen = sol.size ();
    for (int i = 0; i < dimen; i++)
      solx[i] = 0;
    double one[2] = {1,0}, zero[2] = {0, 0};
    Vector Jrhs(*env);
    Jrhs.sdmult (*J, !transp, one, zero, rhs);
    Vector r(Jrhs), p(*env), Jtp (*env), JJtp (*env);
    p.scale(r, -1);
    Jtp.sdmult (*J, transp, one, zero, p);
    JJtp.sdmult (*J, !transp, one, zero, Jtp);
    double rr, rrp, alpha, beta;
    rrp = r.dot(r);
    int cgiter = 1;
    while ( (rrp > 1e-6) || (cgiter < dimen) ) {
      rr = rrp;
      alpha = rr/Jtp.dot(Jtp);
      sol.saxpy(p, alpha);
      r.saxpy(JJtp, alpha);
      rrp = r.dot(r);
      beta = rrp/rr;
      p.scale(beta);
      p.saxpy(r, -1);
      Jtp.sdmult (*J, transp, one, zero, p);
      JJtp.sdmult (*J, !transp, one, zero, Jtp);
      cgiter++;
    }

    res.scale (rhs, -1);
    res.sdmult (*J, transp, one, one, sol);
  }

  void Interface::LinSysCG (Bool transp, const Vector & rhs, Vector & sol) {
    /* Solves
     * A*A' * sol = rhs
     * if transp == 0, or
     * A'*A * sol = rhs
     * if transp == 1.
     */
    sol.scale (rhs, 0);
    int dimen = sol.size ();
    double one[2] = {1, 0}, zero[2] = {0, 0};

    Vector r(*env), p(*env), Jtp(*env), JJtp (*env);
    r.scale (rhs, -1);
    p.scale (rhs, 1);
    Jtp.sdmult (*J, !transp, one, zero, p);
    JJtp.sdmult (*J, transp, one, zero, Jtp);
    double rr, rrp, alpha, beta;
    rrp = r.dot(r);
    int cgiter = 1;
    while ( (rrp > 1e-6) && (cgiter < dimen) ) {
      rr = rrp;
      alpha = rr/Jtp.dot(Jtp);
      sol.saxpy (p, alpha);
      r.saxpy (JJtp, alpha);
      rrp = r.dot(r);
      beta = rrp/rr;
      p.scale (beta);
      p.saxpy (r, -1);
      Jtp.sdmult (*J, !transp, one, zero, p);
      JJtp.sdmult (*J, transp, one, zero, Jtp);
      cgiter++;
    }
  }

  void Interface::update_mu () {
//    Real minmuthisiter = Max (0.1*mu, 1e-24);
    Real minother = 100*Min (rho, rho*rho); 
    minother = Min (minother, gap + lagrgap + infacgap);
    minother = Min (minother, mu);
    minother = Min (minother, 100*normck);

//    mu = Max (minother, minmuthisiter);
    mu = minother;
  }

  void Interface::NAproj (Vector & r, Vector & Pr, Vector & tmp) {
    Real mone[2] = {-1,0}, zero[2] = {0,0};
    Real one[2] = {1,0};

    Pr.sdmult (*J, 0, mone, zero, r); // Pr = -A*r
    tmp.solve (CHOLMOD_A, *LJ, Pr); // A * A' * tmp = Pr
    if (LimLbd) {
      pReal tmpx = tmp.get_doublex();
      for (Int i = 0; i < ncon; i++) {
        Real tmpi = tmpx[i];
        if (tmpi > LbdMax) tmpx[i] = LbdMax;
        else if (tmpi < -LbdMax) tmpx[i] = -LbdMax;
      }
    }
    Pr.scale (r, 1); //Pr = r
    Pr.sdmult (*J, 1, one, one, tmp); //Pr = Pr + J'*tmp => Pr = r - A'*inv(A*A')*A*r

  }

  Int Interface::NAstep (Vector & r, Vector & dr) { //dr = -A'*inv(AA')*r
    Real mone[2] = {-1,0}, zero[2] = {0,0};
    if (!UseCG) {
      dr.solve (CHOLMOD_A, *LJ, r); // dr = inv(AA')r;
    } else {
      LinSysCG (0, r, dr);
    }
    dr.sdmult (*J, 1, mone, zero, dr);

    return 0;
  }

  Int Interface::quadstep (Vector & s, Real & qs) { //Not really used
    if (s.dense == 0) {
      std::cerr << "Allocate s first" << std::endl;
      return -1;
    }

    Real tol = 1e-6;
    Real rtw = 0, dtBd = 0;
    Real alpha = 0;

    Vector r (*g), w (*env), tmp (*env);
    Vector d(*env);
    Vector Bd(*env);

    pReal sx = s.get_doublex();
    for (Int i = 0; i < nvar; i++)
      sx[i] = 0;

    NAproj (r, w, *y);

    d.scale (w, -1);
    Bd.scale (w, -1);
    pReal dx = d.get_doublex();
    pReal Bdx = Bd.get_doublex();

    call_prod (dciFalse, dx, Bdx);
    rtw = r.dot(w);
    dtBd = d.dot(Bd);

    while (rtw >= tol) {
      if (dtBd < tol) {
        // Negative curvature
        Real normd_sqr = d.dot(d);
        Real norms_sqr = s.dot(s);
        Real dts = s.dot(d);
        alpha = sqrt (dts*dts - normd_sqr * (norms_sqr - DeltaH*DeltaH) );
        Real dtr = d.dot(r);
        Real tau = (-dts + alpha)/normd_sqr;
        call_prod (dciTrue, sx, Bdx);
        qs = 0.5*s.dot(Bd) + g->dot(s);
        Real qplus =  qs + tau*dtr + 0.5*tau*tau*dtBd;
        tau = (-dts - alpha)/normd_sqr;
        Real qminus = qs + tau*dtr + 0.5*tau*tau*dtBd;
        if (qplus < qminus) {
          tau = (-dts + alpha)/normd_sqr;
          qminus = qplus;
        }
        qs = qminus;
        s.saxpy (d, tau);
        assert (s.norm () <= DeltaH + 1e-6);
        return 3;
      }
      alpha = rtw/dtBd;
      Vector splus(s);
      splus.saxpy (d, alpha);
      if (splus.norm() > DeltaH) {
        Real normd_sqr = d.dot(d);
        Real norms_sqr = s.dot(s);
        Real dts = s.dot(d);
        alpha = -dts + sqrt (dts*dts - normd_sqr * (norms_sqr - DeltaH*DeltaH) );
        alpha /= normd_sqr;
        s.saxpy (d, alpha);
        call_prod (dciTrue, sx, Bdx);
        qs = 0.5*s.dot(Bd) + g->dot(s);
        assert (s.norm () <= DeltaH + 1e-6);
        return 1;
      }
      s.scale (splus, 1);
      r.saxpy (Bd, alpha);
      NAproj (r, w, tmp);
      Real rtwplus = r.dot(w);
      d.scale(d, rtwplus/rtw);
      d.saxpy(w, -1);
      rtw = rtwplus;
    }

    call_prod (dciTrue, sx, Bdx);
    qs = 0.5*s.dot(Bd) + g->dot(s);
    assert (s.norm () <= DeltaH);

    return 0;

  }

  void Interface::scale_x (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xx[i], Li = blx[i], Ui = bux[i];
      if (Li > -dciInf)
        Vx[i] *= (zi - Li);
      if (Ui < dciInf)
        Vx[i] *= (Ui - zi);
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = sx[i], Li = clx[ineqIdx[i]], Ui = cux[ineqIdx[i]];
      Int j = nvar + i;
      if (Li > -dciInf)
        Vx[j] *= (zi - Li);
      if (Ui < dciInf)
        Vx[j] *= (Ui - zi);
    }
  }

  void Interface::scale_xc (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], Li = blx[i], Ui = bux[i];
      if (Li > -dciInf)
        Vx[i] *= (zi - Li);
      if (Ui < dciInf)
        Vx[i] *= (Ui - zi);
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = scx[i], Li = clx[ineqIdx[i]], Ui = cux[ineqIdx[i]];
      Int j = nvar + i;
      if (Li > -dciInf)
        Vx[j] *= (zi - Li);
      if (Ui < dciInf)
        Vx[j] *= (Ui - zi);
    }
  }

  void Interface::updyineq () {
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i], xi = xx[i];
      yineqx[i] = 0;
      if (bli > -dciInf)
        yineqx[i] -= mu/(xi - bli);
      if (bui < dciInf)
        yineqx[i] -= mu/(xi - bui);
    }
    for (Int i = 0; i < nconI; i++) {
      Real cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], si = sx[i];
      Int j = nvar + i;
      yineqx[j] = 0;
      if (cli > -dciInf)
        yineqx[j] -= mu/(si - cli);
      if (cui < dciInf)
        yineqx[j] -= mu/(si - cui);
    }
  }

  void Interface::updyineq_xc () {
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i], xi = xcx[i];
      yineqx[i] = 0;
      if (bli > -dciInf)
        yineqx[i] -= mu/(xi - bli);
      if (bui < dciInf)
        yineqx[i] -= mu/(xi - bui);
    }
    for (Int i = 0; i < nconI; i++) {
      Real cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], si = scx[i];
      Int j = nvar + i;
      yineqx[j] = 0;
      if (cli > -dciInf)
        yineqx[j] -= mu/(si - cli);
      if (cui < dciInf)
        yineqx[j] -= mu/(si - cui);
    }
  }

  Real Interface::calc_gap () {
    Real gap = 0;
    for (Int i = 0; i < nvar; i++)  {
      Real xi = xx[i], bli = blx[i], bui = bux[i];
      if ( (bli <= -dciInf) && (bui >= dciInf) )
        continue;
      if (bli > -dciInf)
        gap += Max(0.0, -yineqx[i]) * (xi - bli);
      if (bui < dciInf)
        gap += Max(0.0,  yineqx[i]) * (bui - xi);
    }
    for (Int i = 0; i < nconI; i++) {
      Real si = sx[i], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
      Int j = nvar + i;
      if (cli > -dciInf)
        gap += Max(0.0, -yineqx[j]) * (si - cli);
      if (cui < dciInf)
        gap += Max(0.0,  yineqx[j]) * (cui - si);
    }
    gap /= (nvar + nconI);
    
    return gap;
  }

  Real Interface::calc_ydif () {
    Real dif = 0;
    for (Int i = 0; i < nconI; i++)
      dif += fabs (yx[ineqIdx[i]] - yineqx[nvar + i]);
    return dif;
  }

  Real Interface::calc_pen () {
    Real val = 0.0;
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > -dciInf)
        val += log (xx[i] - blx[i]);
      if (bux[i] < dciInf)
        val += log (bux[i] - xx[i]);
    }
    for (Int i = 0; i < nconI; i++) {
      if (clx[ineqIdx[i]] > -dciInf)
        val += log (sx[i] - clx[ineqIdx[i]]);
      if (cux[ineqIdx[i]] < dciInf)
        val += log (cux[ineqIdx[i]] - sx[i]);
    }
    return val;
  }

  Real Interface::calc_pen_xc () {
    Real val = 0.0;
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] > -dciInf)
        val += log (xcx[i] - blx[i]);
      if (bux[i] < dciInf)
        val += log (bux[i] - xcx[i]);
    }
    for (Int i = 0; i < nconI; i++) {
      if (clx[ineqIdx[i]] > -dciInf)
        val += log (scx[i] - clx[ineqIdx[i]]);
      if (cux[ineqIdx[i]] < dciInf)
        val += log (cux[ineqIdx[i]] - scx[i]);
    }
    return val;
  }

  Real Interface::penvthv (const Vector & V) const {
    pReal px = V.get_doublex ();
    Real val = 0.0;

    for (Int i = 0; i < nvar + nconI; i++)
      val -= px[i]*px[i]*yineqx[i]/mu;

    return val;
  }

}
