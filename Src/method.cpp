#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  Bool Interface::calc_feasibilityOpt () {
    if (ccifg == 0)
      return 0;
    if (feasOpt == 0)
      feasOpt = new Vector (*env, nvar);
    pReal fox = feasOpt->get_doublex ();
    Real ci = 0, cli = 0, cui = 0;
    Vector gradc (*env, nvar);
    pReal gcx = gradc.get_doublex ();
    Bool tmpTrue = dciTrue;
    
    for (Int i = 1; i <= ncon; i++) {
      (*ccifg) (&nvar, &i, xx, &ci, gcx, &tmpTrue);
      cli = clx[i-1];
      cui = cux[i-1];

      // If it's an inequality and inside limits, do nothing
      if ( (equatn[i-1] == dciFalse) && (ci >= cli) && (ci <= cui) )
        continue;

      for (Int j = 0; j < nvar; j++)
        fox[j] += ci * gcx[j];
    }

    return (feasOpt->norm () < 1e-6);
  }

  void Interface::update_lambda () {
    if (ncon == 0) {
      *gp = *g;
      return;
    }
    Vector ytmp (*env);
    NAproj (*g, *gp, ytmp);
    *y = ytmp;

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
    sol.reset (rhs.size(), 0);
    int dimen = sol.size ();
    double one[2] = {1, 0}, zero[2] = {0, 0};

    Vector r(*env), p(*env), Jtp(*env), JJtp (*env);
    r.scale (rhs, -1);
    p = rhs;
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
    Real mufactor = 1.0;
    Real minmuthisiter = Max (0.1*mu, 1e-24);
    Real minother = 100*Min (rho, rho*rho); 
    minother = Min (minother, gap + lagrgap + infacgap);
    minother = Min (minother, mufactor * mu);
    minother = Min (minother, 100*normck);

    mu = Max (minother, minmuthisiter);
    mu = minother;
  }

  void Interface::NAproj (Vector & r, Vector & Pr, Vector & tmp) {
    Real mone[2] = {-1,0}, zero[2] = {0,0};
    Real one[2] = {1,0};

    Pr.sdmult (*J, 0, mone, zero, r); // Pr = -A*r
    if (UseMUMPS) {
      Real rhs[nvar + nconI + Pr.size()];
      //rhs(1:nvar+nconI)     = 0
      //rhs(nvar+nconI+1:end) = -A*r
      for (Int i = 0; i < nvar + nconI; i++)
        rhs[i] = 0;
      for (Int i = 0; i < (Int) Pr.size(); i++)
        rhs[nvar + nconI + i] = Pr.get_doublex()[i];
//      std::cout << "rhs = " << std::endl;
//      for (size_t i = 0; i < nvar + nconI + Pr.size(); i++)
//        std::cout << rhs[i] << ' ';
//      std::cout << std::endl;
//      for (int i = 0; i < id.nz; i++) {
//        std::cout << id.irn[i] << ',' << id.jcn[i] << " = " << id.a[i] << std::endl;
//      }
      id.rhs = rhs;
      id.icntl[0] = -1; 
      id.icntl[1] = -1; 
      id.icntl[2] = -1; 
      id.icntl[3] = 0;
      id.job = 6;
      dmumps_c(&id);
      Bool FailedMumps = dciFalse;
      if (id.info[0] == -10)
        FailedMumps = dciTrue;
      while (FailedMumps) {
        std::cout << "Failed MUMPS" << std::endl;
        if (cholCorrection == 0)
          cholCorrection = 1e-12;
        else
          cholCorrection *= 100;
        if (cholCorrection > 1e6)
          break;
        call_ccfsg_xc (dciTrue, ScaleVertical);
        dmumps_c(&id);
        if (id.info[0] != -10)
          FailedMumps = dciFalse;;
      }
//      std::cout << "sol = " << std::endl;
//      for (size_t i = 0; i < nvar + nconI + Pr.size(); i++)
//        std::cout << rhs[i] << ' ';
//      std::cout << std::endl;
      //rhs(1:nvar+nconI)     = -A'*inv(A*A')*A*r
      //rhs(nvar+nconI+1:end) = inv(A*A')*A*r
      tmp.reset(Pr.size());
//      Pr.reset(nvar + nconI);
//      for (Int i = 0; i < nvar + nconI; i++)
//        Pr.get_doublex()[i] = rhs[i];
      for (Int i = 0; i < (Int) tmp.size(); i++)
        tmp.get_doublex()[i] = rhs[nvar + nconI + i];
//      Pr.print_more();
//      tmp.print_more();
//      Pr.sdmult(*J, 1, one, one, tmp); 
//      std::cout << "|rhs[1] + J'*rhs[2]| = " << Pr.norm() << std::endl;
//      assert( Pr.norm() < 1e-3);
      tmp.scale(-1);
    } else {
      tmp.solve (CHOLMOD_A, *LJ, Pr); // A * A' * tmp = Pr
    }
    if (LimLbd) {
      pReal tmpx = tmp.get_doublex();
      for (Int i = 0; i < ncon; i++) {
        Real tmpi = tmpx[i];
        if (tmpi > LbdMax) tmpx[i] = LbdMax;
        else if (tmpi < -LbdMax) tmpx[i] = -LbdMax;
      }
    }
    Pr = r;
    Pr.sdmult (*J, 1, one, one, tmp); //Pr = Pr + J'*tmp => Pr = r - A'*inv(A*A')*A*r

  }

  Int Interface::NAstep (Vector & r, Vector & dr) { //dr = -A'*inv(AA')*r
    Real mone[2] = {-1,0}, zero[2] = {0,0};
    if (UseCG) {
      LinSysCG (0, r, dr);
      dr.sdmult (*J, 1, mone, zero, dr);
    } else if (UseMUMPS) {
      Real rhs[nvar + nconI + r.size()];
      for (Int i = 0; i < nvar + nconI; i++)
        rhs[i] = 0;
      for (Int i = 0; i < (Int) r.size(); i++)
        rhs[nvar + nconI + i] = -r.get_doublex()[i];
      //rhs(1:N)     = 0
      //rhs(N+1:end) = -r
      id.rhs = rhs;
      id.icntl[0] = -1; 
      id.icntl[1] = -1; 
      id.icntl[2] = -1; 
      id.icntl[3] = 0;
      id.job = 6;
      dmumps_c(&id);
      Bool FailedMumps = dciFalse;
      if (id.info[0] == -10)
        FailedMumps = dciTrue;
      while (FailedMumps) {
        std::cout << "Failed MUMPS" << std::endl;
        if (cholCorrection == 0)
          cholCorrection = 1e-12;
        else
          cholCorrection *= 100;
        if (cholCorrection > 1e6)
          break;
        call_ccfsg_xc (dciTrue, ScaleVertical);
        dmumps_c(&id);
        if (id.info[0] != -10)
          FailedMumps = dciFalse;;
      }
      dr.reset(nvar + nconI);
      //rhs(1:nvar+nconI)     = -A'*inv(A*A')*r
      //rhs(nvar+nconI+1:end) = inv(A*A')*r
      for (Int i = 0; i < nvar + nconI; i++)
        dr.get_doublex()[i] = rhs[i];
    } else {
      dr.solve (CHOLMOD_A, *LJ, r); // dr = inv(AA')r;
      dr.sdmult (*J, 1, mone, zero, dr);
    }

    return 0;
  }

  void Interface::scale_x (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xx[i], Li = blx[i], Ui = bux[i];
      if ( (PartialPenal) && (Li > -dciInf) && (Ui < dciInf) ) {
        if ( (zi - Li) < (Ui - zi) )
          Vx[i] *= Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          Vx[i] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
        continue;
      }
      if (Li > -dciInf)
        Vx[i] *= Min(MaxDiag, Max(MinDiag, zi - Li));
      if (Ui < dciInf)
        Vx[i] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = sx[i], Li = clx[ineqIdx[i]], Ui = cux[ineqIdx[i]];
      Int j = nvar + i;
      if ( (PartialPenal) && (Li > -dciInf) && (Ui < dciInf) ) {
        if ( (zi - Li) < (Ui - zi) )
          Vx[j] *= Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          Vx[j] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
        continue;
      }
      if (Li > -dciInf)
        Vx[j] *= Min(MaxDiag, Max(MinDiag, zi - Li));
      if (Ui < dciInf)
        Vx[j] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
    }
  }

  void Interface::scale_xc (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], Li = blx[i], Ui = bux[i];
      if ( (PartialPenal) && (Li > -dciInf) && (Ui < dciInf) ) {
        if ( (zi - Li) < (Ui - zi) )
          Vx[i] *= Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          Vx[i] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
        continue;
      }
      if (Li > -dciInf)
        Vx[i] *= Min(MaxDiag, Max(MinDiag, zi - Li));
      if (Ui < dciInf)
        Vx[i] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = scx[i], Li = clx[ineqIdx[i]], Ui = cux[ineqIdx[i]];
      Int j = nvar + i;
      if ( (PartialPenal) && (Li > -dciInf) && (Ui < dciInf) ) {
        if ( (zi - Li) < (Ui - zi) )
          Vx[j] *= Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          Vx[j] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
        continue;
      }
      if (Li > -dciInf)
        Vx[j] *= Min(MaxDiag, Max(MinDiag, zi - Li));
      if (Ui < dciInf)
        Vx[j] *= Min(MaxDiag, Max(MinDiag, Ui - zi));
    }
  }

  void Interface::updyineq () { //This is -mu*Penalization on obj fun
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i], xi = xx[i];
      yineqx[i] = 0;
      if ( (PartialPenal) && (bli > -dciInf) && (bui < dciInf) ) {
        if ( (xi - bli) < (bui - xi) )
          yineqx[i] -= mu/(xi - bli);
        else
          yineqx[i] -= mu/(xi - bui);
        continue;
      }
      if (bli > -dciInf)
        yineqx[i] -= mu/(xi - bli);
      if (bui < dciInf)
        yineqx[i] -= mu/(xi - bui);
    }
    for (Int i = 0; i < nconI; i++) {
      Real cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], si = sx[i];
      Int j = nvar + i;
      yineqx[j] = 0;
      if ( (PartialPenal) && (cli > -dciInf) && (cui < dciInf) ) {
        if ( (si - cli) < (cui - si) )
          yineqx[j] -= mu/(si - cli);
        else
          yineqx[j] -= mu/(si - cui);
        continue;
      }
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
      if ( (PartialPenal) && (bli > -dciInf) && (bui < dciInf) ) {
        if ( (xi - bli) < (bui - xi) )
          yineqx[i] -= mu/(xi - bli);
        else
          yineqx[i] -= mu/(xi - bui);
        continue;
      }
      if (bli > -dciInf)
        yineqx[i] -= mu/(xi - bli);
      if (bui < dciInf)
        yineqx[i] -= mu/(xi - bui);
    }
    for (Int i = 0; i < nconI; i++) {
      Real cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]], si = scx[i];
      Int j = nvar + i;
      yineqx[j] = 0;
      if ( (PartialPenal) && (cli > -dciInf) && (cui < dciInf) ) {
        if ( (si - cli) < (cui - si) )
          yineqx[j] -= mu/(si - cli);
        else
          yineqx[j] -= mu/(si - cui);
        continue;
      }
      if (cli > -dciInf)
        yineqx[j] -= mu/(si - cli);
      if (cui < dciInf)
        yineqx[j] -= mu/(si - cui);
    }
  }

  /* yineq = -mu*P(x) = theta_U - theta_L
   * theta_L, theta_U >= 0
   * theta_L = max(0.0, -yineq)
   * theta_U = max(0.0, yineq)
   *
   * (Z - L)*theta_L = 0
   * (U - Z)*theta_U = 0
   */
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

  /* -mu*P_s(s) - lambda_I = 0
   * yineq_s = -mu*P_s(s)
   */
  Real Interface::calc_ydif () {
    Real dif = 0;
    for (Int i = 0; i < nconI; i++)
      dif += fabs (yx[ineqIdx[i]] - yineqx[nvar + i]);
    return dif;
  }

  Real Interface::calc_pen () {
    Real val = 0.0;
    for (Int i = 0; i < nvar; i++) {
      if ( (PartialPenal) && (blx[i] > -dciInf) && (bux[i] < dciInf) ) {
        if ( (xx[i] - blx[i]) < (bux[i] - xx[i]) )
          val += log (xx[i] - blx[i]);
        else
          val += log (bux[i] - xx[i]);
        continue;
      }
      if (blx[i] > -dciInf)
        val += log (xx[i] - blx[i]);
      if (bux[i] < dciInf)
        val += log (bux[i] - xx[i]);
    }
    for (Int i = 0; i < nconI; i++) {
      if ( (PartialPenal) && (clx[ineqIdx[i]] > -dciInf) && (cux[ineqIdx[i]] < dciInf) ) {
        if ( (PartialPenal) && (sx[i] - clx[ineqIdx[i]]) < (cux[ineqIdx[i]] - sx[i]) )
          val += log (sx[i] - clx[ineqIdx[i]]);
        else
          val += log (cux[ineqIdx[i]] - sx[i]);
        continue;
      }
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
      if ( (PartialPenal) && (blx[i] > -dciInf) && (bux[i] < dciInf) ) {
        if ( (xcx[i] - blx[i]) < (bux[i] - xcx[i]) )
          val += log (xcx[i] - blx[i]);
        else
          val += log (bux[i] - xcx[i]);
        continue;
      }
      if (blx[i] > -dciInf)
        val += log (xcx[i] - blx[i]);
      if (bux[i] < dciInf)
        val += log (bux[i] - xcx[i]);
    }
    for (Int i = 0; i < nconI; i++) {
      if ( (PartialPenal) && (clx[ineqIdx[i]] > -dciInf) && (cux[ineqIdx[i]] < dciInf) ) {
        if ( (scx[i] - clx[ineqIdx[i]]) < (cux[ineqIdx[i]] - scx[i]) )
          val += log (scx[i] - clx[ineqIdx[i]]);
        else
          val += log (cux[ineqIdx[i]] - scx[i]);
        continue;
      }
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
