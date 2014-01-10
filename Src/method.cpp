#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  // From point x, project g so that l < x + g < u
  void Interface::projectBounds_x (Vector & v) {
    pReal vx = v.get_doublex();
    for (Int i = 0; i < nvar+nconI; i++) {
      Real alpha = 1.0;
      if (vx[i] > 0)
        alpha = Min(alpha, (u_bndx[i] - xx[i])/vx[i]);
      else if (vx[i] < 0)
        alpha = Min(alpha, (l_bndx[i] - xx[i])/vx[i]);
      if (alpha < 1)
        vx[i] *= alpha;
    }
  }

  void Interface::projectBounds_xc (Vector & v) {
    pReal vx = v.get_doublex();
    for (Int i = 0; i < nvar+nconI; i++) {
      Real alpha = 1.0;
      if (vx[i] > 0)
        alpha = Min(alpha, (u_bndx[i] - xcx[i])/vx[i]);
      else if (vx[i] < 0)
        alpha = Min(alpha, (l_bndx[i] - xcx[i])/vx[i]);
      if (alpha < 1)
        vx[i] *= alpha;
    }
  }

  Bool Interface::calcFeasibilityOpt () {
    if (ccifg == 0)
      return 0;
    if (feasOpt == 0)
      feasOpt = new Vector (*env, nvar);
    pReal fox = feasOpt->get_doublex ();
    Real ci = 0, cli = 0, cui = 0;
    Vector gradc (*env, nvar);
    pReal gcx = gradc.get_doublex ();
    Bool tmpTrue = dciTrue;

    feasOpt->reset(nvar, 0.0);
    
    for (Int i = 1; i <= ncon; i++) {
      if (!equatn[i-1])
        continue;
      (*ccifg) (&nvar, &i, xcx, &ci, gcx, &tmpTrue);
      cli = clx[i-1];
      cui = cux[i-1];

      // If it's an inequality and inside limits, do nothing
      if (ci >= cli && ci <= cui)
        continue;

      for (Int j = 0; j < nvar; j++)
        fox[j] += ci * gcx[j];
    }

    return (feasOpt->norm () < 1e-6);
  }

  void Interface::updateMultipliers () {
    if (ncon == 0) {
      *gp = *g;
      return;
    }
    Vector ytmp (*env);
    naProj (*g, *gp, ytmp);
    *y = ytmp;

    bool ForceSign = dciTrue;
    //Force sign
    if (!ForceSign)
      return;
    double coef = 100, expo = 0.5;
    double limiter = coef * pow (mu, expo);
    for (int i = 0; i < nconI; i++) {
      Real cli = l_bndx[nvar+i], cui = u_bndx[nvar+i];
      int j = ineq_index[i];
      if (cli <= -dciInf)
        yx[j] = Max (yx[j], -limiter);
      else if (cui >= dciInf)
        yx[j] = Min (yx[j], limiter);
      else {
        if (xcx[nvar+i] > (cli+cui)/2)
          yx[j] = Max (yx[j], -limiter);
        else
          yx[j] = Min (yx[j], limiter);
      }
    }
  }

  void Interface::leastSquaresCG (Bool transp, const Vector & rhs, Vector & sol, Vector & res) {
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

  void Interface::linearSystemCG (Bool transp, const Vector & rhs, Vector & sol) {
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

  void Interface::updateMu () {
    if (nconI == 0) {
      mu = 0;
      return;
    }
    Real mufactor = 1.0;
    Real minmuthisiter = Max (0.1*mu, 1e-24);
    Real minother = 100*Min (rho, rho*rho); 
//    minother = Min (minother, calcYdif() + lagrgap + infacgap);
    minother = Min (minother, mufactor * mu);
    minother = Min (minother, 100*normck);
    Real dotls = 0.0;
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar+i;
      Real si = xcx[j], li = l_bndx[j], ui = u_bndx[j];
      if (li > -dciInf)
        dotls += Max(0.0, -yx[i]) * (si - li);
      if (ui < dciInf)
        dotls += Max(0.0,  yx[i]) * (ui - si);
    }
    dotls /= nconI;
    minother = Min (minother, dotls);
    mu = Max (minother, minmuthisiter);
    mu = minother;
  }

  void Interface::naProjApprox (Vector & rin, Vector & Pr, Vector & tmp) {
    Real mone[2] = {-1,0}, zero[2] = {0,0}, one[2] = {1,0};

    Sparse M(*J);
    M.band_inplace(-1, 1, 1);
    Factor Mfac(*env);
    Mfac.analyze(M);
    static Real facCorr = 1e-6;
    Mfac.factorize(M, facCorr);
    if (!env->IsPosDef()) {
      facCorr = 1e-3;
      Mfac.factorize(M, facCorr);
    }
#ifdef VERBOSE
    if (verbosity_level > 2) {
      M.print_more();
    }
#endif

    Pr.sdmult (*J, 0, mone, zero, rin);
    tmp.reset (Pr.size(), 0);
    Vector r(*env), p(*env), Jtp(*env), JJtp (*env), t(*env);
    r.scale(Pr, -1);
    // M*t = r
    t.solve (CHOLMOD_A, Mfac, r);
    //t = r
    p.scale(t, -1);
    Jtp.sdmult (*J, 1, one, zero, p);
    JJtp.sdmult (*J, 0, one, zero, Jtp);
    Real rr, rrp, alpha, beta;
    Real rr0;
    rrp = r.dot(t);
    if (rrp < 1e-12) {
      Pr = rin;
      return;
    }
    rr0 = rrp;
    Int cgiter = 1;
    while ( (rrp > ngp*rr0) && (cgiter < nvar + nconI) ) {
      rr = rrp;
      Real gamma = Jtp.dot(Jtp);
      if (gamma < 1e-12) {
        break;
      }
      alpha = rr/gamma;
      tmp.saxpy(p, alpha);
      r.saxpy(JJtp, alpha);
//      t = r;
      t.solve (CHOLMOD_A, Mfac, r);
      rrp = r.dot(t);
      beta = rrp/rr;
      p.scale(beta);
      p.saxpy(t, -1);
      Jtp.sdmult (*J, 1, one, zero, p);
      JJtp.sdmult (*J, 0, one, zero, Jtp);
      cgiter++;
    }
    Pr = rin;
    Pr.sdmult (*J, 1, one, one, tmp);
  }

  void Interface::naProj (Vector & r, Vector & Pr, Vector & tmp) {
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
    Pr = r;
    Pr.sdmult (*J, 1, one, one, tmp); //Pr = Pr + J'*tmp => Pr = r - A'*inv(A*A')*A*r
  }

  Int Interface::naStep (Vector & r, Vector & dr) { //dr = -A'*inv(AA')*r
    Real mone[2] = {-1,0}, zero[2] = {0,0};
    if (use_conjugate_gradient) {
      linearSystemCG (0, r, dr);
      dr.sdmult (*J, 1, mone, zero, dr);
    } else {
      dr.solve (CHOLMOD_A, *LJ, r); // dr = inv(AA')r;
      dr.sdmult (*J, 1, mone, zero, dr);
    }

    return 0;
  }

  void Interface::updateScaling_x () {
    if (scaling_matrix == 0)
      scaling_matrix = new Real[nvar + nconI];

    for (Int i = 0; i < nvar+nconI; i++) {
      Real zi = xx[i], Li = l_bndx[i], Ui = u_bndx[i];
      if (Li > Ui - dciTiny) {
        scaling_matrix[i] = 1.0;
        continue;
      }

      if ( Li > -dciInf && Ui < dciInf ) {
        if ( zi < (Li + Ui)/2 )
          scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      } else if (Li > -dciInf)
        scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, zi - Li));
      else if (Ui < dciInf)
        scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      else
        scaling_matrix[i] = 1.0;

      if (i < nvar)
        scaling_matrix[i] *= variable_scaling[i];
    }
  }

  void Interface::scale_x (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar + nconI; i++)
      Vx[i] *= scaling_matrix[i];
  }

  void Interface::updateScaling_xc () {
    if (scaling_matrix == 0)
      scaling_matrix = new Real[nvar + nconI];

    for (Int i = 0; i < nvar+nconI; i++) {
      Real zi = xcx[i], Li = l_bndx[i], Ui = u_bndx[i];
      if (Li > Ui - dciTiny) {
        scaling_matrix[i] = 1.0;
        continue;
      }
      if ( Li > -dciInf && Ui < dciInf ) {
        if ( zi < (Li + Ui)/2 )
          scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      } else if (Li > -dciInf)
        scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, zi - Li));
      else if (Ui < dciInf)
        scaling_matrix[i] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      else
        scaling_matrix[i] = 1.0;

      if (i < nvar)
        scaling_matrix[i] *= variable_scaling[i];
    }
  }

  void Interface::scale_xc (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar + nconI; i++)
      Vx[i] *= scaling_matrix[i];
  }

  void Interface::update_yineq () { //This is -mu*Penalization on obj fun
    for (Int i = 0; i < nvar+nconI; i++) {
      Real bli = l_bndx[i], bui = u_bndx[i], xi = xx[i];
      yineqx[i] = 0;
      if (bli - bui > - dciTiny)
        continue;
      if ( (partial_penalization) && (bli > -dciInf) && (bui < dciInf) ) {
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
  }

  /* yineq = -mu*P(x) = theta_U - theta_L
   * theta_L, theta_U >= 0
   * theta_L = max(0.0, -yineq)
   * theta_U = max(0.0, yineq)
   *
   * (Z - L)*theta_L = 0
   * (U - Z)*theta_U = 0
   */
  Real Interface::calcGap () {
    Real gap = 0;
    for (Int i = 0; i < nvar+nconI; i++)  {
      Real xi = xx[i], bli = l_bndx[i], bui = u_bndx[i];
      if ( (bli <= -dciInf) && (bui >= dciInf) )
        continue;
      if (bli > -dciInf)
        gap += Max(0.0, -yineqx[i]) * (xi - bli);
      if (bui < dciInf)
        gap += Max(0.0,  yineqx[i]) * (bui - xi);
    }
    gap /= (nvar + nconI);
    
    return gap;
  }

  /* -mu*P_s(s) - lambda_I = 0
   * yineq_s = -mu*P_s(s)
   */
  Real Interface::calcYdif () {
    Real dif = 0;
    for (Int i = 0; i < nconI; i++)
      dif += fabs (yx[ineq_index[i]] - yineqx[nvar + i]);
    return dif;
  }

  Real Interface::calcPen () {
    Real val = 0.0;
    for (Int i = 0; i < nvar+nconI; i++) {
      if (l_bndx[i] - u_bndx[i] > - dciTiny)
        continue;
      if ( (partial_penalization) && 
           (l_bndx[i] > -dciInf) && (u_bndx[i] < dciInf) ) {
        if ( (xx[i] - l_bndx[i]) < (u_bndx[i] - xx[i]) )
          val += log (xx[i] - l_bndx[i]);
        else
          val += log (u_bndx[i] - xx[i]);
        continue;
      }
      if (l_bndx[i] > -dciInf)
        val += log (xx[i] - l_bndx[i]);
      if (u_bndx[i] < dciInf)
        val += log (u_bndx[i] - xx[i]);
    }
    return val;
  }

  Real Interface::calcPen_xc () {
    Real val = 0.0;
    for (Int i = 0; i < nvar+nconI; i++) {
      if (l_bndx[i] - u_bndx[i] > - dciTiny)
        continue;
      if ( (partial_penalization) && (l_bndx[i] > -dciInf) && (u_bndx[i] < dciInf) ) {
        if ( (xcx[i] - l_bndx[i]) < (u_bndx[i] - xcx[i]) )
          val += log (xcx[i] - l_bndx[i]);
        else
          val += log (u_bndx[i] - xcx[i]);
        continue;
      }
      if (l_bndx[i] > -dciInf)
        val += log (xcx[i] - l_bndx[i]);
      if (u_bndx[i] < dciInf)
        val += log (u_bndx[i] - xcx[i]);
    }
    return val;
  }

}
