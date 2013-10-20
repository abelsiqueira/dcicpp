#include "interface.h"
//#include <cassert>
#include <cmath>

namespace DCI {
  // From point x, project g so that l < x + g < u
  void Interface::projectBounds_x (Vector & v) {
    pReal vx = v.get_doublex();
    for (Int i = 0; i < nvar; i++) {
      Real alpha = 1.0;
      if (vx[i] > 0)
        alpha = Min(alpha, (bux[i] - xx[i])/vx[i]);
      else if (vx[i] < 0)
        alpha = Min(alpha, (blx[i] - xx[i])/vx[i]);
      if (alpha < 1)
        vx[i] *= alpha;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real alpha = 1.0;
      if (vx[j] > 0)
        alpha = Min(alpha, (cux[i] - sx[i])/vx[j]);
      else if (vx[j] < 0)
        alpha = Min(alpha, (clx[i] - sx[i])/vx[j]);
      if (alpha < 1)
        vx[j] *= alpha;
    }
  }

  void Interface::projectBounds_xc (Vector & v) {
    pReal vx = v.get_doublex();
    for (Int i = 0; i < nvar; i++) {
      Real alpha = 1.0;
      if (vx[i] > 0)
        alpha = Min(alpha, (bux[i] - xcx[i])/vx[i]);
      else if (vx[i] < 0)
        alpha = Min(alpha, (blx[i] - xcx[i])/vx[i]);
      if (alpha < 1)
        vx[i] *= alpha;
    }
    for (Int i = 0; i < nconI; i++) {
      Int j = nvar + i;
      Real alpha = 1.0;
      if (vx[j] > 0)
        alpha = Min(alpha, (cux[i] - scx[i])/vx[j]);
      else if (vx[j] < 0)
        alpha = Min(alpha, (clx[i] - scx[i])/vx[j]);
      if (alpha < 1)
        vx[j] *= alpha;
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
      (*ccifg) (&nvar, &i, xcx, &ci, gcx, &tmpTrue);
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
      int j = ineq_index[i];
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
      Real si = scx[i], cli = clx[ineq_index[i]], cui = cux[ineq_index[i]];
      if (cli > -dciInf)
        dotls += Max(0.0, -yx[i]) * (si - cli);
      if (cui < dciInf)
        dotls += Max(0.0,  yx[i]) * (cui - si);
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

    for (Int i = 0; i < nvar; i++) {
      Real zi = xx[i], Li = blx[i], Ui = bux[i];
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

      scaling_matrix[i] *= variable_scaling[i];
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = sx[i], Li = clx[ineq_index[i]], Ui = cux[ineq_index[i]];
      Int j = nvar + i;
      if ( Li > -dciInf && Ui < dciInf ) {
        if ( zi < (Li + Ui)/2 )
          scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      } else if (Li > -dciInf)
        scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, zi - Li));
      else if (Ui < dciInf)
        scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      else
        scaling_matrix[j] = 1.0;
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

    for (Int i = 0; i < nvar; i++) {
      Real zi = xcx[i], Li = blx[i], Ui = bux[i];
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

      scaling_matrix[i] *= variable_scaling[i];
    }
    for (Int i = 0; i < nconI; i++) {
      Real zi = scx[i], Li = clx[ineq_index[i]], Ui = cux[ineq_index[i]];
      Int j = nvar + i;
      if ( Li > -dciInf && Ui < dciInf ) {
        if ( zi < (Li + Ui)/2 )
          scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, zi - Li));
        else
          scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      } else if (Li > -dciInf)
        scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, zi - Li));
      else if (Ui < dciInf)
        scaling_matrix[j] = Min(MaxDiag, Max(MinDiag, Ui - zi));
      else
        scaling_matrix[j] = 1.0;
    }
  }

  void Interface::scale_xc (Vector & V) {
    pReal Vx = V.get_doublex ();

    for (Int i = 0; i < nvar + nconI; i++)
      Vx[i] *= scaling_matrix[i];
  }

  void Interface::update_yineq () { //This is -mu*Penalization on obj fun
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i], xi = xx[i];
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
    for (Int i = 0; i < nconI; i++) {
      Real cli = clx[ineq_index[i]], cui = cux[ineq_index[i]], si = sx[i];
      Int j = nvar + i;
      yineqx[j] = 0;
      if ( (partial_penalization) && (cli > -dciInf) && (cui < dciInf) ) {
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
  Real Interface::calcGap () {
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
      Real si = sx[i], cli = clx[ineq_index[i]], cui = cux[ineq_index[i]];
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
  Real Interface::calcYdif () {
    Real dif = 0;
    for (Int i = 0; i < nconI; i++)
      dif += fabs (yx[ineq_index[i]] - yineqx[nvar + i]);
    return dif;
  }

  Real Interface::calcPen () {
    Real val = 0.0;
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] - bux[i] > - dciTiny)
        continue;
      if ( (partial_penalization) && (blx[i] > -dciInf) && (bux[i] < dciInf) ) {
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
      if ( (partial_penalization) && (clx[ineq_index[i]] > -dciInf) && (cux[ineq_index[i]] < dciInf) ) {
        if ( (partial_penalization) && (sx[i] - clx[ineq_index[i]]) < (cux[ineq_index[i]] - sx[i]) )
          val += log (sx[i] - clx[ineq_index[i]]);
        else
          val += log (cux[ineq_index[i]] - sx[i]);
        continue;
      }
      if (clx[ineq_index[i]] > -dciInf)
        val += log (sx[i] - clx[ineq_index[i]]);
      if (cux[ineq_index[i]] < dciInf)
        val += log (cux[ineq_index[i]] - sx[i]);
    }
    return val;
  }

  Real Interface::calcPen_xc () {
    Real val = 0.0;
    for (Int i = 0; i < nvar; i++) {
      if (blx[i] - bux[i] > - dciTiny)
        continue;
      if ( (partial_penalization) && (blx[i] > -dciInf) && (bux[i] < dciInf) ) {
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
      if ( (partial_penalization) && (clx[ineq_index[i]] > -dciInf) && (cux[ineq_index[i]] < dciInf) ) {
        if ( (scx[i] - clx[ineq_index[i]]) < (cux[ineq_index[i]] - scx[i]) )
          val += log (scx[i] - clx[ineq_index[i]]);
        else
          val += log (cux[ineq_index[i]] - scx[i]);
        continue;
      }
      if (clx[ineq_index[i]] > -dciInf)
        val += log (scx[i] - clx[ineq_index[i]]);
      if (cux[ineq_index[i]] < dciInf)
        val += log (cux[ineq_index[i]] - scx[i]);
    }
    return val;
  }

}
