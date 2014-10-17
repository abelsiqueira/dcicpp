#include "interface.h"
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <map>
#include <string>
//#include <cassert>

namespace DCI {

  void Interface::defineParameters () {
    DeltaMax = 1e6;
    maxrest = 200000;
    maxit = 200000;
    maxssmll = 5;
    maxitSteih = 100000;
    minitSteih = 100;
    relitSteih = 10;
    nfailv = 5;
    csic = 1e-6;
    csig = 1e-6;
    rhomin = 1e-8;
    phi1 = 1.0;
    phi2 = 0.99;
    kappa1 = 1e-4;
    kappa2 = 0.25;
    kappa3 = 0.7;
    kappa4 = 2.5;
    beta1 = 0.1;
    beta2 = 0.25;
    zeta1 = 2;
    zeta2 = 1;
    zeta3 = 5;
    alphaR = 0.75;
    alphaI = 2.5;
    alphaS = 6.25e-2;
    eta1 = 1e-3;
    eta2 = 0.2;
    eta3 = 0.9;
    DeltaMin = 1e-4;
    DeltaTiny = 1e-10;
    minstep = 1e-3;
    Delta0 = 1e6;
    thetaR = 0.99;
    LbdMax = 1e6;
    eps1 = 1e-6;
    eps2 = 1e-14;
    eps3 = 1e-8;
    epsmu = 1e-6;
    epsgap = 1e-4;
    c1 = 0.5;
    c2 = 5e-1;
    max_time = 300; // 5 minutes
    minBk = 1e-12;
    display_level = 1;
    debug_level = 0;
    verbosity_level = 0;
    table_print_level = 0;
    print_A_at_end = false;
    MaxDiag = 1e9;
    MinDiag = 0;
    max_objective_scaling = 1e6;
    max_constraint_scaling = 1e6;
    max_variable_scaling = 1e6;
    nvarshowmax = 10;
    nconshowmax = 10;

    readParameters();
  }

  void Interface::readParameters () {
    std::ifstream paramFile("dcicpp.spc");
    bool has_param_file = !(paramFile.fail());

    enum parameters {en_DeltaMax, en_maxrest, en_maxit, en_maxssmll,
      en_maxitSteih, en_minitSteih, en_relitSteih, en_nfailv, en_csic, en_csig,
      en_rhomin, en_phi1, en_phi2, en_kappa1, en_kappa2, en_kappa3, en_kappa4,
      en_zeta1, en_zeta2, en_beta1, en_beta2, en_zeta3, en_alphaR, en_alphaI,
      en_alphaS, en_eta1, en_eta2, en_eta3, en_DeltaMin, en_DeltaTiny,
      en_minstep, en_Delta0, en_thetaR, en_LbdMax, en_eps1, en_eps2, en_eps3,
      en_epsmu, en_epsgap, en_c1, en_c2, en_max_time, en_minBk,
      en_use_conjugate_gradient, en_partial_penalization, en_project_dcp,
      en_trustWorstdn, en_trustConvexBox, en_penal_trust,
      en_scale_normal, en_display_level, en_debug_level,
      en_verbosity_level, en_print_A_at_end, en_MaxDiag, en_MinDiag,
      en_use_objective_scaling, en_objfun_count, en_use_constraint_scaling,
      en_max_objective_scaling, en_use_variable_scaling, en_table_print_level,
      en_max_constraint_scaling, en_max_variable_scaling, en_use_soc,
      en_use_normal_safe_guard, en_nvarshowmax, en_nconshowmax,
      en_normal_fail_reboot, en_chol_correction_increase,
      en_cholesky_base_correction, en_infeasibility_tol, en_use_lsmr
    };
    std::map<std::string, int> paramMap;

    paramMap["normal_fail_reboot"] = en_normal_fail_reboot;
    paramMap["chol_correction_increase"] = en_chol_correction_increase;
    paramMap["cholesky_base_correction"] = en_cholesky_base_correction;
    paramMap["MaxDiag"] = en_MaxDiag;
    paramMap["MinDiag"] = en_MinDiag;
    paramMap["debug_level"] = en_debug_level;
    paramMap["verbosity_level"] = en_verbosity_level;
    paramMap["display_level"] = en_display_level;
    paramMap["table_print_level"] = en_table_print_level;
    paramMap["print_A_at_end"] = en_print_A_at_end;
    paramMap["scale_normal"] = en_scale_normal;
    paramMap["use_conjugate_gradient"] = en_use_conjugate_gradient;
    paramMap["use_objective_scaling"] = en_use_objective_scaling;
    paramMap["use_soc"] = en_use_soc;
    paramMap["use_lsmr"] = en_use_lsmr;
    paramMap["use_normal_safe_guard"] = en_use_normal_safe_guard;
    paramMap["use_constraint_scaling"] = en_use_constraint_scaling;
    paramMap["objfun_count"] = en_objfun_count;
    paramMap["infeasibility_tol"] = en_infeasibility_tol;
    paramMap["max_objective_scaling"] = en_max_objective_scaling;
    paramMap["max_constraint_scaling"] = en_max_constraint_scaling;
    paramMap["max_variable_scaling"] = en_max_variable_scaling;
    paramMap["use_variable_scaling"] = en_use_variable_scaling;
    paramMap["partial_penalization"] = en_partial_penalization;
    paramMap["project_dcp"] = en_project_dcp;
    paramMap["trustWorstdn"] = en_trustWorstdn;
    paramMap["trustConvexBox"] = en_trustConvexBox;
    paramMap["penal_trust"] = en_penal_trust;
    paramMap["DeltaMax"] = en_DeltaMax;
    paramMap["maxrest"] = en_maxrest;
    paramMap["maxit"] = en_maxit;
    paramMap["maxssmll"] = en_maxssmll;
    paramMap["maxitSteih"] = en_maxitSteih;
    paramMap["minitSteih"] = en_minitSteih;
    paramMap["relitSteih"] = en_relitSteih;
    paramMap["nfailv"] = en_nfailv;
    paramMap["csic"] = en_csic;
    paramMap["csig"] = en_csig;
    paramMap["rhomin"] = en_rhomin;
    paramMap["phi1"] = en_phi1;
    paramMap["phi2"] = en_phi2;
    paramMap["kappa1"] = en_kappa1;
    paramMap["kappa2"] = en_kappa2;
    paramMap["kappa3"] = en_kappa3;
    paramMap["kappa4"] = en_kappa4;
    paramMap["beta1"] = en_beta1;
    paramMap["beta2"] = en_beta2;
    paramMap["zeta1"] = en_zeta1;
    paramMap["zeta2"] = en_zeta2;
    paramMap["zeta3"] = en_zeta3;
    paramMap["alphaR"] = en_alphaR;
    paramMap["alphaI"] = en_alphaI;
    paramMap["alphaS"] = en_alphaS;
    paramMap["eta1"] = en_eta1;
    paramMap["eta2"] = en_eta2;
    paramMap["eta3"] = en_eta3;
    paramMap["DeltaMin"] = en_DeltaMin;
    paramMap["DeltaTiny"] = en_DeltaTiny;
    paramMap["minstep"] = en_minstep;
    paramMap["Delta0"] = en_Delta0;
    paramMap["thetaR"] = en_thetaR;
    paramMap["LbdMax"] = en_LbdMax;
    paramMap["eps1"] = en_eps1;
    paramMap["eps2"] = en_eps2;
    paramMap["eps3"] = en_eps3;
    paramMap["epsmu"] = en_epsmu;
    paramMap["epsgap"] = en_epsgap;
    paramMap["c1"] = en_c1;
    paramMap["c2"] = en_c2;
    paramMap["max_time"] = en_max_time;
    paramMap["minBk"] = en_minBk;
    paramMap["nvarshowmax"] = en_nvarshowmax;
    paramMap["nconshowmax"] = en_nconshowmax;

    if (!has_param_file)
      return;

    std::string param, value;

    if (debug_level > 0)
      std::cout << "Parameters loaded:" << std::endl;

    while (getline(paramFile, param, ' ')) {
      getline(paramFile, value, '\n');
      std::stringstream aux;
      aux << value;
      int choice;
      try {
        choice = paramMap.at(param);
      } catch (const std::out_of_range& ex) {
        std::cerr << "Option " << param << " not recognized. Ignoring" <<
          std::endl;
        continue;
      }
      if (debug_level > 0) {
        std::cout << param << " = " << value << std::endl;
      }

      switch (choice) {
        case en_normal_fail_reboot: aux >> normal_fail_reboot; break;
        case en_chol_correction_increase: aux >> chol_correction_increase; break;
        case en_cholesky_base_correction: aux >> cholesky_base_correction; break;
        case en_MaxDiag: aux >> MaxDiag; break;
        case en_MinDiag: aux >> MinDiag; break;
        case en_debug_level: aux >> debug_level; break;
        case en_verbosity_level: aux >> verbosity_level; break;
        case en_display_level: aux >> display_level; break;
        case en_table_print_level: aux >> table_print_level; break;
        case en_print_A_at_end: aux >> print_A_at_end; break;
        case en_DeltaMax: aux >> DeltaMax; break;
        case en_maxrest: aux >> maxrest; break;
        case en_maxit: aux >> maxit; break;
        case en_maxssmll: aux >> maxssmll; break;
        case en_maxitSteih: aux >> maxitSteih; break;
        case en_minitSteih: aux >> minitSteih; break;
        case en_relitSteih: aux >> relitSteih; break;
        case en_nfailv: aux >> nfailv; break;
        case en_csic: aux >> csic; break;
        case en_csig: aux >> csig; break;
        case en_rhomin: aux >> rhomin; break;
        case en_phi1: aux >> phi1; break;
        case en_phi2: aux >> phi2; break;
        case en_kappa1: aux >> kappa1; break;
        case en_kappa2: aux >> kappa2; break;
        case en_kappa3: aux >> kappa3; break;
        case en_kappa4: aux >> kappa4; break;
        case en_beta1: aux >> beta1; break;
        case en_beta2: aux >> beta2; break;
        case en_zeta1: aux >> zeta1; break;
        case en_zeta2: aux >> zeta2; break;
        case en_zeta3: aux >> zeta3; break;
        case en_alphaR: aux >> alphaR; break;
        case en_alphaI: aux >> alphaI; break;
        case en_alphaS: aux >> alphaS; break;
        case en_eta1: aux >> eta1; break;
        case en_eta2: aux >> eta2; break;
        case en_eta3: aux >> eta3; break;
        case en_DeltaMin: aux >> DeltaMin; break;
        case en_DeltaTiny: aux >> DeltaTiny; break;
        case en_minstep: aux >> minstep; break;
        case en_Delta0: aux >> Delta0; break;
        case en_thetaR: aux >> thetaR; break;
        case en_LbdMax: aux >> LbdMax; break;
        case en_eps1: aux >> eps1; break;
        case en_eps2: aux >> eps2; break;
        case en_eps3: aux >> eps3; break;
        case en_epsmu: aux >> epsmu; break;
        case en_epsgap: aux >> epsgap; break;
        case en_c1: aux >> c1; break;
        case en_c2: aux >> c2; break;
        case en_max_time: aux >> max_time; break;
        case en_minBk: aux >> minBk; break;
        case en_use_conjugate_gradient: aux >> use_conjugate_gradient; break;
        case en_use_objective_scaling: aux >> use_objective_scaling; break;
        case en_use_soc: aux >> use_soc; break;
        case en_use_lsmr: aux >> use_lsmr; break;
        case en_use_normal_safe_guard: aux >> use_normal_safe_guard; break;
        case en_use_constraint_scaling: aux >> use_constraint_scaling; break;
        case en_objfun_count: aux >> objfun_count; break;
        case en_infeasibility_tol: aux >> infeasibility_tol; break;
        case en_max_objective_scaling: aux >> max_objective_scaling; break;
        case en_max_constraint_scaling: aux >> max_constraint_scaling; break;
        case en_max_variable_scaling: aux >> max_variable_scaling; break;
        case en_use_variable_scaling: aux >> use_variable_scaling; break;
        case en_partial_penalization: aux >> partial_penalization; break;
        case en_project_dcp: aux >> project_dcp; break;
        case en_trustWorstdn: aux >> trustWorstdn; break;
        case en_trustConvexBox: aux >> trustConvexBox; break;
        case en_penal_trust: aux >> penal_trust; break;
        case en_scale_normal: aux >> scale_normal; break;
        case en_nvarshowmax: aux >> nvarshowmax; break;
        case en_nconshowmax: aux >> nconshowmax; break;
      }
    }
  }

  void Interface::initialization () {
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
    l_bnd = 0;
    u_bnd = 0;
    y = 0;
    yineq = 0;
    cl = 0;
    cu = 0;
    sols = 0;
    xc = 0;
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
    ineq_index = 0;
    current_time = 0;
    max_time = dciInf;
    DLH = dciInf;
    DLV = 0;
    Lref = dciInf;
    start_at_one = dciFalse;
    initialized = dciFalse;
    running = dciFalse;
    solved = dciFalse;
    is_unlimited = dciFalse;
    solved = 0;
    NormalFlag = 0;
    iter = 0;
    itssmll = 0;
    tRest = 0;
    tSteih = 0;
    tRej = 0;
    tSoc = 0;
    cholFacs = 0;
    nSoc = 0;
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
    LimLbd = dciTrue;
    scaling_matrix = 0;
    cholesky_correction = 0;
    cholesky_failed = dciFalse;
    objfun_count = 10;
    infeasibility_tol = 1e-6;

    //Strategy choices
    use_conjugate_gradient = dciFalse;
    use_objective_scaling = dciTrue;
    use_soc = dciTrue;
    use_soc = dciFalse;
    use_normal_safe_guard = dciTrue;
    use_constraint_scaling = dciTrue;
    use_variable_scaling = dciTrue;
    partial_penalization = dciTrue;
    project_dcp = dciFalse;
    project_dn = dciTrue;
    trustWorstdn = dciFalse;
    trustConvexBox = dciFalse;
    penal_trust = dciFalse;
    scale_normal = dciFalse;
    normal_fail_reboot = dciTrue;
    chol_correction_increase = 10;
    cholesky_base_correction = 1e-12;

    //Program properties
    has_ineq = dciFalse; //Has some inequalities
    is_linear = dciFalse; //Only linear constraints
    is_bounded = dciFalse; //Has bounds
  }

  void Interface::initialValues () {

    maxitSteih = Max (minitSteih, Min(maxitSteih, Int(nvar*relitSteih+5) ) );
    minstep *= Min (csig, csic);

    if (use_variable_scaling)
      variable_scaling = new Real[nvar];
    for (Int i = 0; i < nvar; i++) {
      Real bli = l_bndx[i], bui = u_bndx[i];
      variable_scaling[i] = 1.0;
      if ( (xx[i] > bli) && (xx[i] < bui) ) {
        xcx[i] = xx[i];
        continue;
      }
      if (bli - bui > dciTiny) {
        throw("Infeasible bounds");
      } else if (bli - bui > - dciTiny) {
        throw("Fixed variables not allowed.\nUse http://github.com/abelsiqueira/nope");
      } else {
        Real smldelta = Min ( 1e-2, (bui - bli)/100.0);
        xx[i] = Max ( Min ( xx[i], bui - smldelta ), bli + smldelta );
        xcx[i] = xx[i];
        assert (xx[i] > bli);
        assert (xx[i] < bui);
      }
    }

    // Calculating the function value and c without the barrier.
    objective_scaling = 1.0;
    if (ncon > 0) {
      constraint_scaling = new Real[ncon];
      for (Int i = 0; i < ncon; i++)
        constraint_scaling[i] = 1.0;
    }
    call_fn ();
    call_ofg (dciTrue);
    if (use_objective_scaling)
      objective_scaling = Min( Max(Max(1.0, g->norm()), AbsValue(*f)),
          max_objective_scaling );
    if (use_constraint_scaling && ncon > 0) {
      for (Int i = 0; i < ncon; i++) {
        constraint_scaling[i] = Min( Max(1.0, AbsValue(cx[i])),
            max_constraint_scaling);
      }
      for (Int i = 0; i < nconI; i++) {
        Int j = nvar + i;
        if (l_bndx[j] > -dciInf)
          l_bndx[j] /= constraint_scaling[ineq_index[i]];
        if (u_bndx[j] < dciInf)
          u_bndx[j] /= constraint_scaling[ineq_index[i]];
      }
    }
    call_fn();

    if (ncon > 0) {
      for (Int i = 0; i < nconI; i++) {
        Real cxi = cx[ineq_index[i]], cli = l_bndx[nvar+i], cui = u_bndx[nvar+i];
        Real smldelta = Min ( 1e-2, (cui - cli)/100);
        assert (smldelta > 0);
        Int j = nvar + i;
        xx[j] = Max ( Min ( cxi, cui - smldelta ), cli + smldelta );
        xcx[j] = xx[j];
      }
      updateScaling_x();

      running = dciTrue;
      // Now, adding s.
      call_ccfsg ();
      normc = c->norm ();
      call_ofg ();

      this->analyzeJacobian ();
      this->cholesky ();

      if (debug_level > 0) {
        full(*J).print_more();
      }

      if (y == 0) {
        y = new Vector (*env, ncon);
        yx = y->get_doublex();
        for (Int i = 0; i < ncon; i++)
          yx[i] = 0;
      }
      updateMultipliers ();
//      gp->scale (*g, 1);
//      gp->sdmult (*J, 1, one, one, *y);
      Ln = *f + y->dot (*c);
    } else { //No constraints, may have bounds on the variables
      running = dciTrue; //If there are bounds, this will get them.
      updateScaling_x();
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

    running = dciFalse;
  }

}
