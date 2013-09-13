#include "interface.h"
#include <sstream>
#include <fstream>
#include <map>
#include <string>
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
    DeltaTiny = 1e-10;
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
    DisplayLevel = 1;
    VerboseLevel = 0;
    TableLevel = 0;
    MaxDiag = 1e20;
    MinDiag = 0;
    max_objfun_scale = 1e6;

    ReadParameters();
    //Mumps
//    if (UseMUMPS) {
//      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//      id.job=JOB_INIT;
//      id.par=1;
//      id.sym=2;
//      id.comm_fortran=USE_COMM_WORLD;
//      dmumps_c(&id);
//    }
  }

  void Interface::ReadParameters () {
    std::ifstream paramFile("dcicpp.spc");
    if (paramFile.fail())
      return;

    enum parameters {en_DeltaMax, en_maxrest, en_maxit, en_maxssmll, en_maxitSteih, 
      en_minitSteih, en_relitSteih, en_nfailv, en_csic, en_csig, en_rhomin, en_phi1,
      en_phi2, en_kappa1, en_kappa2, en_kappa3, en_kappa4, en_zeta1, en_zeta2, 
      en_zeta3, en_alphaR, en_alphaI, en_alphaS, en_eta1, en_eta2, en_eta3,
      en_DeltaMin, en_DeltaTiny, en_minstep, en_Delta0, en_thetaR, en_LbdMax,
      en_eps1, en_eps2, en_eps3, en_epsmu, en_epsgap, en_bfgsupd, en_c1, en_c2,
      en_MaxTime, en_minBk, en_UseCG, en_PartialPenal, en_project_dcp,
      en_project_bfgs, en_trustWorstdn, en_trustConvexBox, en_penal_trust,
      en_penal_bfgs, en_UseMUMPS, en_ScaleVertical, en_DisplayLevel,
      en_VerboseLevel, en_MaxDiag, en_MinDiag, en_UseVertInteriorPoint,
      en_UseVertSafeguard, en_UsePorcelli, en_UseObjfunScale, 
      en_objfun_count, en_choleskyCorrection, en_max_objfun_scale,
      en_UseVariableScaling, en_TableLevel, en_RebootOnVertFail
    };
    std::map<std::string, int> paramMap;

    paramMap["UseVertSafeguard"] = en_UseVertSafeguard;
    paramMap["choleskyCorrection"] = en_choleskyCorrection;
    paramMap["RebootOnVertFail"] = en_RebootOnVertFail;
    paramMap["UseVertInteriorPoint"] = en_UseVertInteriorPoint;
    paramMap["MaxDiag"] = en_MaxDiag;
    paramMap["MinDiag"] = en_MinDiag;
    paramMap["VerboseLevel"] = en_VerboseLevel;
    paramMap["DisplayLevel"] = en_DisplayLevel;
    paramMap["TableLevel"] = en_TableLevel;
    paramMap["ScaleVertical"] = en_ScaleVertical;
    paramMap["UseMUMPS"] = en_UseMUMPS;
    paramMap["UseCG"] = en_UseCG;
    paramMap["UsePorcelli"] = en_UsePorcelli;
    paramMap["UseObjfunScale"] = en_UseObjfunScale;
    paramMap["objfun_count"] = en_objfun_count;
    paramMap["max_objfun_scale"] = en_max_objfun_scale;
    paramMap["UseVariableScaling"] = en_UseVariableScaling;
    paramMap["PartialPenal"] = en_PartialPenal;
    paramMap["project_dcp"] = en_project_dcp;
    paramMap["project_bfgs"] = en_project_bfgs;
    paramMap["trustWorstdn"] = en_trustWorstdn;
    paramMap["trustConvexBox"] = en_trustConvexBox;
    paramMap["penal_trust"] = en_penal_trust;
    paramMap["penal_bfgs"] = en_penal_bfgs;
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
    paramMap["bfgsupd"] = en_bfgsupd;
    paramMap["c1"] = en_c1;
    paramMap["c2"] = en_c2;
    paramMap["MaxTime"] = en_MaxTime;
    paramMap["minBk"] = en_minBk;

    std::string param, value;
    while (getline(paramFile, param, ' ')) {
      getline(paramFile, value, '\n');
      std::stringstream aux;
      aux << value;
      switch (paramMap[param]) {
        case en_UseVertSafeguard: aux >> UseVertSafeguard; break;
        case en_RebootOnVertFail: aux >> RebootOnVertFail; break;
        case en_choleskyCorrection: aux >> choleskyCorrection; break;
        case en_UseVertInteriorPoint: aux >> UseVertInteriorPoint; break;
        case en_MaxDiag: aux >> MaxDiag; break;
        case en_MinDiag: aux >> MinDiag; break;
        case en_VerboseLevel: aux >> VerboseLevel; break;
        case en_DisplayLevel: aux >> DisplayLevel; break;
        case en_TableLevel: aux >> TableLevel; break;
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
        case en_bfgsupd: aux >> bfgsupd; break;
        case en_c1: aux >> c1; break;
        case en_c2: aux >> c2; break;
        case en_MaxTime: aux >> MaxTime; break;
        case en_minBk: aux >> minBk; break;
        case en_UseCG: aux >> UseCG; break;
        case en_UsePorcelli: aux >> UsePorcelli; break;
        case en_UseObjfunScale: aux >> UseObjfunScale; break;
        case en_objfun_count: aux >> objfun_count; break;
        case en_max_objfun_scale: aux >> max_objfun_scale; break;
        case en_UseVariableScaling: aux >> UseVariableScaling; break;
        case en_UseMUMPS: aux >> UseMUMPS; break;
        case en_PartialPenal: aux >> PartialPenal; break;
        case en_project_dcp: aux >> project_dcp; break;
        case en_project_bfgs: aux >> project_bfgs; break;
        case en_trustWorstdn: aux >> trustWorstdn; break;
        case en_trustConvexBox: aux >> trustConvexBox; break;
        case en_penal_trust: aux >> penal_trust; break;
        case en_penal_bfgs: aux >> penal_bfgs; break;
        case en_ScaleVertical: aux >> ScaleVertical; break;
      }
    }
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
    LimLbd = dciTrue;
    Lambda = 0;
    cholCorrection = 0;
    choleskyCorrection = 1e-6;
    cholFailed = dciFalse;
    objfun_count = 0;

    //Strategy choices
    UseCG = dciFalse;
    UsePorcelli = dciTrue;
    UseObjfunScale = dciTrue;
    UseVariableScaling = dciTrue;
    PartialPenal = dciTrue;
    project_dcp = dciFalse;
    project_dn = dciTrue;
    project_bfgs = dciTrue;
    trustWorstdn = dciFalse;
    trustConvexBox = dciFalse;
    penal_trust = dciTrue;
    penal_bfgs = dciTrue;
    UseMUMPS = dciFalse;
    ScaleVertical = dciTrue;
    UseVertInteriorPoint = dciTrue;
    UseVertSafeguard = dciFalse;
    RebootOnVertFail = dciTrue;

    //Program properties
    Ineq = dciFalse; //Has some inequalities
    Linear = dciFalse; //Only linear constraints
    Bounded = dciFalse; //Has bounds
  }

  void Interface::InitialValues () {

    maxitSteih = Max (minitSteih, Min(maxitSteih, Int(nvar*relitSteih+5) ) );
    minstep *= Min (csig, csic);

    std::vector<Int> fixed_vector;
    for (Int i = 0; i < nvar; i++) {
      Real bli = blx[i], bui = bux[i];
      VariableScaling[i] = 1.0;
      if ( (xx[i] > bli) && (xx[i] < bui) ) {
        xcx[i] = xx[i];
        continue;
      }
      if (bli - bui > dciTiny) {
        throw("Infeasible bounds");
      } else if (bli - bui > - dciTiny) {
        fixed_vector.push_back(i);
        xx[i] = (bli + bui)/2;
        xcx[i] = xx[i];
      } else {
        Real smldelta = Min ( 1e-2, (bui - bli)/100.0);
        xx[i] = Max ( Min ( xx[i], bui - smldelta ), bli + smldelta );
        xcx[i] = xx[i];
        assert (xx[i] > bli);
        assert (xx[i] < bui);
      }
    }

    nfix = fixed_vector.size();
    if (nfix > 0) {
      fixed_index = new Int[nfix];
      for (Int i = 0; i < nfix; i++) {
        fixed_index[i] = fixed_vector[i];
      }
    }

    // Calculating the function value and c without the barrier.
    objfun_scale = 1.0;
    call_fn ();
    call_ofg (dciTrue);
    if (UseObjfunScale)
      objfun_scale = Min( Max(Max(1.0, g->norm()), AbsValue(*f)), max_objfun_scale );

    if (ncon > 0) {
      for (Int i = 0; i < nconI; i++) {
        Real cxi = cx[ineqIdx[i]], cli = clx[ineqIdx[i]], cui = cux[ineqIdx[i]];
        Real smldelta = Min ( 1e-2, (cui - cli)/100);
        assert (smldelta > 0);
        sx[i] = Max ( Min ( cxi, cui - smldelta ), cli + smldelta );
        scx[i] = sx[i];
      }
      UpdateScaling_x();

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
      UpdateScaling_x();
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
