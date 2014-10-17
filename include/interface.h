#ifndef dci_interface_h
#define dci_interface_h

#include "triplet.h"
#include "lsmr.h"

// For LSMR
void Aprod (int *, int *, double *, double *);
void Atprod (int *, int *, double *, double *);

namespace DCI {
  using namespace base_matrices;
  class Interface {
    public:
      // Constructors and Destructor
      Interface ();
      ~Interface ();

      // Problem functions
      int start ();
      int solve ();
      void show (std::ostream & = std::cout);
      void printLatex (char * = 0) const;
      void cuter () { start_at_one = dciTrue; cuter_status = 0;};

      // Sine quibus non
      void unc_setup (Int, Real *, Real *, Real *);
      void con_setup (Int, Real *, Real *, Real *, Int, Real *, Real *,
          Real *, Bool *);
      void set_linear (Int n, Bool *V);

      Real * get_x () const { return xx; }
      Real   get_f () const { return *f; }

      // Internal set nvar, ncon, nmax, mmax
      void set_nvar (Int n) { nvar = n; };
      void set_ncon (Int m) { ncon = m; };
      void set_nmax (Int n) { nmax = n; };
      void set_mmax (Int m) { mmax = m; };
      void set_amax (Int n) { amax = n; };

      // Function Assignment
      // Unconstrained
      void set_ufn (pufn p) { ufn = p; }
      void set_uofg (puofg p) { uofg = p; };
      void set_uprod (puprod p) { uprod = p; };
      void set_unames (punames p) { unames = p; };
      // Constrained
      void set_cfn (pcfn p) { cfn = p; };
      void set_cofg (pcofg p) { cofg = p; };
      void set_ccfsg (pccfsg p) { ccfsg = p; };
      void set_cprod (pcprod p) { cprod = p; };
      void set_ccifg (pccifg p) { ccifg = p; };
      void set_cnames (pcnames p) { cnames = p; };

      // LSMR
      void JacobMult(bool, double *, double *);
      friend void Aprod(int *, int *, double *, double *);
      friend void Atprod(int *, int *, double *, double *);

    protected:
      void GDBSTOP ();
      // These are the pointer to the functions
      // Unconstrained
      pufn     ufn;     // f
      puofg    uofg;    // f and g
      puprod   uprod;   // H*p
      punames  unames;  // UNAMES
      // Constrained
      pcofg    cofg;    // f and g
      pcfn     cfn;     // f and c
      pccfsg   ccfsg;   // c and J
      pccifg   ccifg;   // c(i) and J(i,*)'
      pcprod   cprod;   // H(L)*p
      pcnames  cnames;  // CNAMES

      // Internal function calling
      // Used to call the functions above with the minimum needed parameters
      // since most are members.
      void call_fn ();
      void call_fn_xc ();
      void call_ofg (Bool = dciTrue);
      void call_ofg (pReal, Bool = dciTrue);
      void call_ofg_xc (Bool = dciTrue);
      void call_prod (Bool, pReal, pReal);
      void call_prod_xc (Bool, pReal, pReal);
      void call_ccfsg (Bool = dciTrue, Bool = dciTrue);
      void call_ccfsg_xc (Bool = dciTrue, Bool = dciTrue);
      void call_names ();

      // Internal problem solving functions
      void analyzeJacobian ();
      void cholesky ();
      Int  normalStep ();
      Int  innerNormalDirection (Real &);
      void innerNormalPhase ();
      Int  normalSafeguard ();
      void tangentStep (Real &);
      Int  innerTangentStep (Vector &, Real &, Real &);
      void naProj (Vector &, Vector &, Vector &);
      void naProjApprox (Vector &, Vector &, Vector &);
      Int  naStep (Vector &, Vector &);
      void HiProd (Int, Real, Real, pReal, pReal, pReal, pReal, pReal);
      Int  lineSearch (const Vector &, const Vector &, const Vector &, Real &,
                Real &, Vector &);
      Int  zoom (const Vector &, const Vector &, const Vector &, Vector &, Real,
                Real, Real, Real &, Real &, Real, Real, Real, Real);
      Real interpolate (Real, Real, Real, Real, Real, Real, Bool);
      void defineParameters ();
      void readParameters ();
      void initialValues ();
      void initialization ();
      Real getTime ();
      void updateScaling_x ();
      void updateScaling_xc ();
      void scale_x (Vector &);
      void scale_xc (Vector &);
      void projectBounds_x (Vector &);
      void projectBounds_xc (Vector &);
      void update_yineq ();
      void updateMultipliers ();
      void updateMu ();
      Real calcGap ();
      Real calcYdif ();
      Real calcPen ();
      Real calcPen_xc ();
      Bool calcFeasibilityOpt ();
      void leastSquaresCG (Bool, const Vector &, Vector &, Vector &);
      void linearSystemCG (Bool, const Vector &, Vector &);

#ifndef NDEBUG
      void checkInfactibility ();
#endif
      void assert (Bool);

      // These are the values at the current point x;
      Real * f, * fxc;
      Real mu, mugap, murho;
      Vector * g;
      Sparse * H;
      Triplet * Htrip;
      Vector * c;
      Sparse * J;
      Triplet * Jtrip;
      Factor * LJ;
      Vector * gp;
      Vector * yineq;

      Vector * x, * l_bnd, * u_bnd;
      Vector * solx;
      Vector * y, * cl, * cu;
      Vector * sols;
      Bool * equatn, * linear;
      Vector * xc;
      Vector * feasOpt;

      pReal xx, l_bndx, u_bndx;
      pReal yx, clx, cux;
      pReal gx, cx;
      pReal Jx;
      pInt Jfun, Jvar;
      pReal xcx;
      pReal yineqx;

      Int * ineq_index;
      Real * scaling_matrix; // Scaling matrix
      Real * variable_scaling;

      Real Lref, Lprev, Lc, Ln, Lnew;
      Real DLV, DLH;
      Real gap, lagrgap, infacgap;

      Real normgp;
      Real normg, normc, normy, normx;
      Real normck;
      Real ngp; // normgp / (normg + 1)
      Real engp;

      // The environment variables
      Environment * env;
      Int nmax, mmax, amax;
      Int nvar, ncon, nconE, nconI, nconL, nconNL;
      std::string problemName;
      Bool start_at_one;
      Bool initialized;
      Bool running;
      Bool solved;
      Bool has_ineq, is_linear, is_bounded;
      Bool is_unlimited;
      Int display_level, debug_level, verbosity_level, table_print_level;
      Bool print_A_at_end;
      Int exit_flag;
      Real max_time, current_time, start_time;

      // Parameters and Variables
      Real c1, c2;
      Real rho, rhomax, rhomin;
      Real DeltaH, DeltaV;
      Real DeltaMax, DeltaMin, Delta0, DeltaTiny;
      Real alphaR, alphaI, alphaS;
      Real eps1, eps2, eps3;
      Real epsmu, epsgap;
      Real eta1, eta2, eta3;
      Real phi1, phi2;
      Real csic, csig;
      Real zeta1, zeta2, zeta3;
      Real beta1, beta2;
      Real kappa1, kappa2, kappa3, kappa4;
      Real thetaR, LbdMax;
      Bool GotH, first;
      Int  nSteih, nRej, nSoc, nfailv, nHprod, nRest;
      Int  tSoc, tSteih, tRej, tRest;
      Int  TangentFlag, StepFlag, NormalFlag;
      Int  iter, maxit, maxitSteih, minitSteih, relitSteih;
      Int  minstep, itssmll, maxrest, maxssmll;
      Int  cholFacs;
      Bool Aavail, gavail, LimLbd, FreshA;
      Real minBk;
      Bool use_conjugate_gradient;
      Bool scale_normal;
      Bool normal_fail_reboot;
      Bool partial_penalization, project_dcp, project_dn;
      Bool trustWorstdn, trustConvexBox, penal_trust;
      Real cholesky_correction, chol_correction_increase,
           cholesky_base_correction;
      Bool cholesky_failed;
      Real MaxDiag, MinDiag;
      Real infeasible_gradient;
      Real objective_scaling, max_objective_scaling, max_constraint_scaling,
           max_variable_scaling;
      Real *constraint_scaling;
      Bool use_objective_scaling, use_variable_scaling, use_constraint_scaling;
      Bool use_soc, use_normal_safe_guard;
      Bool use_lsmr;
      Int  objfun_count;
      Real infeasibility_tol;
      Int  cuter_status;
      Int nvarshowmax, nconshowmax;
  };
}

#endif
