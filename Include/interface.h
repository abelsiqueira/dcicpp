#ifndef dci_interface_h
#define dci_interface_h

#include "triplet.h"

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
      void cuter () { start_at_one = dciTrue; };

      // Sine quibus non
      void unc_setup (size_t, Real *, Real *, Real *);
      void con_setup (size_t, Real *, Real *, Real *, size_t, Real *, Real *,
          Bool *);
      void set_linear (size_t n, Bool *V);

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
      Int  verticalStep ();
      Int  innerVerticalStep (Real &);
      Int  verticalSafeguard ();
      void horizontalStep (Real &);
      Int  innerHorizontalStep (Vector &, Real &, Real &);
      Int  leastSquaresTrustRegion (Vector &, pReal, pReal, pReal);
      void naProj (Vector &, Vector &, Vector &);
      void naProjApprox (Vector &, Vector &, Vector &);
      Int  naStep (Vector &, Vector &);
      Int  dcibfgs (const Vector &, Int &);
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

      void copy_sx () {
        for (int i = 0; i < nconI; i++)
          sx[i] = xx[nvar+i];
      }
      void copy_scx () {
        for (int i = 0; i < nconI; i++)
          scx[i] = xcx[nvar+i];
      }

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
      Vector * s; //Slacks s = -c_I(x)
      Vector * sols;
      Bool * equatn, * linear;
      Vector * xc, *sc;
      Vector * feasOpt;

      pReal xx, l_bndx, u_bndx;
      pReal yx, clx, cux;
      pReal sx;
      pReal gx, cx;
      pReal Jx;
      pInt Jfun, Jvar;
      pReal xcx, scx;
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
      Real kappa1, kappa2, kappa3, kappa4;
      Real thetaR, LbdMax;
      Bool GotH, first;
      Int  nSteih, nRej, nSoc, nfailv, nHprod, nRest, nbfgs;
      Int  tSoc, tSteih, tRej, tbfgs, tRest;
      Int  HorzFlag, StepFlag, VertFlag;
      Int  iter, maxit, maxitSteih, minitSteih, relitSteih;
      Int  minstep, itssmll, maxrest, maxssmll;
      Int  bfgsupd;
      Int  cholFacs;
      Bool Aavail, gavail, LimLbd, FreshA;
      Real minBk;
      Bool use_conjugate_gradient;
      Bool scale_vertical;
      Bool vertical_fail_reboot;
      Bool partial_penalization, project_dcp, project_dn, project_bfgs;
      Bool trustWorstdn, trustConvexBox, penal_trust, penal_bfgs;
      Real cholesky_correction;
      Bool cholesky_failed;
      Real MaxDiag, MinDiag;
      Real infeasible_gradient;
      Real objective_scaling, max_objective_scaling, max_constraint_scaling;
      Real *constraint_scaling;
      Bool use_objective_scaling, use_variable_scaling, use_constraint_scaling;
      Bool use_soc;
      Int  objfun_count;
      Int  nfix, *fixed_index;
  };
}

#endif
